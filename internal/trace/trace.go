package trace

import (
	"context"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/mskilab-org/gosh/internal/domain"
)

type Delimiter rune

const (
	DelimiterComma Delimiter = ','
	DelimiterTab   Delimiter = '\t'
)

type RawRecord struct {
	RowOrder int64
	Columns  map[string]string
}

func isHexString(value string) bool {
	if value == "" {
		return false
	}
	for _, r := range value {
		if !((r >= '0' && r <= '9') || (r >= 'a' && r <= 'f') || (r >= 'A' && r <= 'F')) {
			return false
		}
	}
	return true
}

func canonicalTaskIDFromParts(prefix string, rest string) (string, bool) {
	if len(prefix) != 2 || !isHexString(prefix) || !isHexString(rest) {
		return "", false
	}
	return strings.ToLower(prefix) + "/" + strings.ToLower(rest), true
}

func NormalizeCanonicalTaskID(value string) (string, bool) {
	value = strings.TrimSpace(value)
	parts := strings.Split(value, "/")
	if len(parts) != 2 {
		return "", false
	}
	return canonicalTaskIDFromParts(parts[0], parts[1])
}

func stripTaskIDEnclosure(value string) string {
	value = strings.TrimSpace(value)
	for {
		if len(value) < 2 {
			return value
		}

		first := value[0]
		last := value[len(value)-1]
		if (first == '"' && last == '"') ||
			(first == '\'' && last == '\'') ||
			(first == '[' && last == ']') ||
			(first == '(' && last == ')') ||
			(first == '{' && last == '}') {
			value = strings.TrimSpace(value[1 : len(value)-1])
			continue
		}

		return value
	}
}

func DetectDelimiter(path string) (Delimiter, error) {
	if path == "" {
		return 0, fmt.Errorf("detect trace delimiter: empty path (supported extensions: .csv, .tsv, .txt)")
	}

	if len(path) >= 4 {
		switch path[len(path)-4:] {
		case ".csv":
			return DelimiterComma, nil
		case ".tsv", ".txt":
			return DelimiterTab, nil
		}
	}

	extension := ""
	for i := len(path) - 1; i >= 0; i-- {
		switch path[i] {
		case '.':
			extension = path[i:]
			i = -1
		case '/', '\\':
			i = -1
		}
	}
	if extension == "" {
		extension = "<none>"
	}

	return 0, fmt.Errorf("detect trace delimiter: unsupported extension %q for %q (supported extensions: .csv, .tsv, .txt)", extension, path)
}

func ParseTrace(ctx context.Context, runDir domain.RunDir, source domain.SourceFingerprint) ([]domain.Task, error) {
	if ctx == nil {
		return nil, fmt.Errorf("parse trace: nil context")
	}
	if err := ctx.Err(); err != nil {
		return nil, fmt.Errorf("parse trace: %w", err)
	}
	if source.Kind != domain.SourceKindTrace {
		return nil, fmt.Errorf("parse trace: invalid trace source kind %q (want %q)", source.Kind, domain.SourceKindTrace)
	}

	delimiter, err := DetectDelimiter(source.Path)
	if err != nil {
		return nil, fmt.Errorf("parse trace: %w", err)
	}

	file, err := os.Open(source.Path)
	if err != nil {
		return nil, fmt.Errorf("parse trace: open source %q: %w", source.Path, err)
	}
	defer file.Close()

	records, err := ParseTraceRecords(file, delimiter)
	if err != nil {
		return nil, fmt.Errorf("parse trace: %w", err)
	}

	tasks := make([]domain.Task, 0, len(records))
	for _, record := range records {
		if err := ctx.Err(); err != nil {
			return nil, fmt.Errorf("parse trace: %w", err)
		}

		task, err := NormalizeTraceRecord(runDir, record)
		if err != nil {
			return nil, fmt.Errorf("parse trace: %w", err)
		}
		tasks = append(tasks, task)
	}

	return tasks, nil
}

func ParseTraceRecords(reader io.Reader, delimiter Delimiter) ([]RawRecord, error) {
	if reader == nil {
		return nil, fmt.Errorf("parse trace records: nil reader")
	}
	if delimiter != DelimiterComma && delimiter != DelimiterTab {
		return nil, fmt.Errorf("parse trace records: unsupported delimiter %q (supported delimiters: comma %q or tab)", rune(delimiter), rune(DelimiterComma))
	}

	csvReader := csv.NewReader(reader)
	csvReader.Comma = rune(delimiter)

	headers, err := csvReader.Read()
	if err == io.EOF {
		return nil, fmt.Errorf("parse trace records: missing header row")
	}
	if err != nil {
		return nil, fmt.Errorf("parse trace records: malformed header: %w", err)
	}
	if len(headers) == 0 {
		return nil, fmt.Errorf("parse trace records: malformed header: empty header row")
	}

	seenHeaders := make(map[string]int, len(headers))
	for i, header := range headers {
		column := i + 1
		if header == "" {
			return nil, fmt.Errorf("parse trace records: malformed header: empty column name at column %d", column)
		}
		if firstColumn, ok := seenHeaders[header]; ok {
			return nil, fmt.Errorf("parse trace records: malformed header: duplicate column name %q at column %d (first seen at column %d)", header, column, firstColumn)
		}
		seenHeaders[header] = column
	}

	records := make([]RawRecord, 0)
	for {
		values, err := csvReader.Read()
		if err == io.EOF {
			break
		}

		rowOrder := int64(len(records) + 1)
		if err != nil {
			return nil, fmt.Errorf("parse trace records: malformed row %d: %w", rowOrder, err)
		}

		columns := make(map[string]string, len(headers))
		for i, header := range headers {
			columns[header] = values[i]
		}
		records = append(records, RawRecord{RowOrder: rowOrder, Columns: columns})
	}

	return records, nil
}

func NormalizeTraceRecord(runDir domain.RunDir, record RawRecord) (domain.Task, error) {
	column := func(name string) string {
		if record.Columns == nil {
			return ""
		}
		return record.Columns[name]
	}

	hash := column("hash")
	rawWorkdir := column("workdir")

	candidates := []struct {
		name  string
		value string
	}{
		{name: "hash", value: hash},
		{name: "workdir", value: rawWorkdir},
	}

	canonicalID := ""
	parseErrors := make([]string, 0, len(candidates))
	for _, candidate := range candidates {
		if strings.TrimSpace(candidate.value) == "" {
			continue
		}

		id, err := DeriveCanonicalTaskID(candidate.value)
		if err == nil {
			canonicalID = id
			break
		}
		parseErrors = append(parseErrors, fmt.Sprintf("%s %q: %v", candidate.name, candidate.value, err))
	}

	if canonicalID == "" {
		detail := "missing hash and workdir columns"
		if len(parseErrors) > 0 {
			detail = strings.Join(parseErrors, "; ")
		}
		return domain.Task{}, fmt.Errorf("normalize trace record row %d: no parseable hash or workdir (%s)", record.RowOrder, detail)
	}

	workdirInput := rawWorkdir
	if strings.TrimSpace(workdirInput) == "" {
		workdirInput = hash
	}
	workdir, err := ResolveTaskWorkdir(runDir, canonicalID, workdirInput)
	if err != nil {
		return domain.Task{}, fmt.Errorf("normalize trace record row %d: resolve workdir for task %q: %w", record.RowOrder, canonicalID, err)
	}

	exit, err := ParseNullableExit(column("exit"))
	if err != nil {
		return domain.Task{}, fmt.Errorf("normalize trace record row %d: parse exit for task %q: %w", record.RowOrder, canonicalID, err)
	}

	return domain.Task{
		RowOrder: record.RowOrder,
		ID:       canonicalID,
		Status:   NormalizeTaskStatus(column("status")),
		Process:  column("process"),
		Name:     column("name"),
		Tag:      column("tag"),
		Workdir:  workdir,
		Exit:     exit,
		Duration: column("duration"),
		Realtime: column("realtime"),
		CPUs:     column("cpus"),
		Memory:   column("memory"),
	}, nil
}

func NormalizeTaskStatus(raw string) domain.TaskStatus {
	switch strings.ToUpper(strings.TrimSpace(raw)) {
	case string(domain.TaskStatusFailed):
		return domain.TaskStatusFailed
	case string(domain.TaskStatusCompleted):
		return domain.TaskStatusCompleted
	case string(domain.TaskStatusCached):
		return domain.TaskStatusCached
	case string(domain.TaskStatusAborted):
		return domain.TaskStatusAborted
	case string(domain.TaskStatusSubmitted):
		return domain.TaskStatusSubmitted
	case string(domain.TaskStatusRunning):
		return domain.TaskStatusRunning
	default:
		return domain.TaskStatusUnknown
	}
}

func DeriveCanonicalTaskID(rawWorkdir string) (string, error) {
	raw := strings.TrimSpace(rawWorkdir)
	if raw == "" {
		return "", fmt.Errorf("derive canonical task id: empty workdir/hash value")
	}

	value := stripTaskIDEnclosure(raw)
	if value == "" {
		return "", fmt.Errorf("derive canonical task id: empty workdir/hash value")
	}

	normalized := strings.TrimRight(strings.ReplaceAll(value, "\\", "/"), "/")
	if normalized == "" {
		return "", fmt.Errorf("derive canonical task id: empty workdir/hash value")
	}

	parts := strings.Split(normalized, "/")
	foundWorkSegment := false
	for i, part := range parts {
		if !strings.EqualFold(part, "work") {
			continue
		}
		foundWorkSegment = true
		if i+2 >= len(parts) {
			continue
		}
		if id, ok := canonicalTaskIDFromParts(parts[i+1], parts[i+2]); ok {
			return id, nil
		}
	}
	if foundWorkSegment {
		return "", fmt.Errorf("derive canonical task id: workdir value %q contains work segment but no valid xx/rest hash after it", rawWorkdir)
	}

	if len(parts) == 2 {
		if id, ok := canonicalTaskIDFromParts(parts[0], parts[1]); ok {
			return id, nil
		}
		return "", fmt.Errorf("derive canonical task id: workdir value %q is not a valid xx/rest hash", rawWorkdir)
	}

	if len(parts) == 1 && len(normalized) >= 3 && isHexString(normalized) {
		return strings.ToLower(normalized[:2]) + "/" + strings.ToLower(normalized[2:]), nil
	}

	return "", fmt.Errorf("derive canonical task id: unable to parse workdir/hash value %q as full workdir path, direct xx/rest id, or unsplit hash", rawWorkdir)
}

func ResolveTaskWorkdir(runDir domain.RunDir, canonicalID string, rawWorkdir string) (string, error) {
	cleanPath := func(value string) string {
		return filepath.Clean(strings.ReplaceAll(value, "\\", string(filepath.Separator)))
	}

	joinSlashParts := func(parts []string) string {
		return cleanPath(strings.Join(parts, "/"))
	}

	looksHashOnly := func(value string) bool {
		normalized := strings.TrimRight(strings.ReplaceAll(value, "\\", "/"), "/")
		if normalized == "" {
			return false
		}

		parts := strings.Split(normalized, "/")
		if len(parts) == 1 {
			return len(normalized) >= 3 && isHexString(normalized)
		}
		return len(parts) == 2 && len(parts[0]) == 2 && isHexString(parts[0]) && isHexString(parts[1])
	}

	looksFullPathLike := func(value string) bool {
		normalized := strings.TrimRight(strings.ReplaceAll(value, "\\", "/"), "/")
		if strings.HasPrefix(normalized, "/") {
			return true
		}
		return len(strings.Split(normalized, "/")) > 2
	}

	isMissingWorkdirValue := func(value string) bool {
		switch strings.ToLower(strings.TrimSpace(value)) {
		case "", "-", "na", "n/a", "null":
			return true
		default:
			return false
		}
	}

	extractKnownWorkdirPath := func(value string, id string) string {
		normalized := strings.TrimRight(strings.ReplaceAll(value, "\\", "/"), "/")
		parts := strings.Split(normalized, "/")

		for i, part := range parts {
			if !strings.EqualFold(part, "work") || i+2 >= len(parts) {
				continue
			}
			if len(parts[i+1]) == 2 && isHexString(parts[i+1]) && isHexString(parts[i+2]) {
				return joinSlashParts(parts[:i+3])
			}
		}

		if id == "" {
			return ""
		}
		idParts := strings.Split(id, "/")
		if len(idParts) != 2 {
			return ""
		}
		for i := 0; i+1 < len(parts); i++ {
			if strings.EqualFold(parts[i], idParts[0]) && strings.EqualFold(parts[i+1], idParts[1]) {
				return joinSlashParts(parts[:i+2])
			}
		}

		return ""
	}

	normalizeID := func(value string) (string, error) {
		value = strings.TrimSpace(value)
		if value == "" {
			return "", nil
		}
		id, err := DeriveCanonicalTaskID(value)
		if err != nil {
			return "", err
		}
		return id, nil
	}

	canonical, canonicalErr := normalizeID(canonicalID)
	raw := stripTaskIDEnclosure(rawWorkdir)
	if isMissingWorkdirValue(raw) {
		raw = ""
	}

	if raw != "" && !looksHashOnly(raw) {
		if !looksFullPathLike(raw) {
			return "", fmt.Errorf("resolve task workdir: raw workdir value %q is neither a path nor a hash", rawWorkdir)
		}
		if workdir := extractKnownWorkdirPath(raw, canonical); workdir != "" {
			return workdir, nil
		}
		return cleanPath(raw), nil
	}

	if canonicalErr != nil {
		return "", fmt.Errorf("resolve task workdir: invalid canonical id %q: %w", canonicalID, canonicalErr)
	}

	if raw != "" {
		rawID, err := DeriveCanonicalTaskID(raw)
		if err != nil {
			return "", fmt.Errorf("resolve task workdir: invalid hash value %q: %w", rawWorkdir, err)
		}
		if canonical != "" && rawID != canonical {
			return "", fmt.Errorf("resolve task workdir: hash value %q resolves to %q, not canonical id %q", rawWorkdir, rawID, canonical)
		}
		canonical = rawID
	}

	if canonical == "" || strings.TrimSpace(runDir.Path) == "" {
		return "", nil
	}

	idParts := strings.Split(canonical, "/")
	if len(idParts) != 2 {
		return "", fmt.Errorf("resolve task workdir: invalid canonical id %q", canonicalID)
	}

	candidate := filepath.Join(runDir.Path, "work", idParts[0], idParts[1])
	info, err := os.Stat(candidate)
	if err != nil {
		if os.IsNotExist(err) {
			return "", nil
		}
		return "", fmt.Errorf("resolve task workdir: stat derived workdir %q: %w", candidate, err)
	}
	if !info.IsDir() {
		return "", fmt.Errorf("resolve task workdir: derived workdir %q is not a directory", candidate)
	}

	return filepath.Clean(candidate), nil
}

func ParseNullableExit(raw string) (*int, error) {
	value := strings.TrimSpace(raw)
	switch strings.ToLower(value) {
	case "", "-", "na", "n/a", "null", "none":
		return nil, nil
	}

	parsed, err := strconv.Atoi(value)
	if err != nil {
		return nil, fmt.Errorf("parse nullable exit %q: %w", raw, err)
	}
	return &parsed, nil
}
