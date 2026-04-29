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

type ColumnAliasSet struct {
	TaskID   []string
	Hash     []string
	NativeID []string
	Workdir  []string
	Process  []string
	Name     []string
	Tag      []string
	Status   []string
	Exit     []string
	Duration []string
	Realtime []string
	CPU      []string
	Memory   []string
	PeakRSS  []string
	PeakVMem []string
}

type NameParts struct {
	FullName      string
	Process       string
	Tag           string
	SelectorTerms []string
}

type NormalizedRecord struct {
	RowOrder   int64
	TaskID     string
	Hash       string
	NativeID   string
	Workdir    string
	Process    string
	Name       string
	Tag        string
	Status     string
	Exit       string
	Duration   string
	Realtime   string
	CPUDisplay string
	Memory     string
	PeakRSS    string
	PeakVMem   string
	Columns    map[string]string
}

func DefaultColumnAliasSet() ColumnAliasSet {
	return ColumnAliasSet{
		TaskID:   []string{"task_id", "taskid", "task"},
		Hash:     []string{"hash"},
		NativeID: []string{"native_id", "nativeid", "native id"},
		Workdir:  []string{"workdir", "work_dir", "work-dir", "work directory"},
		Process:  []string{"process", "module"},
		Name:     []string{"name"},
		Tag:      []string{"tag"},
		Status:   []string{"status"},
		Exit:     []string{"exit", "exit_status", "exitstatus", "exit_code", "exitcode"},
		Duration: []string{"duration"},
		Realtime: []string{"realtime", "real_time", "real time", "walltime", "wall_time"},
		CPU:      []string{"cpus", "cpu", "%cpu", "pcpu"},
		Memory:   []string{"memory", "mem"},
		PeakRSS:  []string{"peak_rss", "peak-rss", "rss"},
		PeakVMem: []string{"peak_vmem", "peak-vmem", "vmem"},
	}
}

func NormalizeRecordColumns(record RawRecord, aliases ColumnAliasSet) (NormalizedRecord, error) {
	type aliasGroup struct {
		name    string
		aliases []string
		set     func(*NormalizedRecord, string)
	}

	normalizeColumnName := func(value string) string {
		return strings.ToLower(strings.TrimSpace(value))
	}

	groups := []aliasGroup{
		{name: "task_id", aliases: aliases.TaskID, set: func(normalized *NormalizedRecord, value string) { normalized.TaskID = value }},
		{name: "hash", aliases: aliases.Hash, set: func(normalized *NormalizedRecord, value string) { normalized.Hash = value }},
		{name: "native_id", aliases: aliases.NativeID, set: func(normalized *NormalizedRecord, value string) { normalized.NativeID = value }},
		{name: "workdir", aliases: aliases.Workdir, set: func(normalized *NormalizedRecord, value string) { normalized.Workdir = value }},
		{name: "process", aliases: aliases.Process, set: func(normalized *NormalizedRecord, value string) { normalized.Process = value }},
		{name: "name", aliases: aliases.Name, set: func(normalized *NormalizedRecord, value string) { normalized.Name = value }},
		{name: "tag", aliases: aliases.Tag, set: func(normalized *NormalizedRecord, value string) { normalized.Tag = value }},
		{name: "status", aliases: aliases.Status, set: func(normalized *NormalizedRecord, value string) { normalized.Status = value }},
		{name: "exit", aliases: aliases.Exit, set: func(normalized *NormalizedRecord, value string) { normalized.Exit = value }},
		{name: "duration", aliases: aliases.Duration, set: func(normalized *NormalizedRecord, value string) { normalized.Duration = value }},
		{name: "realtime", aliases: aliases.Realtime, set: func(normalized *NormalizedRecord, value string) { normalized.Realtime = value }},
		{name: "cpu", aliases: aliases.CPU, set: func(normalized *NormalizedRecord, value string) { normalized.CPUDisplay = value }},
		{name: "memory", aliases: aliases.Memory, set: func(normalized *NormalizedRecord, value string) { normalized.Memory = value }},
		{name: "peak_rss", aliases: aliases.PeakRSS, set: func(normalized *NormalizedRecord, value string) { normalized.PeakRSS = value }},
		{name: "peak_vmem", aliases: aliases.PeakVMem, set: func(normalized *NormalizedRecord, value string) { normalized.PeakVMem = value }},
	}

	visitUniqueAliases := func(group aliasGroup, visit func(alias string, key string) (bool, error)) error {
		seenInGroup := make(map[string]bool, len(group.aliases))
		for _, alias := range group.aliases {
			key := normalizeColumnName(alias)
			if key == "" {
				return fmt.Errorf("normalize record columns: blank alias for %s", group.name)
			}
			if seenInGroup[key] {
				continue
			}
			seenInGroup[key] = true
			stop, err := visit(alias, key)
			if err != nil {
				return err
			}
			if stop {
				return nil
			}
		}
		return nil
	}

	aliasOwners := make(map[string]string)
	for _, group := range groups {
		if err := visitUniqueAliases(group, func(alias string, key string) (bool, error) {
			if owner, ok := aliasOwners[key]; ok && owner != group.name {
				return false, fmt.Errorf("normalize record columns: alias %q maps to both %s and %s", alias, owner, group.name)
			}
			aliasOwners[key] = group.name
			return false, nil
		}); err != nil {
			return NormalizedRecord{}, err
		}
	}

	normalized := NormalizedRecord{
		RowOrder: record.RowOrder,
		Columns:  make(map[string]string, len(record.Columns)),
	}

	rawColumnsByAliasKey := make(map[string]string, len(record.Columns))
	for columnName, value := range record.Columns {
		normalized.Columns[columnName] = value

		key := normalizeColumnName(columnName)
		if existingColumnName, ok := rawColumnsByAliasKey[key]; ok && existingColumnName != columnName {
			if _, isKnownAlias := aliasOwners[key]; isKnownAlias {
				return NormalizedRecord{}, fmt.Errorf("normalize record columns: column names %q and %q both match alias %q", existingColumnName, columnName, key)
			}
			continue
		}
		rawColumnsByAliasKey[key] = columnName
	}

	for _, group := range groups {
		if err := visitUniqueAliases(group, func(alias string, key string) (bool, error) {
			columnName, ok := rawColumnsByAliasKey[key]
			if !ok {
				return false, nil
			}
			group.set(&normalized, record.Columns[columnName])
			return true, nil
		}); err != nil {
			return NormalizedRecord{}, err
		}
	}

	return normalized, nil
}

func DeriveNameParts(fullName string) NameParts {
	trimmed := strings.TrimSpace(fullName)
	if trimmed == "" {
		return NameParts{}
	}

	parts := NameParts{
		FullName: trimmed,
		Process:  trimmed,
	}

	if strings.HasSuffix(trimmed, ")") {
		open := strings.LastIndex(trimmed, " (")
		if open > 0 && open < len(trimmed)-1 {
			process := strings.TrimSpace(trimmed[:open])
			tag := strings.TrimSpace(trimmed[open+2 : len(trimmed)-1])
			if process != "" && tag != "" && !strings.ContainsAny(tag, "()") {
				parts.Process = process
				parts.Tag = tag
			}
		}
	}

	seen := make(map[string]bool, 3)
	addTerm := func(term string) {
		term = strings.TrimSpace(term)
		if term == "" || seen[term] {
			return
		}
		seen[term] = true
		parts.SelectorTerms = append(parts.SelectorTerms, term)
	}

	addTerm(parts.FullName)
	addTerm(parts.Process)
	addTerm(parts.Tag)

	return parts
}

func TaskFromNormalizedRecord(runDir domain.RunDir, record NormalizedRecord) (domain.Task, error) {
	candidates := []struct {
		name  string
		value string
	}{
		{name: "hash", value: record.Hash},
		{name: "workdir", value: record.Workdir},
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
		return domain.Task{}, fmt.Errorf("task from normalized record row %d: no parseable hash or workdir (%s)", record.RowOrder, detail)
	}

	workdirInput := record.Workdir
	if strings.TrimSpace(workdirInput) == "" {
		workdirInput = record.Hash
	}
	workdir, err := ResolveTaskWorkdir(runDir, canonicalID, workdirInput)
	if err != nil {
		return domain.Task{}, fmt.Errorf("task from normalized record row %d: resolve workdir for task %q: %w", record.RowOrder, canonicalID, err)
	}

	exit, err := ParseNullableExit(record.Exit)
	if err != nil {
		return domain.Task{}, fmt.Errorf("task from normalized record row %d: parse exit for task %q: %w", record.RowOrder, canonicalID, err)
	}

	nameParts := DeriveNameParts(record.Name)
	process := record.Process
	if strings.TrimSpace(process) == "" {
		process = nameParts.Process
	}
	tag := record.Tag
	if strings.TrimSpace(tag) == "" {
		tag = nameParts.Tag
	}

	memoryDisplay := record.Memory
	if strings.TrimSpace(memoryDisplay) == "" {
		peakDisplays := make([]string, 0, 2)
		if peakRSS := strings.TrimSpace(record.PeakRSS); peakRSS != "" {
			peakDisplays = append(peakDisplays, "peak_rss="+peakRSS)
		}
		if peakVMem := strings.TrimSpace(record.PeakVMem); peakVMem != "" {
			peakDisplays = append(peakDisplays, "peak_vmem="+peakVMem)
		}
		memoryDisplay = strings.Join(peakDisplays, "; ")
	}

	return domain.Task{
		RowOrder: record.RowOrder,
		ID:       canonicalID,
		Status:   NormalizeTaskStatus(record.Status),
		Process:  process,
		Name:     record.Name,
		Tag:      tag,
		Workdir:  workdir,
		Exit:     exit,
		Duration: record.Duration,
		Realtime: record.Realtime,
		CPUs:     record.CPUDisplay,
		Memory:   memoryDisplay,
	}, nil
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
	normalized, err := NormalizeRecordColumns(record, DefaultColumnAliasSet())
	if err != nil {
		return domain.Task{}, fmt.Errorf("normalize trace record row %d: %w", record.RowOrder, err)
	}

	task, err := TaskFromNormalizedRecord(runDir, normalized)
	if err != nil {
		return domain.Task{}, fmt.Errorf("normalize trace record row %d: %w", record.RowOrder, err)
	}
	if strings.TrimSpace(normalized.Process) != "" && strings.TrimSpace(normalized.Tag) == "" {
		task.Tag = ""
	}

	return task, nil
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

func resolveHashPrefixWorkdir(runDir domain.RunDir, canonicalID string) (string, error) {
	canonical, ok := NormalizeCanonicalTaskID(canonicalID)
	if !ok {
		return "", fmt.Errorf("resolve hash prefix workdir: invalid canonical id %q", canonicalID)
	}
	if strings.TrimSpace(runDir.Path) == "" {
		return "", nil
	}

	idParts := strings.Split(canonical, "/")
	shard := filepath.Join(runDir.Path, "work", idParts[0])
	entries, err := os.ReadDir(shard)
	if err != nil {
		if os.IsNotExist(err) {
			return "", nil
		}
		return "", fmt.Errorf("resolve hash prefix workdir: read shard %q for canonical id %q: %w", shard, canonical, err)
	}

	matches := make([]string, 0, 1)
	remainingPrefix := idParts[1]
	for _, entry := range entries {
		if !strings.HasPrefix(strings.ToLower(entry.Name()), remainingPrefix) {
			continue
		}

		info, err := entry.Info()
		if err != nil {
			return "", fmt.Errorf("resolve hash prefix workdir: stat entry %q in shard %q for canonical id %q: %w", entry.Name(), shard, canonical, err)
		}
		if !info.IsDir() {
			continue
		}

		matches = append(matches, filepath.Join(shard, entry.Name()))
	}

	if len(matches) == 0 {
		return "", nil
	}
	if len(matches) > 1 {
		return "", fmt.Errorf("resolve hash prefix workdir: ambiguous canonical id %q in shard %q: %d matching directories (%s)", canonical, shard, len(matches), strings.Join(matches, ", "))
	}

	return filepath.Clean(matches[0]), nil
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

	workdir, err := resolveHashPrefixWorkdir(runDir, canonical)
	if err != nil {
		return "", fmt.Errorf("resolve task workdir: %w", err)
	}
	return workdir, nil
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
