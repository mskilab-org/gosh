package nflog

import (
	"context"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/mskilab-org/gosh/internal/domain"
	"github.com/mskilab-org/gosh/internal/trace"
)

type FailureBlock struct {
	Process string
	Name    string
	Workdir string
	Exit    *int
	Block   string
}

func normalizeLogText(value string) string {
	value = strings.ReplaceAll(value, "\r\n", "\n")
	return strings.ReplaceAll(value, "\r", "\n")
}

func cleanSimpleEvidenceValue(value string) string {
	value = strings.TrimSpace(value)
	value = strings.Trim(value, "`\"'")
	value = strings.TrimRight(value, ",;")
	return strings.TrimSpace(value)
}

func firstNonBlankTrimmedLine(lines []string, start int) string {
	for i := start; i < len(lines); i++ {
		candidate := strings.TrimSpace(lines[i])
		if candidate != "" {
			return candidate
		}
	}
	return ""
}

func firstHeaderValue(lines []string, pattern *regexp.Regexp) (string, bool) {
	for i, line := range lines {
		matches := pattern.FindStringSubmatch(line)
		if len(matches) == 0 {
			continue
		}

		value := strings.TrimSpace(matches[1])
		if value == "" {
			value = firstNonBlankTrimmedLine(lines, i+1)
		}
		return value, true
	}
	return "", false
}

func readNormalizedLogText(reader io.Reader, context string) (string, error) {
	if reader == nil {
		return "", fmt.Errorf("%s: nil reader", context)
	}

	data, err := io.ReadAll(reader)
	if err != nil {
		return "", fmt.Errorf("%s: read log: %w", context, err)
	}

	return normalizeLogText(string(data)), nil
}

func selectedLogPath(artifacts domain.ArtifactSet) string {
	if artifacts.Log != nil && artifacts.Log.Path != "" {
		return artifacts.Log.Path
	}
	return "selected Nextflow log"
}

func searchedPatternsText(artifacts domain.ArtifactSet) string {
	if len(artifacts.SearchedPatterns) == 0 {
		return ""
	}
	return strings.Join(artifacts.SearchedPatterns, ", ")
}

func splitProcessLabel(label string) (string, string) {
	label = strings.TrimSpace(label)
	if label == "" {
		return "", ""
	}

	if strings.HasSuffix(label, ")") {
		open := strings.LastIndex(label, " (")
		if open > 0 && open < len(label)-1 {
			process := strings.TrimSpace(label[:open])
			name := strings.TrimSpace(label[open+2 : len(label)-1])
			if process != "" && name != "" && !strings.ContainsAny(name, "()") {
				return process, name
			}
		}
	}

	return label, ""
}

func ParseLogOnlyFailures(ctx context.Context, runDir domain.RunDir, source domain.SourceFingerprint) ([]domain.LogOnlyFailure, error) {
	if ctx == nil {
		return nil, fmt.Errorf("parse log-only failures: nil context")
	}
	if err := ctx.Err(); err != nil {
		return nil, fmt.Errorf("parse log-only failures: %w", err)
	}
	if source.Kind != domain.SourceKindLog {
		return nil, fmt.Errorf("parse log-only failures: invalid log source kind %q (want %q)", source.Kind, domain.SourceKindLog)
	}

	file, err := os.Open(source.Path)
	if err != nil {
		return nil, fmt.Errorf("parse log-only failures: open source %q: %w", source.Path, err)
	}
	defer file.Close()

	blocks, err := ExtractFailureBlocks(file)
	if err != nil {
		return nil, fmt.Errorf("parse log-only failures: %w", err)
	}

	failures := make([]domain.LogOnlyFailure, 0, len(blocks))
	for _, block := range blocks {
		if err := ctx.Err(); err != nil {
			return nil, fmt.Errorf("parse log-only failures: %w", err)
		}

		failure, err := NormalizeFailureBlock(runDir, block)
		if err != nil {
			return nil, fmt.Errorf("parse log-only failures: %w", err)
		}
		failures = append(failures, failure)
	}

	return failures, nil
}

func ParseLogOnlyTaskEvidence(ctx context.Context, runDir domain.RunDir, source domain.SourceFingerprint) ([]domain.LogOnlyTaskEvidence, error) {
	const workdirEvidenceMaxBytes int64 = 4096

	if ctx == nil {
		return nil, fmt.Errorf("parse log-only task evidence: nil context")
	}
	if err := ctx.Err(); err != nil {
		return nil, fmt.Errorf("parse log-only task evidence: %w", err)
	}
	if source.Kind != domain.SourceKindLog {
		return nil, fmt.Errorf("parse log-only task evidence: invalid log source kind %q (want %q)", source.Kind, domain.SourceKindLog)
	}

	file, err := os.Open(source.Path)
	if err != nil {
		return nil, fmt.Errorf("parse log-only task evidence: open source %q: %w", source.Path, err)
	}
	defer file.Close()

	data, err := io.ReadAll(file)
	if err != nil {
		return nil, fmt.Errorf("parse log-only task evidence: read source %q: %w", source.Path, err)
	}
	text := string(data)

	lifecycleEvidence, err := ExtractLifecycleEvidence(strings.NewReader(text))
	if err != nil {
		return nil, fmt.Errorf("parse log-only task evidence: %w", err)
	}
	if err := ctx.Err(); err != nil {
		return nil, fmt.Errorf("parse log-only task evidence: %w", err)
	}

	failureBlocks, err := ExtractFailureBlocks(strings.NewReader(text))
	if err != nil {
		return nil, fmt.Errorf("parse log-only task evidence: %w", err)
	}

	copyExit := func(exit *int) *int {
		if exit == nil {
			return nil
		}
		value := *exit
		return &value
	}

	copySources := func(sources []domain.LogOnlyEvidenceSource) []domain.LogOnlyEvidenceSource {
		if len(sources) == 0 {
			return nil
		}
		return append([]domain.LogOnlyEvidenceSource(nil), sources...)
	}

	markLogSourcePath := func(item domain.LogOnlyTaskEvidence) domain.LogOnlyTaskEvidence {
		for index := range item.Sources {
			if item.Sources[index].Path == "" {
				item.Sources[index].Path = source.Path
			}
		}
		return item
	}

	resolveReferencedWorkdir := func(item domain.LogOnlyTaskEvidence) (domain.LogOnlyTaskEvidence, error) {
		if strings.TrimSpace(item.Workdir) != "" || strings.TrimSpace(item.ID) == "" {
			return item, nil
		}

		workdir, err := trace.ResolveTaskWorkdir(runDir, item.ID, item.ID)
		if err != nil {
			return item, err
		}
		item.Workdir = workdir
		return item, nil
	}

	mergeEvidence := func(existing *domain.LogOnlyTaskEvidence, next domain.LogOnlyTaskEvidence) {
		if existing.ID == "" && next.ID != "" {
			existing.ID = next.ID
		}
		if next.Workdir != "" {
			existing.Workdir = next.Workdir
		}
		if next.Process != "" {
			existing.Process = next.Process
		}
		if next.Name != "" {
			existing.Name = next.Name
		}
		if next.ObservedStatus != "" {
			existing.ObservedStatus = next.ObservedStatus
		}
		if next.Exit != nil {
			existing.Exit = copyExit(next.Exit)
		}
		if next.ErrorSummary != "" {
			existing.ErrorSummary = next.ErrorSummary
		}
		if next.ErrorBlock != "" {
			existing.ErrorBlock = next.ErrorBlock
		}
		if len(next.Sources) > 0 {
			existing.Sources = append(existing.Sources, copySources(next.Sources)...)
		}
		if next.Completeness != "" {
			existing.Completeness = next.Completeness
		} else if existing.Completeness == "" {
			existing.Completeness = domain.LogOnlyEvidencePartial
		}
		if next.CommandFilesAvailable {
			existing.CommandFilesAvailable = true
		}
	}

	evidence := make([]domain.LogOnlyTaskEvidence, 0, len(lifecycleEvidence)+len(failureBlocks))
	byID := make(map[string]int)
	appendOrMerge := func(item domain.LogOnlyTaskEvidence) {
		if item.Completeness == "" {
			item.Completeness = domain.LogOnlyEvidencePartial
		}
		item.Sources = copySources(item.Sources)
		item.Exit = copyExit(item.Exit)

		id := strings.TrimSpace(item.ID)
		if id != "" {
			item.ID = id
			if index, ok := byID[id]; ok {
				mergeEvidence(&evidence[index], item)
				return
			}
			byID[id] = len(evidence)
		}

		evidence = append(evidence, item)
	}

	for _, item := range lifecycleEvidence {
		if err := ctx.Err(); err != nil {
			return nil, fmt.Errorf("parse log-only task evidence: %w", err)
		}

		item = markLogSourcePath(item)
		item, err = resolveReferencedWorkdir(item)
		if err != nil {
			return nil, fmt.Errorf("parse log-only task evidence: resolve lifecycle workdir %q: %w", item.ID, err)
		}
		appendOrMerge(item)
	}

	for _, block := range failureBlocks {
		if err := ctx.Err(); err != nil {
			return nil, fmt.Errorf("parse log-only task evidence: %w", err)
		}

		failure, err := NormalizeFailureBlock(runDir, block)
		if err != nil {
			return nil, fmt.Errorf("parse log-only task evidence: %w", err)
		}
		appendOrMerge(LogOnlyEvidenceFromFailure(failure, source))
	}

	enriched, err := EnrichLogOnlyEvidenceFromWorkdirs(ctx, evidence, workdirEvidenceMaxBytes)
	if err != nil {
		return nil, fmt.Errorf("parse log-only task evidence: %w", err)
	}
	return enriched, nil
}

func ExtractLifecycleEvidence(reader io.Reader) ([]domain.LogOnlyTaskEvidence, error) {
	text, err := readNormalizedLogText(reader, "extract lifecycle evidence")
	if err != nil {
		return nil, err
	}
	if strings.TrimSpace(text) == "" {
		return []domain.LogOnlyTaskEvidence{}, nil
	}

	lifecycleLine := regexp.MustCompile(`(?i)\[([0-9a-f]{2}/[0-9a-f]+|[0-9a-f]{3,})\]\s+((?:re[-\s]*)?submitted|submit|cached|completed|complete|failed|aborted|started|running|launch(?:ed|ing)?)\s+process\s*>\s*(.+?)\s*$`)

	statusFromEvent := func(event string) domain.TaskStatus {
		normalized := strings.ToLower(strings.TrimSpace(event))
		normalized = strings.ReplaceAll(normalized, "-", "")
		normalized = strings.ReplaceAll(normalized, " ", "")

		switch {
		case strings.Contains(normalized, "cache"):
			return domain.TaskStatusCached
		case strings.Contains(normalized, "complete"):
			return domain.TaskStatusCompleted
		case strings.Contains(normalized, "fail"):
			return domain.TaskStatusFailed
		case strings.Contains(normalized, "abort"):
			return domain.TaskStatusAborted
		case strings.Contains(normalized, "submit"):
			return domain.TaskStatusSubmitted
		case strings.Contains(normalized, "start") || strings.Contains(normalized, "run") || strings.Contains(normalized, "launch"):
			return domain.TaskStatusRunning
		default:
			return domain.TaskStatusUnknown
		}
	}

	evidence := make([]domain.LogOnlyTaskEvidence, 0)
	seen := make(map[string]int)
	lines := strings.Split(text, "\n")
	for i, line := range lines {
		matches := lifecycleLine.FindStringSubmatch(line)
		if len(matches) != 4 {
			continue
		}

		id, err := trace.DeriveCanonicalTaskID(matches[1])
		if err != nil {
			continue
		}

		status := statusFromEvent(matches[2])
		if status == domain.TaskStatusUnknown {
			continue
		}

		process, name := splitProcessLabel(cleanSimpleEvidenceValue(matches[3]))
		source := domain.LogOnlyEvidenceSource{
			Kind:   domain.LogOnlyEvidenceSourceLog,
			Detail: fmt.Sprintf("lifecycle line %d: %s process", i+1, strings.ToLower(strings.TrimSpace(matches[2]))),
		}

		if index, ok := seen[id]; ok {
			if process != "" {
				evidence[index].Process = process
			}
			if name != "" {
				evidence[index].Name = name
			}
			evidence[index].ObservedStatus = status
			evidence[index].Sources = append(evidence[index].Sources, source)
			continue
		}

		seen[id] = len(evidence)
		evidence = append(evidence, domain.LogOnlyTaskEvidence{
			ID:             id,
			Process:        process,
			Name:           name,
			ObservedStatus: status,
			Sources:        []domain.LogOnlyEvidenceSource{source},
			Completeness:   domain.LogOnlyEvidencePartial,
		})
	}

	return evidence, nil
}

func LogOnlyEvidenceFromFailure(failure domain.LogOnlyFailure, source domain.SourceFingerprint) domain.LogOnlyTaskEvidence {
	return domain.LogOnlyTaskEvidence{
		ID:             failure.ID,
		Workdir:        failure.Workdir,
		Process:        failure.Process,
		Name:           failure.Name,
		ObservedStatus: domain.TaskStatusFailed,
		Exit:           failure.Exit,
		ErrorSummary:   failure.ErrorSummary,
		ErrorBlock:     failure.ErrorBlock,
		Sources: []domain.LogOnlyEvidenceSource{
			{
				Kind:   domain.LogOnlyEvidenceSourceLog,
				Path:   source.Path,
				Detail: "failure block parsed from selected log",
			},
		},
		Completeness: domain.LogOnlyEvidencePartial,
	}
}

func EnrichLogOnlyEvidenceFromWorkdirs(ctx context.Context, evidence []domain.LogOnlyTaskEvidence, maxBytes int64) ([]domain.LogOnlyTaskEvidence, error) {
	if ctx == nil {
		return nil, fmt.Errorf("enrich log-only evidence from workdirs: nil context")
	}
	if err := ctx.Err(); err != nil {
		return nil, fmt.Errorf("enrich log-only evidence from workdirs: %w", err)
	}
	if maxBytes <= 0 {
		return nil, fmt.Errorf("enrich log-only evidence from workdirs: max bytes must be positive")
	}

	type workdirEvidence struct {
		exit                  *int
		errorSummary          string
		errorBlock            string
		sources               []domain.LogOnlyEvidenceSource
		commandFilesAvailable bool
	}

	maxBytesAsInt := func() int {
		maxInt := int(^uint(0) >> 1)
		if maxBytes > int64(maxInt) {
			return maxInt
		}
		return int(maxBytes)
	}

	statRegularFile := func(path string) (bool, error) {
		if err := ctx.Err(); err != nil {
			return false, err
		}

		info, err := os.Stat(path)
		if err != nil {
			if os.IsNotExist(err) {
				return false, nil
			}
			return false, fmt.Errorf("stat %q: %w", path, err)
		}
		if info.IsDir() {
			return false, fmt.Errorf("%q is a directory", path)
		}
		return true, nil
	}

	readBoundedText := func(path string) (string, bool, error) {
		if err := ctx.Err(); err != nil {
			return "", false, err
		}

		file, err := os.Open(path)
		if err != nil {
			if os.IsNotExist(err) {
				return "", false, nil
			}
			return "", false, fmt.Errorf("open %q: %w", path, err)
		}
		defer file.Close()

		limit := maxBytes + 1
		if limit <= maxBytes {
			limit = maxBytes
		}
		data, err := io.ReadAll(io.LimitReader(file, limit))
		if err != nil {
			return "", false, fmt.Errorf("read %q: %w", path, err)
		}

		truncated := int64(len(data)) > maxBytes
		if truncated {
			data = data[:len(data)-1]
		}
		return normalizeLogText(string(data)), truncated, nil
	}

	parseExitCode := func(path string, content string) (*int, error) {
		value := cleanSimpleEvidenceValue(content)
		if value == "" {
			return nil, fmt.Errorf("read .exitcode %q: invalid exit value: missing value", path)
		}

		parsed, err := strconv.Atoi(value)
		if err != nil {
			return nil, fmt.Errorf("read .exitcode %q: invalid exit value %q: %w", path, value, err)
		}
		return &parsed, nil
	}

	enrichWorkdir := func(workdir string) (workdirEvidence, error) {
		result := workdirEvidence{}
		if err := ctx.Err(); err != nil {
			return result, err
		}

		cleanWorkdir := filepath.Clean(workdir)
		info, err := os.Stat(cleanWorkdir)
		if err != nil {
			if os.IsNotExist(err) {
				return result, nil
			}
			return result, fmt.Errorf("stat referenced workdir %q: %w", cleanWorkdir, err)
		}
		if !info.IsDir() {
			return result, fmt.Errorf("referenced workdir %q is not a directory", cleanWorkdir)
		}

		result.sources = append(result.sources, domain.LogOnlyEvidenceSource{
			Kind:   domain.LogOnlyEvidenceSourceWorkdir,
			Path:   cleanWorkdir,
			Detail: "referenced workdir checked for .exitcode, .command.err, and .command.log evidence",
		})

		exitPath := filepath.Join(cleanWorkdir, ".exitcode")
		exitExists, err := statRegularFile(exitPath)
		if err != nil {
			return result, fmt.Errorf("enrich referenced workdir %q: %w", cleanWorkdir, err)
		}
		if exitExists {
			content, truncated, err := readBoundedText(exitPath)
			if err != nil {
				return result, fmt.Errorf("enrich referenced workdir %q: %w", cleanWorkdir, err)
			}
			if truncated {
				return result, fmt.Errorf("read .exitcode %q: invalid exit value: file exceeds max bytes %d", exitPath, maxBytes)
			}

			exit, err := parseExitCode(exitPath, content)
			if err != nil {
				return result, fmt.Errorf("enrich referenced workdir %q: %w", cleanWorkdir, err)
			}
			result.exit = exit
			result.commandFilesAvailable = true
			result.sources = append(result.sources, domain.LogOnlyEvidenceSource{
				Kind:   domain.LogOnlyEvidenceSourceCommand,
				Path:   exitPath,
				Detail: "parsed exit code from .exitcode",
			})
		}

		for _, kind := range []domain.CommandFileKind{domain.CommandFileErr, domain.CommandFileLog} {
			path := filepath.Join(cleanWorkdir, string(kind))
			exists, err := statRegularFile(path)
			if err != nil {
				return result, fmt.Errorf("enrich referenced workdir %q: %w", cleanWorkdir, err)
			}
			if !exists {
				continue
			}

			result.commandFilesAvailable = true
			content, truncated, err := readBoundedText(path)
			if err != nil {
				return result, fmt.Errorf("enrich referenced workdir %q: %w", cleanWorkdir, err)
			}

			detail := fmt.Sprintf("read bounded %s evidence (max %d bytes)", kind, maxBytes)
			if truncated {
				detail = fmt.Sprintf("read bounded %s evidence (max %d bytes; truncated)", kind, maxBytes)
			}
			result.sources = append(result.sources, domain.LogOnlyEvidenceSource{
				Kind:   domain.LogOnlyEvidenceSourceCommand,
				Path:   path,
				Detail: detail,
			})

			if strings.TrimSpace(content) == "" {
				continue
			}
			result.errorBlock = content
			result.errorSummary = SummarizeErrorBlock(content, maxBytesAsInt())
			break
		}

		return result, nil
	}

	enriched := make([]domain.LogOnlyTaskEvidence, len(evidence))
	copy(enriched, evidence)

	byWorkdir := make(map[string]workdirEvidence)
	for index := range enriched {
		if err := ctx.Err(); err != nil {
			return nil, fmt.Errorf("enrich log-only evidence from workdirs: %w", err)
		}

		workdir := strings.TrimSpace(enriched[index].Workdir)
		if workdir == "" {
			continue
		}
		workdir = filepath.Clean(workdir)

		details, ok := byWorkdir[workdir]
		if !ok {
			var err error
			details, err = enrichWorkdir(workdir)
			if err != nil {
				return nil, fmt.Errorf("enrich log-only evidence from workdirs: %w", err)
			}
			byWorkdir[workdir] = details
		}

		changed := false
		if details.exit != nil && enriched[index].Exit == nil {
			exit := *details.exit
			enriched[index].Exit = &exit
			changed = true
		}
		if details.errorBlock != "" {
			if strings.TrimSpace(enriched[index].ErrorBlock) == "" {
				enriched[index].ErrorBlock = details.errorBlock
				changed = true
			}
			if strings.TrimSpace(enriched[index].ErrorSummary) == "" {
				enriched[index].ErrorSummary = details.errorSummary
				changed = true
			}
		}
		if details.commandFilesAvailable && !enriched[index].CommandFilesAvailable {
			enriched[index].CommandFilesAvailable = true
			changed = true
		}
		if len(details.sources) > 0 {
			enriched[index].Sources = append(enriched[index].Sources, details.sources...)
			changed = true
		}
		if changed && enriched[index].Completeness == "" {
			enriched[index].Completeness = domain.LogOnlyEvidencePartial
		}
	}

	return enriched, nil
}

func BuildLogOnlyEvidenceStatus(runDir domain.RunDir, artifacts domain.ArtifactSet, evidence []domain.LogOnlyTaskEvidence) (domain.StatusSummary, error) {
	logOnlyEvidence := make([]domain.LogOnlyTaskEvidence, len(evidence))
	for index := range evidence {
		logOnlyEvidence[index] = evidence[index]
		if evidence[index].Exit != nil {
			exit := *evidence[index].Exit
			logOnlyEvidence[index].Exit = &exit
		}
		if len(evidence[index].Sources) > 0 {
			logOnlyEvidence[index].Sources = append([]domain.LogOnlyEvidenceSource(nil), evidence[index].Sources...)
		}
	}

	selectedLog := selectedLogPath(artifacts)
	searchedPatterns := searchedPatternsText(artifacts)

	isFailureLikeStatus := func(status domain.TaskStatus) bool {
		return status == domain.TaskStatusFailed || status == domain.TaskStatusAborted
	}

	countsByStatus := make(map[domain.TaskStatus]int)
	firstSeenStatuses := make([]domain.TaskStatus, 0)
	failedCount := 0
	logOnlyFailures := make([]domain.LogOnlyFailure, 0)
	for _, item := range logOnlyEvidence {
		if item.ObservedStatus != "" {
			if _, ok := countsByStatus[item.ObservedStatus]; !ok {
				firstSeenStatuses = append(firstSeenStatuses, item.ObservedStatus)
			}
			countsByStatus[item.ObservedStatus]++
		}

		if isFailureLikeStatus(item.ObservedStatus) {
			failedCount++
			logOnlyFailures = append(logOnlyFailures, domain.LogOnlyFailure{
				ID:           item.ID,
				Workdir:      item.Workdir,
				Process:      item.Process,
				Name:         item.Name,
				Exit:         item.Exit,
				ErrorSummary: item.ErrorSummary,
				ErrorBlock:   item.ErrorBlock,
			})
		}
	}

	orderedKnownStatuses := []domain.TaskStatus{
		domain.TaskStatusAborted,
		domain.TaskStatusCached,
		domain.TaskStatusCompleted,
		domain.TaskStatusFailed,
		domain.TaskStatusRunning,
		domain.TaskStatusSubmitted,
		domain.TaskStatusUnknown,
	}
	counts := make([]domain.StatusCount, 0, len(countsByStatus))
	usedStatus := make(map[domain.TaskStatus]bool, len(countsByStatus))
	for _, status := range orderedKnownStatuses {
		if count, ok := countsByStatus[status]; ok {
			counts = append(counts, domain.StatusCount{Status: status, Count: count})
			usedStatus[status] = true
		}
	}
	for _, status := range firstSeenStatuses {
		if usedStatus[status] {
			continue
		}
		counts = append(counts, domain.StatusCount{Status: status, Count: countsByStatus[status]})
		usedStatus[status] = true
	}

	degradedDetail := []string{
		"Selected log: " + selectedLog,
		fmt.Sprintf("Observed log-only evidence rows: %d", len(logOnlyEvidence)),
		"Observed status counts are incomplete because log-only evidence is not a complete task table.",
		"Complete task counts, per-status totals, durations, CPU, and memory data require a Nextflow trace file.",
	}
	if searchedPatterns != "" {
		degradedDetail = append(degradedDetail, "Searched patterns: "+searchedPatterns)
	}

	diagnostics := make([]domain.Diagnostic, 0, len(artifacts.Diagnostics)+3)
	diagnostics = append(diagnostics, artifacts.Diagnostics...)
	diagnostics = append(diagnostics, domain.Diagnostic{
		Severity: domain.DiagnosticWarning,
		Code:     "log_only_degraded",
		Message:  "log-only status is degraded; observed counts are incomplete and complete task/resource/status data is unavailable",
		Detail:   strings.Join(degradedDetail, "\n"),
	})

	if len(logOnlyEvidence) == 0 {
		missingDetail := []string{
			"Selected log: " + selectedLog,
			"No parseable task, lifecycle, failure, or workdir evidence was found in the selected Nextflow log.",
			"The run may have failed before task evidence was emitted, or this log format is unsupported.",
			"complete task/resource/status data is unavailable without a Nextflow trace file.",
		}
		if searchedPatterns != "" {
			missingDetail = append(missingDetail, "Searched patterns: "+searchedPatterns)
		}

		diagnostics = append(diagnostics, domain.Diagnostic{
			Severity: domain.DiagnosticError,
			Code:     "log_only_no_parseable_evidence",
			Message:  "No parseable task evidence found in selected Nextflow log",
			Detail:   strings.Join(missingDetail, "\n"),
		})
	}

	hasTraceRecommendation := false
	for _, diagnostic := range diagnostics {
		if diagnostic.Code == "nextflow_with_trace_recommended" {
			hasTraceRecommendation = true
			break
		}
	}
	if !hasTraceRecommendation {
		diagnostics = append(diagnostics, domain.NextflowTraceRecommendationDiagnostic())
	}

	return domain.StatusSummary{
		RunDir:          runDir,
		Mode:            domain.IndexModeLogOnly,
		Freshness:       domain.IndexFreshnessUnsupported,
		Sources:         artifacts,
		Counts:          counts,
		FailedCount:     failedCount,
		FailedPreview:   []domain.FailedTaskPreview{},
		LogOnlyFailures: logOnlyFailures,
		LogOnlyEvidence: logOnlyEvidence,
		Diagnostics:     diagnostics,
	}, nil
}

func ExtractFailureBlocks(reader io.Reader) ([]FailureBlock, error) {
	text, err := readNormalizedLogText(reader, "extract failure blocks")
	if err != nil {
		return nil, err
	}
	if strings.TrimSpace(text) == "" {
		return []FailureBlock{}, nil
	}

	lines := strings.Split(text, "\n")
	if len(lines) > 0 && lines[len(lines)-1] == "" {
		lines = lines[:len(lines)-1]
	}

	errorExecutingSingle := regexp.MustCompile(`(?i)Error executing process\s*>\s*'([^']+)'`)
	errorExecutingDouble := regexp.MustCompile(`(?i)Error executing process\s*>\s*"([^"]+)"`)
	errorExecutingBacktick := regexp.MustCompile("(?i)Error executing process\\s*>\\s*`([^`]+)`")
	processTerminated := regexp.MustCompile("(?i)Process\\s+[`'\"]([^`'\"]+)[`'\"]\\s+terminated\\s+with\\s+an\\s+error")
	logRecord := regexp.MustCompile(`^(?:[A-Z][a-z]{2}-\d{2}|\d{4}-\d{2}-\d{2})\s+\d{2}:\d{2}:\d{2}(?:\.\d+)?\s+\[[^\]]+\]\s+[A-Z]+(?:\s+\S+)?\s+-\s+`)
	workdirHeader := regexp.MustCompile(`(?i)^\s*(?:work\s*dir|working\s+directory)\s*:\s*(.*)$`)
	exitHeader := regexp.MustCompile(`(?i)^\s*(?:command\s+)?exit\s+(?:status|code)\s*:\s*(.*)$`)
	exitParen := regexp.MustCompile(`(?i)error\s+exit\s+status\s*\(\s*(-?\d+)\s*\)`)

	isFailureStart := func(line string) bool {
		lower := strings.ToLower(line)
		return strings.Contains(lower, "error executing process") || processTerminated.MatchString(line)
	}

	isNextUntimestampedErrorBlock := func(line string) bool {
		trimmed := strings.TrimSpace(line)
		return strings.HasPrefix(trimmed, "ERROR ~") && isFailureStart(line)
	}

	labelFromBlock := func(block string) string {
		for _, pattern := range []*regexp.Regexp{errorExecutingSingle, errorExecutingDouble, errorExecutingBacktick, processTerminated} {
			matches := pattern.FindStringSubmatch(block)
			if len(matches) >= 2 {
				label := strings.TrimSpace(matches[1])
				if label != "" {
					return label
				}
			}
		}
		return ""
	}

	extractWorkdir := func(blockLines []string) string {
		value, ok := firstHeaderValue(blockLines, workdirHeader)
		if !ok {
			return ""
		}
		return cleanSimpleEvidenceValue(value)
	}

	extractExit := func(block string, blockLines []string) *int {
		if matches := exitParen.FindStringSubmatch(block); len(matches) == 2 {
			if value, err := strconv.Atoi(matches[1]); err == nil {
				return &value
			}
		}

		value, ok := firstHeaderValue(blockLines, exitHeader)
		if !ok {
			return nil
		}

		parsed, err := strconv.Atoi(value)
		if err == nil {
			return &parsed
		}
		return nil
	}

	trimOuterBlankLines := func(blockLines []string) []string {
		for len(blockLines) > 0 && strings.TrimSpace(blockLines[0]) == "" {
			blockLines = blockLines[1:]
		}
		for len(blockLines) > 0 && strings.TrimSpace(blockLines[len(blockLines)-1]) == "" {
			blockLines = blockLines[:len(blockLines)-1]
		}
		return blockLines
	}

	blocks := make([]FailureBlock, 0)
	for i := 0; i < len(lines); {
		if !isFailureStart(lines[i]) {
			i++
			continue
		}

		start := i
		end := len(lines)
		for j := i + 1; j < len(lines); j++ {
			if logRecord.MatchString(lines[j]) || isNextUntimestampedErrorBlock(lines[j]) {
				end = j
				break
			}
		}

		blockLines := trimOuterBlankLines(lines[start:end])
		blockText := strings.Join(blockLines, "\n")
		label := labelFromBlock(blockText)
		process, name := splitProcessLabel(label)

		blocks = append(blocks, FailureBlock{
			Process: process,
			Name:    name,
			Workdir: extractWorkdir(blockLines),
			Exit:    extractExit(blockText, blockLines),
			Block:   blockText,
		})

		i = end
	}

	return blocks, nil
}

func NormalizeFailureBlock(runDir domain.RunDir, block FailureBlock) (domain.LogOnlyFailure, error) {
	const errorSummaryMaxBytes = 4096

	blockText := block.Block
	canonicalID := ""
	workdir := ""

	if strings.TrimSpace(block.Workdir) != "" {
		id, err := trace.DeriveCanonicalTaskID(block.Workdir)
		if err != nil {
			return domain.LogOnlyFailure{}, fmt.Errorf("normalize failure block: invalid parsed workdir %q: %w", block.Workdir, err)
		}
		canonicalID = id

		resolvedWorkdir, err := trace.ResolveTaskWorkdir(runDir, canonicalID, block.Workdir)
		if err != nil {
			return domain.LogOnlyFailure{}, fmt.Errorf("normalize failure block: resolve parsed workdir %q: %w", block.Workdir, err)
		}
		workdir = resolvedWorkdir
	} else {
		extractedID, extractedWorkdir, err := ExtractWorkdirEvidence(blockText)
		if err != nil {
			return domain.LogOnlyFailure{}, fmt.Errorf("normalize failure block: extract workdir evidence: %w", err)
		}
		canonicalID = extractedID

		workdirInput := extractedWorkdir
		if strings.TrimSpace(workdirInput) == "" {
			workdirInput = canonicalID
		}
		if canonicalID != "" || strings.TrimSpace(workdirInput) != "" {
			resolvedWorkdir, err := trace.ResolveTaskWorkdir(runDir, canonicalID, workdirInput)
			if err != nil {
				return domain.LogOnlyFailure{}, fmt.Errorf("normalize failure block: resolve extracted workdir %q: %w", workdirInput, err)
			}
			workdir = resolvedWorkdir
		}
	}

	exit := block.Exit
	if exit == nil {
		extractedExit, err := ExtractExitEvidence(blockText)
		if err != nil {
			return domain.LogOnlyFailure{}, fmt.Errorf("normalize failure block: extract exit evidence: %w", err)
		}
		exit = extractedExit
	}

	return domain.LogOnlyFailure{
		ID:           canonicalID,
		Workdir:      workdir,
		Process:      strings.TrimSpace(block.Process),
		Name:         strings.TrimSpace(block.Name),
		Exit:         exit,
		ErrorSummary: SummarizeErrorBlock(blockText, errorSummaryMaxBytes),
		ErrorBlock:   blockText,
	}, nil
}

func ExtractWorkdirEvidence(block string) (id string, workdir string, err error) {
	text := normalizeLogText(block)
	if strings.TrimSpace(text) == "" {
		return "", "", nil
	}

	isHex := func(value string) bool {
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

	cleanValue := func(value string) string {
		value = strings.TrimSpace(value)
		for {
			previous := value
			value = strings.TrimSpace(value)
			value = strings.TrimRight(value, ",;.:")
			value = strings.TrimSpace(value)

			if len(value) >= 2 {
				first := value[0]
				last := value[len(value)-1]
				if (first == '`' && last == '`') ||
					(first == '"' && last == '"') ||
					(first == '\'' && last == '\'') ||
					(first == '[' && last == ']') ||
					(first == '(' && last == ')') ||
					(first == '{' && last == '}') ||
					(first == '<' && last == '>') {
					value = strings.TrimSpace(value[1 : len(value)-1])
				}
			}

			if value == previous {
				return value
			}
		}
	}

	extractWorkdirPath := func(value string) string {
		value = cleanValue(value)
		if value == "" {
			return ""
		}

		normalized := strings.TrimRight(strings.ReplaceAll(value, "\\", "/"), "/")
		if normalized == "" {
			return ""
		}

		parts := strings.Split(normalized, "/")
		for i, part := range parts {
			if !strings.EqualFold(part, "work") || i+2 >= len(parts) {
				continue
			}
			if len(parts[i+1]) == 2 && isHex(parts[i+1]) && isHex(parts[i+2]) {
				return strings.Join(parts[:i+3], "/")
			}
		}

		return ""
	}

	parseEvidence := func(value string, direct bool) (string, string, error) {
		value = cleanValue(value)
		if value == "" {
			if direct {
				return "", "", fmt.Errorf("extract workdir evidence: invalid workdir value: missing value")
			}
			return "", "", nil
		}

		canonicalID, deriveErr := trace.DeriveCanonicalTaskID(value)
		if deriveErr != nil {
			if direct {
				return "", "", fmt.Errorf("extract workdir evidence: invalid workdir value %q: %w", value, deriveErr)
			}
			return "", "", nil
		}

		return canonicalID, extractWorkdirPath(value), nil
	}

	lines := strings.Split(text, "\n")
	workdirHeader := regexp.MustCompile(`(?i)^\s*(?:work\s*dir|working\s+directory)\s*:\s*(.*)$`)
	for i, line := range lines {
		matches := workdirHeader.FindStringSubmatch(line)
		if len(matches) == 0 {
			continue
		}

		value := cleanValue(matches[1])
		if value == "" {
			for j := i + 1; j < len(lines); j++ {
				candidate := cleanValue(lines[j])
				if candidate != "" {
					value = candidate
					break
				}
			}
		}

		return parseEvidence(value, true)
	}

	for _, token := range strings.Fields(text) {
		candidate := cleanValue(token)
		if candidate == "" {
			continue
		}
		if strings.Contains(candidate, "=") {
			candidate = candidate[strings.LastIndex(candidate, "=")+1:]
		}
		if extractWorkdirPath(candidate) == "" {
			continue
		}

		canonicalID, workdirPath, parseErr := parseEvidence(candidate, false)
		if parseErr != nil {
			return "", "", parseErr
		}
		if canonicalID != "" {
			return canonicalID, workdirPath, nil
		}
	}

	return "", "", nil
}

func ExtractExitEvidence(block string) (*int, error) {
	text := normalizeLogText(block)
	if strings.TrimSpace(text) == "" {
		return nil, nil
	}

	lines := strings.Split(text, "\n")
	exitHeader := regexp.MustCompile(`(?i)^\s*(?:command\s+)?exit\s+(?:status|code)\s*:\s*(.*)$`)

	for i, line := range lines {
		matches := exitHeader.FindStringSubmatch(line)
		if len(matches) == 0 {
			continue
		}

		value := cleanSimpleEvidenceValue(matches[1])
		if value == "" {
			for j := i + 1; j < len(lines); j++ {
				candidate := cleanSimpleEvidenceValue(lines[j])
				if candidate != "" {
					value = candidate
					break
				}
			}
		}
		if value == "" {
			return nil, fmt.Errorf("extract exit evidence: invalid exit value: missing value after %q", strings.TrimSpace(line))
		}

		parsed, err := strconv.Atoi(value)
		if err != nil {
			return nil, fmt.Errorf("extract exit evidence: invalid exit value %q: %w", value, err)
		}
		return &parsed, nil
	}

	exitParen := regexp.MustCompile(`(?i)\berror\s+exit\s+status\s*\(\s*(-?\d+)\s*\)`)
	if matches := exitParen.FindStringSubmatch(text); len(matches) == 2 {
		parsed, err := strconv.Atoi(matches[1])
		if err != nil {
			return nil, fmt.Errorf("extract exit evidence: invalid exit value %q: %w", matches[1], err)
		}
		return &parsed, nil
	}

	return nil, nil
}

func SummarizeErrorBlock(block string, maxBytes int) string {
	if maxBytes <= 0 {
		return ""
	}

	text := normalizeLogText(block)
	if strings.TrimSpace(text) == "" {
		return ""
	}

	lines := strings.Split(text, "\n")
	commandErrorHeader := regexp.MustCompile(`(?i)^\s*command\s+error\s*:\s*(.*)$`)
	isMetadataHeader := func(line string) bool {
		lower := strings.ToLower(strings.TrimSpace(line))
		return strings.HasPrefix(lower, "work dir:") ||
			strings.HasPrefix(lower, "working directory:") ||
			strings.HasPrefix(lower, "tip:") ||
			strings.HasPrefix(lower, "command executed:") ||
			strings.HasPrefix(lower, "command output:") ||
			strings.HasPrefix(lower, "command exit status:") ||
			strings.HasPrefix(lower, "command exit code:") ||
			strings.HasPrefix(lower, "exit status:") ||
			strings.HasPrefix(lower, "exit code:")
	}

	collectCommandError := func() []string {
		for i, line := range lines {
			matches := commandErrorHeader.FindStringSubmatch(line)
			if len(matches) == 0 {
				continue
			}

			summaryLines := make([]string, 0)
			if value := strings.TrimSpace(matches[1]); value != "" {
				summaryLines = append(summaryLines, value)
			}
			for j := i + 1; j < len(lines); j++ {
				trimmed := strings.TrimSpace(lines[j])
				if trimmed == "" {
					continue
				}
				if isMetadataHeader(trimmed) {
					break
				}
				summaryLines = append(summaryLines, trimmed)
			}
			return summaryLines
		}
		return nil
	}

	summaryLines := collectCommandError()
	if len(summaryLines) == 0 {
		skipNextWorkdirValue := false
		for _, line := range lines {
			trimmed := strings.TrimSpace(line)
			if trimmed == "" {
				skipNextWorkdirValue = false
				continue
			}

			lower := strings.ToLower(trimmed)
			if skipNextWorkdirValue {
				skipNextWorkdirValue = false
				continue
			}
			if strings.HasPrefix(lower, "tip:") {
				continue
			}
			if strings.HasPrefix(lower, "work dir:") || strings.HasPrefix(lower, "working directory:") {
				parts := strings.SplitN(trimmed, ":", 2)
				if len(parts) == 2 && strings.TrimSpace(parts[1]) == "" {
					skipNextWorkdirValue = true
				}
				continue
			}

			summaryLines = append(summaryLines, trimmed)
		}
	}

	summary := strings.Join(summaryLines, "\n")
	if len(summary) <= maxBytes {
		return summary
	}

	cut := 0
	for i := range summary {
		if i > maxBytes {
			break
		}
		cut = i
	}
	return summary[:cut]
}

func BuildLogOnlyStatus(runDir domain.RunDir, artifacts domain.ArtifactSet, failures []domain.LogOnlyFailure) (domain.StatusSummary, error) {
	logOnlyFailures := make([]domain.LogOnlyFailure, len(failures))
	copy(logOnlyFailures, failures)

	selectedLog := selectedLogPath(artifacts)
	searchedPatterns := searchedPatternsText(artifacts)

	degradedDetail := []string{
		"Selected log: " + selectedLog,
		"Only deterministic failure evidence parsed from the Nextflow log is available.",
		"Complete task counts, per-status totals, durations, CPU, and memory data require a trace file.",
	}
	if searchedPatterns != "" {
		degradedDetail = append(degradedDetail, "Searched patterns: "+searchedPatterns)
	}

	diagnostics := make([]domain.Diagnostic, 0, len(artifacts.Diagnostics)+3)
	diagnostics = append(diagnostics, artifacts.Diagnostics...)
	diagnostics = append(diagnostics, domain.Diagnostic{
		Severity: domain.DiagnosticWarning,
		Code:     "log_only_degraded",
		Message:  "log-only status is degraded; complete task/resource/status data is unavailable",
		Detail:   strings.Join(degradedDetail, "\n"),
	})

	if len(logOnlyFailures) == 0 {
		missingDetail := []string{
			"Selected log: " + selectedLog,
			"The run may have failed before a task failure block was emitted, or this log format is unsupported.",
			"complete task/resource/status data is unavailable without a Nextflow trace file.",
		}
		if searchedPatterns != "" {
			missingDetail = append(missingDetail, "Searched patterns: "+searchedPatterns)
		}

		diagnostics = append(diagnostics,
			domain.Diagnostic{
				Severity: domain.DiagnosticError,
				Code:     "log_only_no_parseable_failures",
				Message:  "No parseable task failure evidence found in selected Nextflow log",
				Detail:   strings.Join(missingDetail, "\n"),
			},
			domain.NextflowTraceRecommendationDiagnostic(),
		)
	}

	return domain.StatusSummary{
		RunDir:          runDir,
		Mode:            domain.IndexModeLogOnly,
		Freshness:       domain.IndexFreshnessUnsupported,
		Sources:         artifacts,
		Counts:          []domain.StatusCount{},
		FailedCount:     len(logOnlyFailures),
		FailedPreview:   []domain.FailedTaskPreview{},
		LogOnlyFailures: logOnlyFailures,
		Diagnostics:     diagnostics,
	}, nil
}
