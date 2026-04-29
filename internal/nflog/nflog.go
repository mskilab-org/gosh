package nflog

import (
	"context"
	"fmt"
	"io"
	"os"
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

func ExtractFailureBlocks(reader io.Reader) ([]FailureBlock, error) {
	if reader == nil {
		return nil, fmt.Errorf("extract failure blocks: nil reader")
	}

	data, err := io.ReadAll(reader)
	if err != nil {
		return nil, fmt.Errorf("extract failure blocks: read log: %w", err)
	}

	text := normalizeLogText(string(data))
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

	splitProcessLabel := func(label string) (string, string) {
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

	selectedLog := "selected Nextflow log"
	if artifacts.Log != nil && artifacts.Log.Path != "" {
		selectedLog = artifacts.Log.Path
	}

	searchedPatterns := ""
	if len(artifacts.SearchedPatterns) > 0 {
		searchedPatterns = strings.Join(artifacts.SearchedPatterns, ", ")
	}

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
			domain.Diagnostic{
				Severity: domain.DiagnosticInfo,
				Code:     "nextflow_with_trace_recommended",
				Message:  "Run future Nextflow workflows with -with-trace",
				Detail:   "Use `nextflow run ... -with-trace` for future runs so gosh can build a complete trace-backed task index.",
			},
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
