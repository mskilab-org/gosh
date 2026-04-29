package render

import (
	"encoding/json"
	"fmt"
	"io"
	"sort"
	"strings"
	"time"

	"github.com/mskilab-org/gosh/internal/domain"
)

type diagnosticJSON struct {
	Severity string `json:"severity"`
	Code     string `json:"code"`
	Message  string `json:"message"`
	Detail   string `json:"detail"`
}

type sourceFingerprintJSON struct {
	Kind    string  `json:"kind"`
	Path    string  `json:"path"`
	ModTime *string `json:"mod_time"`
	Size    int64   `json:"size"`
}

type artifactSetJSON struct {
	RunDir           string                 `json:"run_dir"`
	Mode             string                 `json:"mode"`
	Trace            *sourceFingerprintJSON `json:"trace"`
	Log              *sourceFingerprintJSON `json:"log"`
	SelectedAt       *string                `json:"selected_at"`
	SearchedPatterns []string               `json:"searched_patterns"`
	Diagnostics      []diagnosticJSON       `json:"diagnostics"`
}

type indexMetadataJSON struct {
	SchemaVersion int                    `json:"schema_version"`
	RunDir        string                 `json:"run_dir"`
	IndexPath     string                 `json:"index_path"`
	Mode          string                 `json:"mode"`
	Trace         *sourceFingerprintJSON `json:"trace"`
	Log           *sourceFingerprintJSON `json:"log"`
	BuiltAt       *string                `json:"built_at"`
	Freshness     string                 `json:"freshness"`
	StaleReason   string                 `json:"stale_reason"`
	TaskCount     int                    `json:"task_count"`
}

type taskJSON struct {
	ID           string `json:"id"`
	RowOrder     int64  `json:"row_order"`
	Status       string `json:"status"`
	Process      string `json:"process"`
	Name         string `json:"name"`
	Tag          string `json:"tag"`
	Workdir      string `json:"workdir"`
	Exit         *int   `json:"exit"`
	Duration     string `json:"duration"`
	Realtime     string `json:"realtime"`
	CPUs         string `json:"cpus"`
	Memory       string `json:"memory"`
	ErrorSummary string `json:"error_summary"`
}

func jsonTime(value time.Time) *string {
	if value.IsZero() {
		return nil
	}
	formatted := value.UTC().Format(time.RFC3339Nano)
	return &formatted
}

func jsonTimePtr(value *time.Time) *string {
	if value == nil {
		return nil
	}
	return jsonTime(*value)
}

func diagnosticsJSON(diagnostics []domain.Diagnostic) []diagnosticJSON {
	out := make([]diagnosticJSON, len(diagnostics))
	for i, diagnostic := range diagnostics {
		out[i] = diagnosticJSON{
			Severity: string(diagnostic.Severity),
			Code:     diagnostic.Code,
			Message:  diagnostic.Message,
			Detail:   diagnostic.Detail,
		}
	}
	return out
}

func sourceFingerprintJSONFor(source *domain.SourceFingerprint) *sourceFingerprintJSON {
	if source == nil {
		return nil
	}
	return &sourceFingerprintJSON{
		Kind:    string(source.Kind),
		Path:    source.Path,
		ModTime: jsonTime(source.ModTime),
		Size:    source.Size,
	}
}

func stringsJSON(values []string) []string {
	out := make([]string, len(values))
	copy(out, values)
	return out
}

func artifactSetJSONFor(artifacts domain.ArtifactSet) artifactSetJSON {
	return artifactSetJSON{
		RunDir:           artifacts.RunDir.Path,
		Mode:             string(artifacts.Mode),
		Trace:            sourceFingerprintJSONFor(artifacts.Trace),
		Log:              sourceFingerprintJSONFor(artifacts.Log),
		SelectedAt:       jsonTime(artifacts.SelectedAt),
		SearchedPatterns: stringsJSON(artifacts.SearchedPatterns),
		Diagnostics:      diagnosticsJSON(artifacts.Diagnostics),
	}
}

func indexMetadataJSONFor(metadata *domain.IndexMetadata) *indexMetadataJSON {
	if metadata == nil {
		return nil
	}
	return &indexMetadataJSON{
		SchemaVersion: metadata.SchemaVersion,
		RunDir:        metadata.RunDir,
		IndexPath:     metadata.IndexPath,
		Mode:          string(metadata.Mode),
		Trace:         sourceFingerprintJSONFor(metadata.Trace),
		Log:           sourceFingerprintJSONFor(metadata.Log),
		BuiltAt:       jsonTime(metadata.BuiltAt),
		Freshness:     string(metadata.Freshness),
		StaleReason:   metadata.StaleReason,
		TaskCount:     metadata.TaskCount,
	}
}

func taskValueJSON(task domain.Task) taskJSON {
	return taskJSON{
		ID:           task.ID,
		RowOrder:     task.RowOrder,
		Status:       string(task.Status),
		Process:      task.Process,
		Name:         task.Name,
		Tag:          task.Tag,
		Workdir:      task.Workdir,
		Exit:         task.Exit,
		Duration:     task.Duration,
		Realtime:     task.Realtime,
		CPUs:         task.CPUs,
		Memory:       task.Memory,
		ErrorSummary: task.ErrorSummary,
	}
}

func taskPointerJSON(task *domain.Task) *taskJSON {
	if task == nil {
		return nil
	}
	out := taskValueJSON(*task)
	return &out
}

func tasksJSON(tasks []domain.Task) []taskJSON {
	out := make([]taskJSON, len(tasks))
	for i, task := range tasks {
		out[i] = taskValueJSON(task)
	}
	return out
}

func normalizeLineBreaks(value string) string {
	value = strings.ReplaceAll(value, "\r\n", "\n")
	return strings.ReplaceAll(value, "\r", "\n")
}

func inlineText(value string) string {
	fields := strings.Fields(normalizeLineBreaks(value))
	if len(fields) == 0 {
		return ""
	}
	return strings.Join(fields, " ")
}

func humanField(value string) string {
	text := inlineText(value)
	if text == "" {
		return "-"
	}
	return text
}

func formatOptionalInt(value *int) string {
	if value == nil {
		return "-"
	}
	return fmt.Sprintf("%d", *value)
}

func snippetLineRangeText(snippet *domain.Snippet) string {
	if snippet == nil || snippet.StartLine <= 0 || snippet.EndLine <= 0 {
		return "-"
	}
	return fmt.Sprintf("%d-%d", snippet.StartLine, snippet.EndLine)
}

func commandFileKindRank(kind domain.CommandFileKind) int {
	switch kind {
	case domain.CommandFileShell:
		return 0
	case domain.CommandFileLog:
		return 1
	case domain.CommandFileErr:
		return 2
	case domain.CommandFileOut:
		return 3
	case domain.CommandFileRun:
		return 4
	default:
		return 5
	}
}

func orderedCommandFiles(files []domain.CommandFile) []domain.CommandFile {
	ordered := append([]domain.CommandFile(nil), files...)
	sort.SliceStable(ordered, func(i, j int) bool {
		leftRank := commandFileKindRank(ordered[i].Kind)
		rightRank := commandFileKindRank(ordered[j].Kind)
		if leftRank != rightRank {
			return leftRank < rightRank
		}
		leftKind := string(ordered[i].Kind)
		rightKind := string(ordered[j].Kind)
		if leftKind != rightKind {
			return leftKind < rightKind
		}
		return ordered[i].Path < ordered[j].Path
	})
	return ordered
}

func collectDiagnostics(parts ...[]domain.Diagnostic) []domain.Diagnostic {
	diagnostics := make([]domain.Diagnostic, 0)
	seen := make(map[string]bool)
	for _, part := range parts {
		for _, diagnostic := range part {
			key := string(diagnostic.Severity) + "\x00" + diagnostic.Code + "\x00" + diagnostic.Message + "\x00" + diagnostic.Detail
			if seen[key] {
				continue
			}
			seen[key] = true
			diagnostics = append(diagnostics, diagnostic)
		}
	}
	return diagnostics
}

func appendDiagnosticBlocksHuman(builder *strings.Builder, diagnostics []domain.Diagnostic, context string) error {
	blocks := BuildDiagnosticBlocks(diagnostics)
	if len(blocks) == 0 {
		return nil
	}

	var diagnosticsBuilder strings.Builder
	if err := RenderDiagnosticBlocksHuman(&diagnosticsBuilder, blocks); err != nil {
		if context == "" {
			return err
		}
		return fmt.Errorf("%s: %w", context, err)
	}
	text := diagnosticsBuilder.String()
	if strings.TrimSpace(text) == "" {
		return nil
	}

	if builder.Len() > 0 {
		current := builder.String()
		if !strings.HasSuffix(current, "\n") {
			builder.WriteString("\n")
		}
		if !strings.HasSuffix(current, "\n\n") {
			builder.WriteString("\n")
		}
	}
	builder.WriteString(text)
	return nil
}

func writeIndentedHumanContent(builder *strings.Builder, indent string, content string) {
	content = normalizeLineBreaks(content)
	content = strings.TrimSuffix(content, "\n")
	if content == "" {
		fmt.Fprintf(builder, "%s-\n", indent)
		return
	}
	for _, line := range strings.Split(content, "\n") {
		fmt.Fprintf(builder, "%s%s\n", indent, line)
	}
}

func writeCommandFilesHuman(builder *strings.Builder, files []domain.CommandFile, snippetIndent string) {
	if len(files) == 0 {
		builder.WriteString("command_files: none\n")
		return
	}

	ordered := orderedCommandFiles(files)

	builder.WriteString("command_files:\n")
	for _, file := range ordered {
		path := file.Path
		if strings.TrimSpace(path) == "" && file.Snippet != nil {
			path = file.Snippet.Path
		}
		fmt.Fprintf(builder, "  - kind=%s path=%s exists=%t size=%d\n",
			humanField(string(file.Kind)),
			humanField(path),
			file.Exists,
			file.Size,
		)
		if file.Snippet == nil {
			continue
		}
		fmt.Fprintf(builder, "    snippet: strategy=%s lines=%s truncated=%t max_bytes=%d\n",
			humanField(string(file.Snippet.Strategy)),
			snippetLineRangeText(file.Snippet),
			file.Snippet.Truncated,
			file.Snippet.MaxBytes,
		)
		builder.WriteString("    content:\n")
		writeIndentedHumanContent(builder, snippetIndent, file.Snippet.Content)
	}
}

func writeIndentedJSON(writer io.Writer, payload any, context string) error {
	encoder := json.NewEncoder(writer)
	encoder.SetEscapeHTML(false)
	encoder.SetIndent("", "  ")
	if err := encoder.Encode(payload); err != nil {
		return fmt.Errorf("%s: write: %w", context, err)
	}
	return nil
}

func BuildDiagnosticBlocks(diagnostics []domain.Diagnostic) []domain.DiagnosticBlock {
	if len(diagnostics) == 0 {
		return []domain.DiagnosticBlock{}
	}

	inline := inlineText
	severityRank := func(severity domain.DiagnosticSeverity) int {
		switch severity {
		case domain.DiagnosticError:
			return 0
		case domain.DiagnosticWarning:
			return 1
		case domain.DiagnosticInfo:
			return 2
		default:
			return 3
		}
	}
	appendUnique := func(values []string, additions ...string) []string {
		for _, addition := range additions {
			addition = inline(addition)
			if addition == "" {
				continue
			}
			seen := false
			for _, value := range values {
				if value == addition {
					seen = true
					break
				}
			}
			if !seen {
				values = append(values, addition)
			}
		}
		return values
	}
	isHintText := func(value string) bool {
		lower := strings.ToLower(strings.TrimSpace(value))
		return strings.HasPrefix(lower, "use ") ||
			strings.HasPrefix(lower, "try ") ||
			strings.HasPrefix(lower, "run ") ||
			strings.HasPrefix(lower, "refresh ") ||
			strings.HasPrefix(lower, "re-run ") ||
			strings.HasPrefix(lower, "rerun ")
	}
	splitContextLine := func(line string) (domain.DiagnosticContextLine, bool) {
		index := strings.Index(line, ":")
		if index <= 0 || index >= len(line)-1 {
			return domain.DiagnosticContextLine{}, false
		}
		label := inline(line[:index])
		value := inline(line[index+1:])
		if label == "" || value == "" {
			return domain.DiagnosticContextLine{}, false
		}
		return domain.DiagnosticContextLine{Label: label, Value: value}, true
	}
	buildBlock := func(diagnostic domain.Diagnostic) domain.DiagnosticBlock {
		title := inline(diagnostic.Message)
		if title == "" {
			title = inline(diagnostic.Code)
		}
		if title == "" {
			title = inline(string(diagnostic.Severity))
		}

		block := domain.DiagnosticBlock{
			Severity: diagnostic.Severity,
			Code:     diagnostic.Code,
			Title:    title,
		}

		for _, rawLine := range strings.Split(normalizeLineBreaks(diagnostic.Detail), "\n") {
			line := inline(rawLine)
			if line == "" {
				continue
			}
			if labeled, ok := splitContextLine(line); ok {
				switch strings.ToLower(labeled.Label) {
				case "hint", "try":
					block.Hints = appendUnique(block.Hints, labeled.Value)
				default:
					block.Context = append(block.Context, labeled)
				}
				continue
			}
			if isHintText(line) {
				block.Hints = appendUnique(block.Hints, line)
				continue
			}
			block.Details = append(block.Details, line)
		}

		return block
	}
	isStandaloneHintBlock := func(block domain.DiagnosticBlock) bool {
		if block.Severity != domain.DiagnosticInfo {
			return false
		}
		if strings.EqualFold(strings.TrimSpace(block.Code), "nextflow_with_trace_recommended") {
			return true
		}
		return len(block.Context) == 0 && len(block.Details) == 0 && (len(block.Hints) > 0 || isHintText(block.Title))
	}
	hintsFromBlock := func(block domain.DiagnosticBlock) []string {
		hints := appendUnique(nil, block.Hints...)
		if len(hints) == 0 {
			hints = appendUnique(hints, block.Title)
		}
		return hints
	}

	blocks := make([]domain.DiagnosticBlock, 0, len(diagnostics))
	hintBlocks := make([]domain.DiagnosticBlock, 0)
	for _, diagnostic := range diagnostics {
		block := buildBlock(diagnostic)
		if isStandaloneHintBlock(block) {
			hintBlocks = append(hintBlocks, block)
			continue
		}
		blocks = append(blocks, block)
	}

	sort.SliceStable(blocks, func(i, j int) bool {
		return severityRank(blocks[i].Severity) < severityRank(blocks[j].Severity)
	})
	if len(blocks) == 0 {
		return hintBlocks
	}
	for _, hintBlock := range hintBlocks {
		blocks[0].Hints = appendUnique(blocks[0].Hints, hintsFromBlock(hintBlock)...)
	}
	return blocks
}

func RenderDiagnosticBlocksHuman(writer io.Writer, blocks []domain.DiagnosticBlock) error {
	if writer == nil {
		return fmt.Errorf("render diagnostic blocks human: nil writer")
	}

	inline := inlineText
	severityText := func(severity domain.DiagnosticSeverity) string {
		switch severity {
		case domain.DiagnosticError:
			return "error"
		case domain.DiagnosticWarning:
			return "warning"
		case domain.DiagnosticInfo:
			return "info"
		default:
			text := inline(string(severity))
			if text == "" {
				return "info"
			}
			return text
		}
	}
	primaryTitle := func(block domain.DiagnosticBlock) string {
		if title := inline(block.Title); title != "" {
			return title
		}
		if code := inline(block.Code); code != "" {
			return code
		}
		return "-"
	}

	var builder strings.Builder
	for index, block := range blocks {
		if index > 0 {
			builder.WriteString("\n")
		}

		title := primaryTitle(block)
		fmt.Fprintf(&builder, "%s: %s\n", severityText(block.Severity), title)

		code := inline(block.Code)
		if code != "" && code != title {
			fmt.Fprintf(&builder, "  code: %s\n", code)
		}

		for _, context := range block.Context {
			label := inline(context.Label)
			value := inline(context.Value)
			if label == "" || value == "" {
				continue
			}
			fmt.Fprintf(&builder, "  %s: %s\n", label, value)
		}

		for _, detail := range block.Details {
			text := inline(detail)
			if text == "" {
				continue
			}
			fmt.Fprintf(&builder, "  %s\n", text)
		}

		for _, hint := range block.Hints {
			text := inline(hint)
			if text == "" {
				continue
			}
			fmt.Fprintf(&builder, "hint: %s\n", text)
		}
	}

	if builder.Len() == 0 {
		return nil
	}
	if _, err := io.WriteString(writer, builder.String()); err != nil {
		return fmt.Errorf("render diagnostic blocks human: write: %w", err)
	}
	return nil
}

func RenderLogOnlyInspectHuman(writer io.Writer, view domain.InspectView) error {
	if writer == nil {
		return fmt.Errorf("render log-only inspect human: nil writer")
	}

	inline := humanField
	exitText := formatOptionalInt
	writeSources := func(builder *strings.Builder, sources []domain.LogOnlyEvidenceSource) {
		if len(sources) == 0 {
			builder.WriteString("sources: none\n")
			return
		}

		builder.WriteString("sources:\n")
		for _, source := range sources {
			fmt.Fprintf(builder, "  - kind=%s path=%s detail=%s\n",
				inline(string(source.Kind)),
				inline(source.Path),
				inline(source.Detail),
			)
		}
	}
	commandFilesAvailable := func(evidence domain.LogOnlyTaskEvidence, inventory domain.CommandFileInventory) bool {
		if evidence.CommandFilesAvailable {
			return true
		}
		for _, file := range inventory.Files {
			if file.Exists {
				return true
			}
		}
		return false
	}
	workdirPath := func(evidence domain.LogOnlyTaskEvidence, inventory domain.CommandFileInventory) string {
		if strings.TrimSpace(evidence.Workdir) != "" {
			return evidence.Workdir
		}
		return inventory.Workdir
	}
	appendDiagnostics := func(builder *strings.Builder, diagnostics []domain.Diagnostic) error {
		return appendDiagnosticBlocksHuman(builder, diagnostics, "render log-only inspect human diagnostics")
	}

	selector := view.Resolution.Selector
	kind := view.Resolution.Kind
	var evidence *domain.LogOnlyTaskEvidence
	matches := []domain.LogOnlyTaskEvidence(nil)
	resolutionDiagnostics := view.Resolution.Diagnostics
	if view.LogOnlyResolution != nil {
		selector = view.LogOnlyResolution.Selector
		kind = view.LogOnlyResolution.Kind
		evidence = view.LogOnlyResolution.Evidence
		matches = view.LogOnlyResolution.Matches
		resolutionDiagnostics = view.LogOnlyResolution.Diagnostics
	}
	if kind == "" && view.LogOnlyDossier != nil {
		kind = domain.SelectorResolutionExact
	}

	evidenceKind := view.EvidenceKind
	if evidenceKind == "" {
		evidenceKind = domain.InspectEvidenceLogOnly
	}

	var builder strings.Builder
	fmt.Fprintf(&builder, "selector: %s\n", inline(selector))
	fmt.Fprintf(&builder, "resolution: %s\n", inline(string(kind)))
	fmt.Fprintf(&builder, "evidence_kind: %s\n", inline(string(evidenceKind)))

	switch kind {
	case domain.SelectorResolutionExact:
		inventory := domain.CommandFileInventory{}
		dossierDiagnostics := []domain.Diagnostic(nil)
		hasEvidence := false
		var selectedEvidence domain.LogOnlyTaskEvidence
		if view.LogOnlyDossier != nil {
			selectedEvidence = view.LogOnlyDossier.Evidence
			inventory = view.LogOnlyDossier.Inventory
			dossierDiagnostics = view.LogOnlyDossier.Diagnostics
			hasEvidence = true
		} else if evidence != nil {
			selectedEvidence = *evidence
			hasEvidence = true
		}

		if !hasEvidence {
			builder.WriteString("dossier: unavailable\n")
			if err := appendDiagnostics(&builder, collectDiagnostics(resolutionDiagnostics, view.Diagnostics)); err != nil {
				return err
			}
			break
		}

		path := workdirPath(selectedEvidence, inventory)
		available := strings.TrimSpace(path) != ""
		builder.WriteString("evidence:\n")
		fmt.Fprintf(&builder, "  id: %s\n", inline(selectedEvidence.ID))
		fmt.Fprintf(&builder, "  observed_status: %s\n", inline(string(selectedEvidence.ObservedStatus)))
		fmt.Fprintf(&builder, "  process: %s\n", inline(selectedEvidence.Process))
		fmt.Fprintf(&builder, "  name: %s\n", inline(selectedEvidence.Name))
		fmt.Fprintf(&builder, "  workdir: %s\n", inline(selectedEvidence.Workdir))
		fmt.Fprintf(&builder, "  exit: %s\n", exitText(selectedEvidence.Exit))
		fmt.Fprintf(&builder, "  completeness: %s\n", inline(string(selectedEvidence.Completeness)))
		fmt.Fprintf(&builder, "  command_files_available: %t\n", commandFilesAvailable(selectedEvidence, inventory))
		fmt.Fprintf(&builder, "  error_summary: %s\n", inline(selectedEvidence.ErrorSummary))
		builder.WriteString("  error_block:\n")
		writeIndentedHumanContent(&builder, "    ", selectedEvidence.ErrorBlock)
		writeSources(&builder, selectedEvidence.Sources)
		builder.WriteString("workdir:\n")
		fmt.Fprintf(&builder, "  path: %s\n", inline(path))
		fmt.Fprintf(&builder, "  available: %t\n", available)
		writeCommandFilesHuman(&builder, inventory.Files, "      ")

		if err := appendDiagnostics(&builder, collectDiagnostics(dossierDiagnostics, resolutionDiagnostics, view.Diagnostics)); err != nil {
			return err
		}
	case domain.SelectorResolutionAmbiguous:
		if len(matches) == 0 {
			builder.WriteString("matches: none\n")
		} else {
			builder.WriteString("matches:\n")
			builder.WriteString("id\tstatus\tprocess\tname\tworkdir\tcommand_files_available\n")
			for _, match := range matches {
				fmt.Fprintf(&builder, "%s\t%s\t%s\t%s\t%s\t%t\n",
					inline(match.ID),
					inline(string(match.ObservedStatus)),
					inline(match.Process),
					inline(match.Name),
					inline(match.Workdir),
					match.CommandFilesAvailable,
				)
			}
		}
		if err := appendDiagnostics(&builder, collectDiagnostics(resolutionDiagnostics, view.Diagnostics)); err != nil {
			return err
		}
	case domain.SelectorResolutionNotFound:
		builder.WriteString("matches: none\n")
		if err := appendDiagnostics(&builder, collectDiagnostics(resolutionDiagnostics, view.Diagnostics)); err != nil {
			return err
		}
	default:
		if err := appendDiagnostics(&builder, collectDiagnostics(resolutionDiagnostics, view.Diagnostics)); err != nil {
			return err
		}
	}

	if _, err := io.WriteString(writer, builder.String()); err != nil {
		return fmt.Errorf("render log-only inspect human: write: %w", err)
	}
	return nil
}

func RenderStatusHuman(writer io.Writer, view domain.StatusView) error {
	if writer == nil {
		return fmt.Errorf("render status human: nil writer")
	}

	summary := view.Summary
	inline := humanField
	timestamp := func(value time.Time) string {
		if value.IsZero() {
			return "-"
		}
		return value.UTC().Format(time.RFC3339)
	}
	exitText := formatOptionalInt
	orderedCounts := func(counts []domain.StatusCount) []domain.StatusCount {
		ordered := append([]domain.StatusCount(nil), counts...)
		sort.SliceStable(ordered, func(i, j int) bool {
			left := string(ordered[i].Status)
			right := string(ordered[j].Status)
			if left == right {
				return ordered[i].Count < ordered[j].Count
			}
			return left < right
		})
		return ordered
	}
	appendDiagnostics := func(builder *strings.Builder, diagnostics []domain.Diagnostic) error {
		return appendDiagnosticBlocksHuman(builder, diagnostics, "render status human diagnostics")
	}
	writeObservedCounts := func(builder *strings.Builder, counts []domain.StatusCount) {
		for _, count := range counts {
			fmt.Fprintf(builder, "  %s: %d\n", inline(string(count.Status)), count.Count)
		}
	}
	writeLogOnlyEvidence := func(builder *strings.Builder) {
		if len(summary.LogOnlyEvidence) > 0 {
			builder.WriteString("log_only_evidence:\n")
			for _, evidence := range summary.LogOnlyEvidence {
				fmt.Fprintf(builder, "  - id=%s status=%s process=%s name=%s workdir=%s exit=%s completeness=%s command_files_available=%t error=%s\n",
					inline(evidence.ID),
					inline(string(evidence.ObservedStatus)),
					inline(evidence.Process),
					inline(evidence.Name),
					inline(evidence.Workdir),
					exitText(evidence.Exit),
					inline(string(evidence.Completeness)),
					evidence.CommandFilesAvailable,
					inline(evidence.ErrorSummary),
				)
			}
			return
		}

		if len(summary.LogOnlyFailures) == 0 {
			builder.WriteString("log_only_evidence: none\n")
			return
		}

		builder.WriteString("log_only_evidence:\n")
		for _, failure := range summary.LogOnlyFailures {
			fmt.Fprintf(builder, "  - id=%s status=FAILED process=%s name=%s workdir=%s exit=%s completeness=partial command_files_available=false error=%s\n",
				inline(failure.ID),
				inline(failure.Process),
				inline(failure.Name),
				inline(failure.Workdir),
				exitText(failure.Exit),
				inline(failure.ErrorSummary),
			)
		}
	}

	var builder strings.Builder
	fmt.Fprintf(&builder, "run_dir: %s\n", inline(summary.RunDir.Path))
	fmt.Fprintf(&builder, "mode: %s\n", inline(string(summary.Mode)))
	if strings.TrimSpace(summary.IndexPath) != "" {
		fmt.Fprintf(&builder, "index: %s\n", strings.TrimSpace(summary.IndexPath))
	}
	if summary.Freshness != "" {
		fmt.Fprintf(&builder, "freshness: %s\n", summary.Freshness)
	}
	if summary.BuiltAt != nil {
		fmt.Fprintf(&builder, "built_at: %s\n", timestamp(*summary.BuiltAt))
	}

	writeSource := func(label string, source *domain.SourceFingerprint) {
		if source == nil || strings.TrimSpace(source.Path) == "" {
			fmt.Fprintf(&builder, "  %s: none\n", label)
			return
		}
		fmt.Fprintf(&builder, "  %s: %s (mtime=%s size=%d)\n", label, strings.TrimSpace(source.Path), timestamp(source.ModTime), source.Size)
	}

	builder.WriteString("sources:\n")
	writeSource("trace", summary.Sources.Trace)
	writeSource("log", summary.Sources.Log)
	if len(summary.Sources.SearchedPatterns) > 0 {
		patterns := make([]string, 0, len(summary.Sources.SearchedPatterns))
		for _, pattern := range summary.Sources.SearchedPatterns {
			if strings.TrimSpace(pattern) != "" {
				patterns = append(patterns, strings.TrimSpace(pattern))
			}
		}
		if len(patterns) > 0 {
			fmt.Fprintf(&builder, "  searched_patterns: %s\n", strings.Join(patterns, ", "))
		}
	}

	counts := orderedCounts(summary.Counts)
	if summary.Mode == domain.IndexModeLogOnly {
		if len(counts) == 0 {
			builder.WriteString("counts: unavailable (log-only mode; trace file required)\n")
		} else {
			builder.WriteString("counts: incomplete (log-only evidence; trace file required)\n")
			builder.WriteString("observed_counts:\n")
			writeObservedCounts(&builder, counts)
		}
	} else if len(counts) == 0 {
		builder.WriteString("counts: none\n")
	} else {
		builder.WriteString("counts:\n")
		writeObservedCounts(&builder, counts)
	}

	fmt.Fprintf(&builder, "failed_count: %d\n", summary.FailedCount)
	if summary.Mode == domain.IndexModeLogOnly {
		writeLogOnlyEvidence(&builder)
	} else if len(summary.FailedPreview) == 0 {
		builder.WriteString("failed_preview: none\n")
	} else {
		builder.WriteString("failed_preview:\n")
		for _, failure := range summary.FailedPreview {
			fmt.Fprintf(&builder, "  - id=%s status=%s process=%s name=%s tag=%s workdir=%s exit=%s error=%s\n",
				inline(failure.ID),
				inline(string(failure.Status)),
				inline(failure.Process),
				inline(failure.Name),
				inline(failure.Tag),
				inline(failure.Workdir),
				exitText(failure.Exit),
				inline(failure.ErrorSummary),
			)
		}
	}

	if err := appendDiagnostics(&builder, summary.Diagnostics); err != nil {
		return err
	}

	if _, err := io.WriteString(writer, builder.String()); err != nil {
		return fmt.Errorf("render status human: write: %w", err)
	}
	return nil
}

func RenderStatusJSON(writer io.Writer, view domain.StatusView) error {
	if writer == nil {
		return fmt.Errorf("render status json: nil writer")
	}

	type statusCountJSON struct {
		Status string `json:"status"`
		Count  int    `json:"count"`
	}
	type failedTaskPreviewJSON struct {
		ID           string `json:"id"`
		Status       string `json:"status"`
		Process      string `json:"process"`
		Name         string `json:"name"`
		Tag          string `json:"tag"`
		Workdir      string `json:"workdir"`
		Exit         *int   `json:"exit"`
		ErrorSummary string `json:"error_summary"`
	}
	type logOnlyFailureJSON struct {
		ID           string `json:"id"`
		Workdir      string `json:"workdir"`
		Process      string `json:"process"`
		Name         string `json:"name"`
		Exit         *int   `json:"exit"`
		ErrorSummary string `json:"error_summary"`
		ErrorBlock   string `json:"error_block"`
	}
	type logOnlyEvidenceSourceJSON struct {
		Kind   string `json:"kind"`
		Path   string `json:"path"`
		Detail string `json:"detail"`
	}
	type logOnlyTaskEvidenceJSON struct {
		ID                    string                      `json:"id"`
		Workdir               string                      `json:"workdir"`
		Process               string                      `json:"process"`
		Name                  string                      `json:"name"`
		ObservedStatus        string                      `json:"observed_status"`
		Exit                  *int                        `json:"exit"`
		ErrorSummary          string                      `json:"error_summary"`
		ErrorBlock            string                      `json:"error_block"`
		Sources               []logOnlyEvidenceSourceJSON `json:"sources"`
		Completeness          string                      `json:"completeness"`
		CommandFilesAvailable bool                        `json:"command_files_available"`
	}
	type statusSummaryJSON struct {
		RunDir          string                    `json:"run_dir"`
		Mode            string                    `json:"mode"`
		IndexPath       string                    `json:"index_path"`
		Freshness       string                    `json:"freshness"`
		BuiltAt         *string                   `json:"built_at"`
		Sources         artifactSetJSON           `json:"sources"`
		Counts          []statusCountJSON         `json:"counts"`
		FailedCount     int                       `json:"failed_count"`
		FailedPreview   []failedTaskPreviewJSON   `json:"failed_preview"`
		LogOnlyFailures []logOnlyFailureJSON      `json:"log_only_failures"`
		LogOnlyEvidence []logOnlyTaskEvidenceJSON `json:"log_only_evidence"`
		Diagnostics     []diagnosticJSON          `json:"diagnostics"`
	}
	type statusViewJSON struct {
		Format  string            `json:"format"`
		Summary statusSummaryJSON `json:"summary"`
	}

	countsJSON := func(counts []domain.StatusCount) []statusCountJSON {
		ordered := append([]domain.StatusCount(nil), counts...)
		sort.SliceStable(ordered, func(i, j int) bool {
			left := string(ordered[i].Status)
			right := string(ordered[j].Status)
			if left == right {
				return ordered[i].Count < ordered[j].Count
			}
			return left < right
		})

		out := make([]statusCountJSON, len(ordered))
		for i, count := range ordered {
			out[i] = statusCountJSON{Status: string(count.Status), Count: count.Count}
		}
		return out
	}
	failedPreviewJSON := func(failures []domain.FailedTaskPreview) []failedTaskPreviewJSON {
		out := make([]failedTaskPreviewJSON, len(failures))
		for i, failure := range failures {
			out[i] = failedTaskPreviewJSON{
				ID:           failure.ID,
				Status:       string(failure.Status),
				Process:      failure.Process,
				Name:         failure.Name,
				Tag:          failure.Tag,
				Workdir:      failure.Workdir,
				Exit:         failure.Exit,
				ErrorSummary: failure.ErrorSummary,
			}
		}
		return out
	}
	logOnlyFailuresJSON := func(failures []domain.LogOnlyFailure) []logOnlyFailureJSON {
		out := make([]logOnlyFailureJSON, len(failures))
		for i, failure := range failures {
			out[i] = logOnlyFailureJSON{
				ID:           failure.ID,
				Workdir:      failure.Workdir,
				Process:      failure.Process,
				Name:         failure.Name,
				Exit:         failure.Exit,
				ErrorSummary: failure.ErrorSummary,
				ErrorBlock:   failure.ErrorBlock,
			}
		}
		return out
	}
	logOnlyEvidenceSourcesJSON := func(sources []domain.LogOnlyEvidenceSource) []logOnlyEvidenceSourceJSON {
		out := make([]logOnlyEvidenceSourceJSON, len(sources))
		for i, source := range sources {
			out[i] = logOnlyEvidenceSourceJSON{
				Kind:   string(source.Kind),
				Path:   source.Path,
				Detail: source.Detail,
			}
		}
		return out
	}
	logOnlyEvidenceJSON := func(evidence []domain.LogOnlyTaskEvidence) []logOnlyTaskEvidenceJSON {
		out := make([]logOnlyTaskEvidenceJSON, len(evidence))
		for i, item := range evidence {
			out[i] = logOnlyTaskEvidenceJSON{
				ID:                    item.ID,
				Workdir:               item.Workdir,
				Process:               item.Process,
				Name:                  item.Name,
				ObservedStatus:        string(item.ObservedStatus),
				Exit:                  item.Exit,
				ErrorSummary:          item.ErrorSummary,
				ErrorBlock:            item.ErrorBlock,
				Sources:               logOnlyEvidenceSourcesJSON(item.Sources),
				Completeness:          string(item.Completeness),
				CommandFilesAvailable: item.CommandFilesAvailable,
			}
		}
		return out
	}

	summary := view.Summary
	format := string(view.Format)
	if format == "" {
		format = string(domain.OutputFormatJSON)
	}
	payload := statusViewJSON{
		Format: format,
		Summary: statusSummaryJSON{
			RunDir:          summary.RunDir.Path,
			Mode:            string(summary.Mode),
			IndexPath:       summary.IndexPath,
			Freshness:       string(summary.Freshness),
			BuiltAt:         jsonTimePtr(summary.BuiltAt),
			Sources:         artifactSetJSONFor(summary.Sources),
			Counts:          countsJSON(summary.Counts),
			FailedCount:     summary.FailedCount,
			FailedPreview:   failedPreviewJSON(summary.FailedPreview),
			LogOnlyFailures: logOnlyFailuresJSON(summary.LogOnlyFailures),
			LogOnlyEvidence: logOnlyEvidenceJSON(summary.LogOnlyEvidence),
			Diagnostics:     diagnosticsJSON(summary.Diagnostics),
		},
	}

	return writeIndentedJSON(writer, payload, "render status json")
}

func RenderTasksHuman(writer io.Writer, view domain.TasksView) error {
	if writer == nil {
		return fmt.Errorf("render tasks human: nil writer")
	}

	inline := humanField
	exitText := formatOptionalInt
	nameTagText := func(task domain.Task) string {
		name := inline(task.Name)
		tag := inline(task.Tag)
		switch {
		case name == "-" && tag == "-":
			return "-"
		case name == "-":
			return tag
		case tag == "-":
			return name
		default:
			return name + "/" + tag
		}
	}

	hasErrorDiagnostic := func(diagnostics []domain.Diagnostic) bool {
		for _, diagnostic := range diagnostics {
			if diagnostic.Severity == domain.DiagnosticError {
				return true
			}
		}
		return false
	}
	appendDiagnostics := func(builder *strings.Builder, diagnostics []domain.Diagnostic) error {
		return appendDiagnosticBlocksHuman(builder, diagnostics, "render tasks human diagnostics")
	}

	var builder strings.Builder
	if len(view.Tasks) == 0 {
		if !hasErrorDiagnostic(view.Diagnostics) {
			builder.WriteString("tasks: none\n")
		}
	} else {
		builder.WriteString("id\tstatus\tprocess\tname/tag\tworkdir\texit\tduration\trealtime\tcpus\tmemory\n")
		for _, task := range view.Tasks {
			fmt.Fprintf(&builder, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				inline(task.ID),
				inline(string(task.Status)),
				inline(task.Process),
				nameTagText(task),
				inline(task.Workdir),
				exitText(task.Exit),
				inline(task.Duration),
				inline(task.Realtime),
				inline(task.CPUs),
				inline(task.Memory),
			)
		}
	}

	if err := appendDiagnostics(&builder, view.Diagnostics); err != nil {
		return err
	}

	if _, err := io.WriteString(writer, builder.String()); err != nil {
		return fmt.Errorf("render tasks human: write: %w", err)
	}
	return nil
}

func RenderTasksJSON(writer io.Writer, view domain.TasksView) error {
	if writer == nil {
		return fmt.Errorf("render tasks json: nil writer")
	}

	type taskQueryJSON struct {
		ProcessSubstring string `json:"process_substring"`
		NameSubstring    string `json:"name_substring"`
		SampleSubstring  string `json:"sample_substring"`
		Status           string `json:"status"`
		StatusRaw        string `json:"status_raw"`
	}
	type tasksViewJSON struct {
		Format      string             `json:"format"`
		Query       taskQueryJSON      `json:"query"`
		Metadata    *indexMetadataJSON `json:"metadata"`
		Tasks       []taskJSON         `json:"tasks"`
		Diagnostics []diagnosticJSON   `json:"diagnostics"`
	}

	format := string(view.Format)
	if format == "" {
		format = string(domain.OutputFormatJSON)
	}
	payload := tasksViewJSON{
		Format: format,
		Query: taskQueryJSON{
			ProcessSubstring: view.Query.ProcessSubstring,
			NameSubstring:    view.Query.NameSubstring,
			SampleSubstring:  view.Query.SampleSubstring,
			Status:           string(view.Query.Status),
			StatusRaw:        view.Query.StatusRaw,
		},
		Metadata:    indexMetadataJSONFor(view.Metadata),
		Tasks:       tasksJSON(view.Tasks),
		Diagnostics: diagnosticsJSON(view.Diagnostics),
	}

	return writeIndentedJSON(writer, payload, "render tasks json")
}

func RenderInspectHuman(writer io.Writer, view domain.InspectView) error {
	if writer == nil {
		return fmt.Errorf("render inspect human: nil writer")
	}
	if view.EvidenceKind == domain.InspectEvidenceLogOnly || view.LogOnlyResolution != nil || view.LogOnlyDossier != nil {
		return RenderLogOnlyInspectHuman(writer, view)
	}

	inline := humanField
	exitText := formatOptionalInt
	nameTagText := func(task domain.Task) string {
		name := inline(task.Name)
		tag := inline(task.Tag)
		switch {
		case name == "-" && tag == "-":
			return "-"
		case name == "-":
			return tag
		case tag == "-":
			return name
		default:
			return name + "/" + tag
		}
	}
	appendDiagnostics := func(builder *strings.Builder, diagnostics []domain.Diagnostic) error {
		return appendDiagnosticBlocksHuman(builder, diagnostics, "render inspect human diagnostics")
	}
	resolution := view.Resolution
	var builder strings.Builder
	fmt.Fprintf(&builder, "selector: %s\n", inline(resolution.Selector))
	fmt.Fprintf(&builder, "resolution: %s\n", inline(string(resolution.Kind)))

	switch resolution.Kind {
	case domain.SelectorResolutionExact:
		var task domain.Task
		var files []domain.CommandFile
		var dossierDiagnostics []domain.Diagnostic
		hasTask := false
		if view.Dossier != nil {
			task = view.Dossier.Task
			files = view.Dossier.Inventory.Files
			dossierDiagnostics = view.Dossier.Diagnostics
			hasTask = true
		} else if resolution.Task != nil {
			task = *resolution.Task
			hasTask = true
		}

		if hasTask {
			builder.WriteString("task:\n")
			fmt.Fprintf(&builder, "  id: %s\n", inline(task.ID))
			fmt.Fprintf(&builder, "  status: %s\n", inline(string(task.Status)))
			fmt.Fprintf(&builder, "  process: %s\n", inline(task.Process))
			fmt.Fprintf(&builder, "  name: %s\n", inline(task.Name))
			fmt.Fprintf(&builder, "  tag: %s\n", inline(task.Tag))
			fmt.Fprintf(&builder, "  workdir: %s\n", inline(task.Workdir))
			fmt.Fprintf(&builder, "  exit: %s\n", exitText(task.Exit))
			fmt.Fprintf(&builder, "  duration: %s\n", inline(task.Duration))
			fmt.Fprintf(&builder, "  realtime: %s\n", inline(task.Realtime))
			fmt.Fprintf(&builder, "  cpus: %s\n", inline(task.CPUs))
			fmt.Fprintf(&builder, "  memory: %s\n", inline(task.Memory))
			fmt.Fprintf(&builder, "  error_summary: %s\n", inline(task.ErrorSummary))
			writeCommandFilesHuman(&builder, files, "      ")
		} else {
			builder.WriteString("dossier: unavailable\n")
		}
		if err := appendDiagnostics(&builder, collectDiagnostics(dossierDiagnostics, resolution.Diagnostics, view.Diagnostics)); err != nil {
			return err
		}
	case domain.SelectorResolutionAmbiguous:
		if len(resolution.Matches) == 0 {
			builder.WriteString("matches: none\n")
		} else {
			builder.WriteString("matches:\n")
			builder.WriteString("id\tstatus\tprocess\tname/tag\tworkdir\n")
			for _, task := range resolution.Matches {
				fmt.Fprintf(&builder, "%s\t%s\t%s\t%s\t%s\n",
					inline(task.ID),
					inline(string(task.Status)),
					inline(task.Process),
					nameTagText(task),
					inline(task.Workdir),
				)
			}
		}
		if err := appendDiagnostics(&builder, collectDiagnostics(resolution.Diagnostics, view.Diagnostics)); err != nil {
			return err
		}
	case domain.SelectorResolutionNotFound:
		builder.WriteString("matches: none\n")
		if err := appendDiagnostics(&builder, collectDiagnostics(resolution.Diagnostics, view.Diagnostics)); err != nil {
			return err
		}
	default:
		if err := appendDiagnostics(&builder, collectDiagnostics(resolution.Diagnostics, view.Diagnostics)); err != nil {
			return err
		}
	}

	if _, err := io.WriteString(writer, builder.String()); err != nil {
		return fmt.Errorf("render inspect human: write: %w", err)
	}
	return nil
}

func RenderInspectJSON(writer io.Writer, view domain.InspectView) error {
	if writer == nil {
		return fmt.Errorf("render inspect json: nil writer")
	}

	type snippetJSON struct {
		Path      string `json:"path"`
		Strategy  string `json:"strategy"`
		StartLine int    `json:"start_line"`
		EndLine   int    `json:"end_line"`
		Content   string `json:"content"`
		Truncated bool   `json:"truncated"`
		MaxBytes  int64  `json:"max_bytes"`
	}
	type commandFileJSON struct {
		Kind    string       `json:"kind"`
		Path    string       `json:"path"`
		Exists  bool         `json:"exists"`
		Size    int64        `json:"size"`
		Snippet *snippetJSON `json:"snippet"`
	}
	type commandFileInventoryJSON struct {
		Workdir string            `json:"workdir"`
		Files   []commandFileJSON `json:"files"`
	}
	type logOnlyEvidenceSourceJSON struct {
		Kind   string `json:"kind"`
		Path   string `json:"path"`
		Detail string `json:"detail"`
	}
	type logOnlyTaskEvidenceJSON struct {
		ID                    string                      `json:"id"`
		Workdir               string                      `json:"workdir"`
		Process               string                      `json:"process"`
		Name                  string                      `json:"name"`
		ObservedStatus        string                      `json:"observed_status"`
		Exit                  *int                        `json:"exit"`
		ErrorSummary          string                      `json:"error_summary"`
		ErrorBlock            string                      `json:"error_block"`
		Sources               []logOnlyEvidenceSourceJSON `json:"sources"`
		Completeness          string                      `json:"completeness"`
		CommandFilesAvailable bool                        `json:"command_files_available"`
	}
	type taskDossierJSON struct {
		Task        taskJSON                 `json:"task"`
		Inventory   commandFileInventoryJSON `json:"inventory"`
		Diagnostics []diagnosticJSON         `json:"diagnostics"`
	}
	type selectorResolutionJSON struct {
		Kind        string           `json:"kind"`
		Selector    string           `json:"selector"`
		Task        *taskJSON        `json:"task"`
		Matches     []taskJSON       `json:"matches"`
		Diagnostics []diagnosticJSON `json:"diagnostics"`
	}
	type logOnlySelectorResolutionJSON struct {
		Kind        string                    `json:"kind"`
		Selector    string                    `json:"selector"`
		Evidence    *logOnlyTaskEvidenceJSON  `json:"evidence"`
		Matches     []logOnlyTaskEvidenceJSON `json:"matches"`
		Diagnostics []diagnosticJSON          `json:"diagnostics"`
	}
	type logOnlyTaskDossierJSON struct {
		Evidence    logOnlyTaskEvidenceJSON  `json:"evidence"`
		Inventory   commandFileInventoryJSON `json:"inventory"`
		Diagnostics []diagnosticJSON         `json:"diagnostics"`
	}
	type inspectViewJSON struct {
		Format            string                 `json:"format"`
		EvidenceKind      string                 `json:"evidence_kind,omitempty"`
		Resolution        selectorResolutionJSON `json:"resolution"`
		Dossier           *taskDossierJSON       `json:"dossier"`
		LogOnlyResolution any                    `json:"log_only_resolution,omitempty"`
		LogOnlyDossier    any                    `json:"log_only_dossier,omitempty"`
		Diagnostics       []diagnosticJSON       `json:"diagnostics"`
	}

	snippetJSONFor := func(snippet *domain.Snippet) *snippetJSON {
		if snippet == nil {
			return nil
		}
		return &snippetJSON{
			Path:      snippet.Path,
			Strategy:  string(snippet.Strategy),
			StartLine: snippet.StartLine,
			EndLine:   snippet.EndLine,
			Content:   snippet.Content,
			Truncated: snippet.Truncated,
			MaxBytes:  snippet.MaxBytes,
		}
	}
	commandFilesJSON := func(files []domain.CommandFile) []commandFileJSON {
		ordered := orderedCommandFiles(files)

		out := make([]commandFileJSON, len(ordered))
		for i, file := range ordered {
			out[i] = commandFileJSON{
				Kind:    string(file.Kind),
				Path:    file.Path,
				Exists:  file.Exists,
				Size:    file.Size,
				Snippet: snippetJSONFor(file.Snippet),
			}
		}
		return out
	}
	dossierJSON := func(dossier *domain.TaskDossier) *taskDossierJSON {
		if dossier == nil {
			return nil
		}
		return &taskDossierJSON{
			Task: taskValueJSON(dossier.Task),
			Inventory: commandFileInventoryJSON{
				Workdir: dossier.Inventory.Workdir,
				Files:   commandFilesJSON(dossier.Inventory.Files),
			},
			Diagnostics: diagnosticsJSON(dossier.Diagnostics),
		}
	}
	logOnlyEvidenceSourcesJSON := func(sources []domain.LogOnlyEvidenceSource) []logOnlyEvidenceSourceJSON {
		out := make([]logOnlyEvidenceSourceJSON, len(sources))
		for i, source := range sources {
			out[i] = logOnlyEvidenceSourceJSON{
				Kind:   string(source.Kind),
				Path:   source.Path,
				Detail: source.Detail,
			}
		}
		return out
	}
	logOnlyEvidenceJSON := func(evidence domain.LogOnlyTaskEvidence) logOnlyTaskEvidenceJSON {
		return logOnlyTaskEvidenceJSON{
			ID:                    evidence.ID,
			Workdir:               evidence.Workdir,
			Process:               evidence.Process,
			Name:                  evidence.Name,
			ObservedStatus:        string(evidence.ObservedStatus),
			Exit:                  evidence.Exit,
			ErrorSummary:          evidence.ErrorSummary,
			ErrorBlock:            evidence.ErrorBlock,
			Sources:               logOnlyEvidenceSourcesJSON(evidence.Sources),
			Completeness:          string(evidence.Completeness),
			CommandFilesAvailable: evidence.CommandFilesAvailable,
		}
	}
	logOnlyEvidencePointerJSON := func(evidence *domain.LogOnlyTaskEvidence) *logOnlyTaskEvidenceJSON {
		if evidence == nil {
			return nil
		}
		out := logOnlyEvidenceJSON(*evidence)
		return &out
	}
	logOnlyEvidenceSliceJSON := func(evidence []domain.LogOnlyTaskEvidence) []logOnlyTaskEvidenceJSON {
		out := make([]logOnlyTaskEvidenceJSON, len(evidence))
		for i, item := range evidence {
			out[i] = logOnlyEvidenceJSON(item)
		}
		return out
	}
	logOnlyResolutionJSON := func(resolution *domain.LogOnlySelectorResolution) *logOnlySelectorResolutionJSON {
		if resolution == nil {
			return nil
		}
		return &logOnlySelectorResolutionJSON{
			Kind:        string(resolution.Kind),
			Selector:    resolution.Selector,
			Evidence:    logOnlyEvidencePointerJSON(resolution.Evidence),
			Matches:     logOnlyEvidenceSliceJSON(resolution.Matches),
			Diagnostics: diagnosticsJSON(resolution.Diagnostics),
		}
	}
	logOnlyDossierJSON := func(dossier *domain.LogOnlyTaskDossier) *logOnlyTaskDossierJSON {
		if dossier == nil {
			return nil
		}
		return &logOnlyTaskDossierJSON{
			Evidence: logOnlyEvidenceJSON(dossier.Evidence),
			Inventory: commandFileInventoryJSON{
				Workdir: dossier.Inventory.Workdir,
				Files:   commandFilesJSON(dossier.Inventory.Files),
			},
			Diagnostics: diagnosticsJSON(dossier.Diagnostics),
		}
	}

	format := string(view.Format)
	if format == "" {
		format = string(domain.OutputFormatJSON)
	}
	hasLogOnlyInspect := view.EvidenceKind == domain.InspectEvidenceLogOnly || view.LogOnlyResolution != nil || view.LogOnlyDossier != nil
	evidenceKind := string(view.EvidenceKind)
	if evidenceKind == "" && hasLogOnlyInspect {
		evidenceKind = string(domain.InspectEvidenceLogOnly)
	}

	payload := inspectViewJSON{
		Format:       format,
		EvidenceKind: evidenceKind,
		Resolution: selectorResolutionJSON{
			Kind:        string(view.Resolution.Kind),
			Selector:    view.Resolution.Selector,
			Task:        taskPointerJSON(view.Resolution.Task),
			Matches:     tasksJSON(view.Resolution.Matches),
			Diagnostics: diagnosticsJSON(view.Resolution.Diagnostics),
		},
		Dossier:     dossierJSON(view.Dossier),
		Diagnostics: diagnosticsJSON(view.Diagnostics),
	}
	if hasLogOnlyInspect {
		if view.LogOnlyResolution != nil {
			payload.LogOnlyResolution = logOnlyResolutionJSON(view.LogOnlyResolution)
		} else {
			payload.LogOnlyResolution = (*logOnlySelectorResolutionJSON)(nil)
		}
		if view.LogOnlyDossier != nil {
			payload.LogOnlyDossier = logOnlyDossierJSON(view.LogOnlyDossier)
		} else {
			payload.LogOnlyDossier = (*logOnlyTaskDossierJSON)(nil)
		}
	}

	return writeIndentedJSON(writer, payload, "render inspect json")
}

func RenderIndexHuman(writer io.Writer, view domain.IndexView) error {
	if writer == nil {
		return fmt.Errorf("render index human: nil writer")
	}

	inline := humanField
	timestamp := func(value time.Time) string {
		if value.IsZero() {
			return "-"
		}
		return value.UTC().Format(time.RFC3339)
	}
	cleanPatterns := func(values []string) []string {
		patterns := make([]string, 0, len(values))
		for _, value := range values {
			value = strings.TrimSpace(value)
			if value == "" {
				continue
			}
			patterns = append(patterns, value)
		}
		return patterns
	}
	appendDiagnostics := func(builder *strings.Builder, diagnostics []domain.Diagnostic) error {
		return appendDiagnosticBlocksHuman(builder, diagnostics, "render index human diagnostics")
	}
	writeSearchedLocations := func(builder *strings.Builder, locations []domain.ArtifactSearchLocation, fallbackPatterns []string) {
		wroteLocations := false
		for _, location := range locations {
			patterns := cleanPatterns(location.Patterns)
			kind := inline(string(location.Kind))
			baseDir := strings.TrimSpace(location.BaseDir)
			description := inline(location.Description)
			if description == "-" {
				description = ""
			}
			if kind == "-" && baseDir == "" && description == "" && len(patterns) == 0 {
				continue
			}

			if !wroteLocations {
				builder.WriteString("  searched_locations:\n")
				wroteLocations = true
			}

			if baseDir == "" {
				baseDir = "-"
			}
			patternText := "none"
			if len(patterns) > 0 {
				patternText = strings.Join(patterns, ", ")
			}
			if description == "" {
				fmt.Fprintf(builder, "    - %s: %s (patterns: %s)\n", kind, baseDir, patternText)
				continue
			}
			fmt.Fprintf(builder, "    - %s: %s (%s; patterns: %s)\n", kind, baseDir, description, patternText)
		}
		if wroteLocations {
			return
		}

		patterns := cleanPatterns(fallbackPatterns)
		if len(patterns) > 0 {
			fmt.Fprintf(builder, "  searched_patterns: %s\n", strings.Join(patterns, ", "))
		}
	}

	diagnostics := view.Diagnostics
	metadata := diagnostics.Metadata
	runDir := diagnostics.RunDir.Path
	if strings.TrimSpace(runDir) == "" {
		runDir = diagnostics.Artifacts.RunDir.Path
	}
	if strings.TrimSpace(runDir) == "" && metadata != nil {
		runDir = metadata.RunDir
	}

	mode := diagnostics.Artifacts.Mode
	if mode == "" && metadata != nil {
		mode = metadata.Mode
	}

	freshness := domain.IndexFreshnessUnknown
	if metadata != nil && metadata.Freshness != "" {
		freshness = metadata.Freshness
	} else if mode == domain.IndexModeUnsupported {
		freshness = domain.IndexFreshnessUnsupported
	}

	useMetadataSources := metadata != nil && diagnostics.Artifacts.Mode == "" && diagnostics.Artifacts.Trace == nil && diagnostics.Artifacts.Log == nil && strings.TrimSpace(diagnostics.Artifacts.RunDir.Path) == "" && diagnostics.Artifacts.SelectedAt.IsZero() && len(diagnostics.Artifacts.SearchedPatterns) == 0 && len(diagnostics.Artifacts.Diagnostics) == 0
	trace := diagnostics.Artifacts.Trace
	log := diagnostics.Artifacts.Log
	if useMetadataSources {
		trace = metadata.Trace
		log = metadata.Log
	}

	var builder strings.Builder
	fmt.Fprintf(&builder, "run_dir: %s\n", inline(runDir))
	fmt.Fprintf(&builder, "mode: %s\n", inline(string(mode)))
	if metadata != nil && strings.TrimSpace(metadata.IndexPath) != "" {
		fmt.Fprintf(&builder, "index: %s\n", strings.TrimSpace(metadata.IndexPath))
	}
	fmt.Fprintf(&builder, "freshness: %s\n", inline(string(freshness)))
	if metadata != nil && strings.TrimSpace(metadata.StaleReason) != "" {
		fmt.Fprintf(&builder, "stale_reason: %s\n", inline(metadata.StaleReason))
	}
	if metadata != nil && !metadata.BuiltAt.IsZero() {
		fmt.Fprintf(&builder, "built_at: %s\n", timestamp(metadata.BuiltAt))
	}
	if metadata == nil {
		builder.WriteString("task_count: unknown\n")
	} else {
		fmt.Fprintf(&builder, "task_count: %d\n", metadata.TaskCount)
	}

	writeSource := func(label string, source *domain.SourceFingerprint) {
		if source == nil || strings.TrimSpace(source.Path) == "" {
			fmt.Fprintf(&builder, "  %s: none\n", label)
			return
		}
		fmt.Fprintf(&builder, "  %s: %s (mtime=%s size=%d)\n", label, strings.TrimSpace(source.Path), timestamp(source.ModTime), source.Size)
	}

	builder.WriteString("sources:\n")
	writeSource("trace", trace)
	writeSource("log", log)
	writeSearchedLocations(&builder, diagnostics.Artifacts.SearchLocations, diagnostics.Artifacts.SearchedPatterns)

	if err := appendDiagnostics(&builder, collectDiagnostics(diagnostics.Artifacts.Diagnostics, diagnostics.Diagnostics)); err != nil {
		return err
	}

	if _, err := io.WriteString(writer, builder.String()); err != nil {
		return fmt.Errorf("render index human: write: %w", err)
	}
	return nil
}

func RenderIndexJSON(writer io.Writer, view domain.IndexView) error {
	if writer == nil {
		return fmt.Errorf("render index json: nil writer")
	}

	type indexViewJSON struct {
		Format      string             `json:"format"`
		RunDir      string             `json:"run_dir"`
		Artifacts   artifactSetJSON    `json:"artifacts"`
		Metadata    *indexMetadataJSON `json:"metadata"`
		Diagnostics []diagnosticJSON   `json:"diagnostics"`
	}

	format := string(view.Format)
	if format == "" {
		format = string(domain.OutputFormatJSON)
	}
	diagnostics := view.Diagnostics
	payload := indexViewJSON{
		Format:      format,
		RunDir:      diagnostics.RunDir.Path,
		Artifacts:   artifactSetJSONFor(diagnostics.Artifacts),
		Metadata:    indexMetadataJSONFor(diagnostics.Metadata),
		Diagnostics: diagnosticsJSON(diagnostics.Diagnostics),
	}

	return writeIndentedJSON(writer, payload, "render index json")
}

func RenderUnsupportedDiagnostics(writer io.Writer, diagnostics []domain.Diagnostic, format domain.OutputFormat) error {
	renderHuman := func() error {
		blocks := BuildDiagnosticBlocks(diagnostics)
		if err := RenderDiagnosticBlocksHuman(writer, blocks); err != nil {
			return fmt.Errorf("render unsupported diagnostics human: %w", err)
		}
		return nil
	}

	renderJSON := func() error {
		if writer == nil {
			return fmt.Errorf("render unsupported diagnostics json: nil writer")
		}

		type unsupportedDiagnosticsJSON struct {
			Format      string           `json:"format"`
			Diagnostics []diagnosticJSON `json:"diagnostics"`
		}

		payload := unsupportedDiagnosticsJSON{
			Format:      string(domain.OutputFormatJSON),
			Diagnostics: diagnosticsJSON(diagnostics),
		}

		return writeIndentedJSON(writer, payload, "render unsupported diagnostics json")
	}

	switch format {
	case domain.OutputFormatJSON:
		return renderJSON()
	case "", domain.OutputFormatHuman:
		return renderHuman()
	default:
		return fmt.Errorf("render unsupported diagnostics: unsupported format %q", format)
	}
}
