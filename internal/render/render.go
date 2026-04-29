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

func writeIndentedJSON(writer io.Writer, payload any, context string) error {
	encoder := json.NewEncoder(writer)
	encoder.SetEscapeHTML(false)
	encoder.SetIndent("", "  ")
	if err := encoder.Encode(payload); err != nil {
		return fmt.Errorf("%s: write: %w", context, err)
	}
	return nil
}

func RenderStatusHuman(writer io.Writer, view domain.StatusView) error {
	if writer == nil {
		return fmt.Errorf("render status human: nil writer")
	}

	summary := view.Summary
	inline := func(value string) string {
		value = strings.ReplaceAll(value, "\r\n", "\n")
		value = strings.ReplaceAll(value, "\r", "\n")
		fields := strings.Fields(value)
		if len(fields) == 0 {
			return "-"
		}
		return strings.Join(fields, " ")
	}
	timestamp := func(value time.Time) string {
		if value.IsZero() {
			return "-"
		}
		return value.UTC().Format(time.RFC3339)
	}
	exitText := func(value *int) string {
		if value == nil {
			return "-"
		}
		return fmt.Sprintf("%d", *value)
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

	if summary.Mode == domain.IndexModeLogOnly {
		builder.WriteString("counts: unavailable (log-only mode; trace file required)\n")
	} else {
		counts := append([]domain.StatusCount(nil), summary.Counts...)
		sort.SliceStable(counts, func(i, j int) bool {
			left := string(counts[i].Status)
			right := string(counts[j].Status)
			if left == right {
				return counts[i].Count < counts[j].Count
			}
			return left < right
		})

		if len(counts) == 0 {
			builder.WriteString("counts: none\n")
		} else {
			builder.WriteString("counts:\n")
			for _, count := range counts {
				fmt.Fprintf(&builder, "  %s: %d\n", inline(string(count.Status)), count.Count)
			}
		}
	}

	fmt.Fprintf(&builder, "failed_count: %d\n", summary.FailedCount)
	if summary.Mode == domain.IndexModeLogOnly {
		if len(summary.LogOnlyFailures) == 0 {
			builder.WriteString("log_only_failures: none\n")
		} else {
			builder.WriteString("log_only_failures:\n")
			for _, failure := range summary.LogOnlyFailures {
				fmt.Fprintf(&builder, "  - id=%s process=%s name=%s workdir=%s exit=%s error=%s\n",
					inline(failure.ID),
					inline(failure.Process),
					inline(failure.Name),
					inline(failure.Workdir),
					exitText(failure.Exit),
					inline(failure.ErrorSummary),
				)
			}
		}
	} else {
		if len(summary.FailedPreview) == 0 {
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
	}

	if len(summary.Diagnostics) > 0 {
		builder.WriteString("diagnostics:\n")
		for _, diagnostic := range summary.Diagnostics {
			severity := inline(string(diagnostic.Severity))
			code := inline(diagnostic.Code)
			message := inline(diagnostic.Message)
			if code == "-" {
				fmt.Fprintf(&builder, "  - %s: %s\n", severity, message)
			} else {
				fmt.Fprintf(&builder, "  - %s %s: %s\n", severity, code, message)
			}

			detail := strings.ReplaceAll(diagnostic.Detail, "\r\n", "\n")
			detail = strings.ReplaceAll(detail, "\r", "\n")
			for _, line := range strings.Split(detail, "\n") {
				if strings.TrimSpace(line) == "" {
					continue
				}
				fmt.Fprintf(&builder, "    detail: %s\n", inline(line))
			}
		}
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
	type statusSummaryJSON struct {
		RunDir          string                  `json:"run_dir"`
		Mode            string                  `json:"mode"`
		IndexPath       string                  `json:"index_path"`
		Freshness       string                  `json:"freshness"`
		BuiltAt         *string                 `json:"built_at"`
		Sources         artifactSetJSON         `json:"sources"`
		Counts          []statusCountJSON       `json:"counts"`
		FailedCount     int                     `json:"failed_count"`
		FailedPreview   []failedTaskPreviewJSON `json:"failed_preview"`
		LogOnlyFailures []logOnlyFailureJSON    `json:"log_only_failures"`
		Diagnostics     []diagnosticJSON        `json:"diagnostics"`
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
			Diagnostics:     diagnosticsJSON(summary.Diagnostics),
		},
	}

	return writeIndentedJSON(writer, payload, "render status json")
}

func RenderTasksHuman(writer io.Writer, view domain.TasksView) error {
	if writer == nil {
		return fmt.Errorf("render tasks human: nil writer")
	}

	inline := func(value string) string {
		value = strings.ReplaceAll(value, "\r\n", "\n")
		value = strings.ReplaceAll(value, "\r", "\n")
		fields := strings.Fields(value)
		if len(fields) == 0 {
			return "-"
		}
		return strings.Join(fields, " ")
	}
	exitText := func(value *int) string {
		if value == nil {
			return "-"
		}
		return fmt.Sprintf("%d", *value)
	}
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

	var builder strings.Builder
	if len(view.Tasks) == 0 {
		builder.WriteString("tasks: none\n")
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

	if len(view.Diagnostics) > 0 {
		builder.WriteString("diagnostics:\n")
		for _, diagnostic := range view.Diagnostics {
			severity := inline(string(diagnostic.Severity))
			code := inline(diagnostic.Code)
			message := inline(diagnostic.Message)
			if code == "-" {
				fmt.Fprintf(&builder, "  - %s: %s\n", severity, message)
			} else {
				fmt.Fprintf(&builder, "  - %s %s: %s\n", severity, code, message)
			}

			detail := strings.ReplaceAll(diagnostic.Detail, "\r\n", "\n")
			detail = strings.ReplaceAll(detail, "\r", "\n")
			for _, line := range strings.Split(detail, "\n") {
				if strings.TrimSpace(line) == "" {
					continue
				}
				fmt.Fprintf(&builder, "    detail: %s\n", inline(line))
			}
		}
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

	inline := func(value string) string {
		value = strings.ReplaceAll(value, "\r\n", "\n")
		value = strings.ReplaceAll(value, "\r", "\n")
		fields := strings.Fields(value)
		if len(fields) == 0 {
			return "-"
		}
		return strings.Join(fields, " ")
	}
	exitText := func(value *int) string {
		if value == nil {
			return "-"
		}
		return fmt.Sprintf("%d", *value)
	}
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
	lineRangeText := func(snippet *domain.Snippet) string {
		if snippet == nil || snippet.StartLine <= 0 || snippet.EndLine <= 0 {
			return "-"
		}
		return fmt.Sprintf("%d-%d", snippet.StartLine, snippet.EndLine)
	}
	commandFileRank := func(kind domain.CommandFileKind) int {
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
	writeDiagnostics := func(builder *strings.Builder, diagnostics []domain.Diagnostic) {
		if len(diagnostics) == 0 {
			return
		}
		builder.WriteString("diagnostics:\n")
		for _, diagnostic := range diagnostics {
			severity := inline(string(diagnostic.Severity))
			code := inline(diagnostic.Code)
			message := inline(diagnostic.Message)
			if code == "-" {
				fmt.Fprintf(builder, "  - %s: %s\n", severity, message)
			} else {
				fmt.Fprintf(builder, "  - %s %s: %s\n", severity, code, message)
			}

			detail := strings.ReplaceAll(diagnostic.Detail, "\r\n", "\n")
			detail = strings.ReplaceAll(detail, "\r", "\n")
			for _, line := range strings.Split(detail, "\n") {
				if strings.TrimSpace(line) == "" {
					continue
				}
				fmt.Fprintf(builder, "    detail: %s\n", inline(line))
			}
		}
	}
	collectDiagnostics := func(parts ...[]domain.Diagnostic) []domain.Diagnostic {
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
	writeSnippetContent := func(builder *strings.Builder, content string) {
		content = strings.ReplaceAll(content, "\r\n", "\n")
		content = strings.ReplaceAll(content, "\r", "\n")
		content = strings.TrimSuffix(content, "\n")
		if content == "" {
			builder.WriteString("      -\n")
			return
		}
		for _, line := range strings.Split(content, "\n") {
			fmt.Fprintf(builder, "      %s\n", line)
		}
	}
	writeCommandFiles := func(builder *strings.Builder, files []domain.CommandFile) {
		if len(files) == 0 {
			builder.WriteString("command_files: none\n")
			return
		}

		ordered := append([]domain.CommandFile(nil), files...)
		sort.SliceStable(ordered, func(i, j int) bool {
			leftRank := commandFileRank(ordered[i].Kind)
			rightRank := commandFileRank(ordered[j].Kind)
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

		builder.WriteString("command_files:\n")
		for _, file := range ordered {
			path := file.Path
			if strings.TrimSpace(path) == "" && file.Snippet != nil {
				path = file.Snippet.Path
			}
			fmt.Fprintf(builder, "  - kind=%s path=%s exists=%t size=%d\n",
				inline(string(file.Kind)),
				inline(path),
				file.Exists,
				file.Size,
			)
			if file.Snippet == nil {
				continue
			}
			fmt.Fprintf(builder, "    snippet: strategy=%s lines=%s truncated=%t max_bytes=%d\n",
				inline(string(file.Snippet.Strategy)),
				lineRangeText(file.Snippet),
				file.Snippet.Truncated,
				file.Snippet.MaxBytes,
			)
			builder.WriteString("    content:\n")
			writeSnippetContent(builder, file.Snippet.Content)
		}
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
			writeCommandFiles(&builder, files)
		} else {
			builder.WriteString("dossier: unavailable\n")
		}
		writeDiagnostics(&builder, collectDiagnostics(dossierDiagnostics, resolution.Diagnostics, view.Diagnostics))
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
		writeDiagnostics(&builder, collectDiagnostics(resolution.Diagnostics, view.Diagnostics))
	case domain.SelectorResolutionNotFound:
		builder.WriteString("matches: none\n")
		writeDiagnostics(&builder, collectDiagnostics(resolution.Diagnostics, view.Diagnostics))
	default:
		writeDiagnostics(&builder, collectDiagnostics(resolution.Diagnostics, view.Diagnostics))
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
	type inspectViewJSON struct {
		Format      string                 `json:"format"`
		Resolution  selectorResolutionJSON `json:"resolution"`
		Dossier     *taskDossierJSON       `json:"dossier"`
		Diagnostics []diagnosticJSON       `json:"diagnostics"`
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
	commandFileRank := func(kind domain.CommandFileKind) int {
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
	commandFilesJSON := func(files []domain.CommandFile) []commandFileJSON {
		ordered := append([]domain.CommandFile(nil), files...)
		sort.SliceStable(ordered, func(i, j int) bool {
			leftRank := commandFileRank(ordered[i].Kind)
			rightRank := commandFileRank(ordered[j].Kind)
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

	format := string(view.Format)
	if format == "" {
		format = string(domain.OutputFormatJSON)
	}
	payload := inspectViewJSON{
		Format: format,
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

	return writeIndentedJSON(writer, payload, "render inspect json")
}

func RenderIndexHuman(writer io.Writer, view domain.IndexView) error {
	if writer == nil {
		return fmt.Errorf("render index human: nil writer")
	}

	inline := func(value string) string {
		value = strings.ReplaceAll(value, "\r\n", "\n")
		value = strings.ReplaceAll(value, "\r", "\n")
		fields := strings.Fields(value)
		if len(fields) == 0 {
			return "-"
		}
		return strings.Join(fields, " ")
	}
	timestamp := func(value time.Time) string {
		if value.IsZero() {
			return "-"
		}
		return value.UTC().Format(time.RFC3339)
	}
	writeDiagnostics := func(builder *strings.Builder, diagnostics []domain.Diagnostic) {
		if len(diagnostics) == 0 {
			return
		}
		builder.WriteString("diagnostics:\n")
		for _, diagnostic := range diagnostics {
			severity := inline(string(diagnostic.Severity))
			code := inline(diagnostic.Code)
			message := inline(diagnostic.Message)
			if code == "-" {
				fmt.Fprintf(builder, "  - %s: %s\n", severity, message)
			} else {
				fmt.Fprintf(builder, "  - %s %s: %s\n", severity, code, message)
			}

			detail := strings.ReplaceAll(diagnostic.Detail, "\r\n", "\n")
			detail = strings.ReplaceAll(detail, "\r", "\n")
			for _, line := range strings.Split(detail, "\n") {
				if strings.TrimSpace(line) == "" {
					continue
				}
				fmt.Fprintf(builder, "    detail: %s\n", inline(line))
			}
		}
	}
	collectDiagnostics := func(parts ...[]domain.Diagnostic) []domain.Diagnostic {
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
	if len(diagnostics.Artifacts.SearchedPatterns) > 0 {
		patterns := make([]string, 0, len(diagnostics.Artifacts.SearchedPatterns))
		for _, pattern := range diagnostics.Artifacts.SearchedPatterns {
			if strings.TrimSpace(pattern) != "" {
				patterns = append(patterns, strings.TrimSpace(pattern))
			}
		}
		if len(patterns) > 0 {
			fmt.Fprintf(&builder, "  searched_patterns: %s\n", strings.Join(patterns, ", "))
		}
	}

	writeDiagnostics(&builder, collectDiagnostics(diagnostics.Artifacts.Diagnostics, diagnostics.Diagnostics))

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
		if writer == nil {
			return fmt.Errorf("render unsupported diagnostics human: nil writer")
		}

		inline := func(value string) string {
			value = strings.ReplaceAll(value, "\r\n", "\n")
			value = strings.ReplaceAll(value, "\r", "\n")
			fields := strings.Fields(value)
			if len(fields) == 0 {
				return "-"
			}
			return strings.Join(fields, " ")
		}

		var builder strings.Builder
		if len(diagnostics) == 0 {
			builder.WriteString("diagnostics: none\n")
		} else {
			builder.WriteString("diagnostics:\n")
			for _, diagnostic := range diagnostics {
				severity := inline(string(diagnostic.Severity))
				code := inline(diagnostic.Code)
				message := inline(diagnostic.Message)
				if code == "-" {
					fmt.Fprintf(&builder, "  - %s: %s\n", severity, message)
				} else {
					fmt.Fprintf(&builder, "  - %s %s: %s\n", severity, code, message)
				}

				detail := strings.ReplaceAll(diagnostic.Detail, "\r\n", "\n")
				detail = strings.ReplaceAll(detail, "\r", "\n")
				for _, line := range strings.Split(detail, "\n") {
					if strings.TrimSpace(line) == "" {
						continue
					}
					fmt.Fprintf(&builder, "    detail: %s\n", inline(line))
				}
			}
		}

		if _, err := io.WriteString(writer, builder.String()); err != nil {
			return fmt.Errorf("render unsupported diagnostics human: write: %w", err)
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
