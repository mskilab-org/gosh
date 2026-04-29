package tasks

import (
	"fmt"
	"path/filepath"
	"sort"
	"strings"

	"github.com/mskilab-org/gosh/internal/domain"
	"github.com/mskilab-org/gosh/internal/trace"
)

func NormalizeTaskQuery(query domain.TaskQuery) (domain.TaskQuery, error) {
	normalized := query
	normalized.ProcessSubstring = strings.TrimSpace(query.ProcessSubstring)
	normalized.NameSubstring = strings.TrimSpace(query.NameSubstring)
	normalized.SampleSubstring = strings.TrimSpace(query.SampleSubstring)
	normalized.StatusRaw = strings.TrimSpace(query.StatusRaw)

	if normalized.StatusRaw != "" {
		status := trace.NormalizeTaskStatus(normalized.StatusRaw)
		if status == domain.TaskStatusUnknown {
			return domain.TaskQuery{}, fmt.Errorf("unknown task status %q", normalized.StatusRaw)
		}
		normalized.Status = status
		return normalized, nil
	}

	if normalized.Status != "" {
		rawStatus := strings.TrimSpace(string(normalized.Status))
		if rawStatus == "" {
			normalized.Status = ""
			return normalized, nil
		}

		status := trace.NormalizeTaskStatus(rawStatus)
		if status == domain.TaskStatusUnknown {
			if !strings.EqualFold(rawStatus, string(domain.TaskStatusUnknown)) {
				return domain.TaskQuery{}, fmt.Errorf("unknown task status %q", rawStatus)
			}
			normalized.Status = domain.TaskStatusUnknown
		} else {
			normalized.Status = status
		}
	}

	return normalized, nil
}

func ApplyTaskQuery(taskList []domain.Task, query domain.TaskQuery) ([]domain.Task, error) {
	normalized, err := NormalizeTaskQuery(query)
	if err != nil {
		return nil, err
	}

	containsFold := func(value string, needle string) bool {
		if needle == "" {
			return true
		}
		return strings.Contains(strings.ToLower(value), strings.ToLower(needle))
	}

	matches := make([]domain.Task, 0, len(taskList))
	for _, task := range taskList {
		if !containsFold(task.Process, normalized.ProcessSubstring) {
			continue
		}
		if !containsFold(task.Name, normalized.NameSubstring) {
			continue
		}
		if normalized.SampleSubstring != "" &&
			!containsFold(task.Name, normalized.SampleSubstring) &&
			!containsFold(task.Tag, normalized.SampleSubstring) {
			continue
		}
		if normalized.Status != "" && task.Status != normalized.Status {
			continue
		}

		matches = append(matches, task)
	}

	sort.SliceStable(matches, func(i, j int) bool {
		return matches[i].RowOrder < matches[j].RowOrder
	})

	return matches, nil
}

func BuildStatusSummary(runDir domain.RunDir, metadata domain.IndexMetadata, counts []domain.StatusCount, failedTasks []domain.Task) (domain.StatusSummary, error) {
	const failedPreviewLimit = 3

	copySource := func(source *domain.SourceFingerprint) *domain.SourceFingerprint {
		if source == nil {
			return nil
		}
		copied := *source
		return &copied
	}

	countsCopy := make([]domain.StatusCount, len(counts))
	copy(countsCopy, counts)

	builtAt := &metadata.BuiltAt
	if metadata.BuiltAt.IsZero() {
		builtAt = nil
	}

	if runDir.Path == "" && metadata.RunDir != "" {
		runDir = domain.RunDir{Path: metadata.RunDir}
	}

	return domain.StatusSummary{
		RunDir:    runDir,
		Mode:      metadata.Mode,
		IndexPath: metadata.IndexPath,
		Freshness: metadata.Freshness,
		BuiltAt:   builtAt,
		Sources: domain.ArtifactSet{
			RunDir: runDir,
			Mode:   metadata.Mode,
			Trace:  copySource(metadata.Trace),
			Log:    copySource(metadata.Log),
		},
		Counts:          countsCopy,
		FailedCount:     len(failedTasks),
		FailedPreview:   SelectFailedPreview(failedTasks, failedPreviewLimit),
		LogOnlyFailures: []domain.LogOnlyFailure{},
		Diagnostics:     []domain.Diagnostic{},
	}, nil
}

func SelectFailedPreview(failedTasks []domain.Task, limit int) []domain.FailedTaskPreview {
	if limit <= 0 || len(failedTasks) == 0 {
		return []domain.FailedTaskPreview{}
	}

	ordered := make([]domain.Task, len(failedTasks))
	copy(ordered, failedTasks)
	sort.SliceStable(ordered, func(i, j int) bool {
		return ordered[i].RowOrder < ordered[j].RowOrder
	})

	if limit > len(ordered) {
		limit = len(ordered)
	}

	previews := make([]domain.FailedTaskPreview, 0, limit)
	for _, task := range ordered[:limit] {
		var exit *int
		if task.Exit != nil {
			exitValue := *task.Exit
			exit = &exitValue
		}

		previews = append(previews, domain.FailedTaskPreview{
			ID:           task.ID,
			Status:       task.Status,
			Process:      task.Process,
			Name:         task.Name,
			Tag:          task.Tag,
			Workdir:      task.Workdir,
			Exit:         exit,
			ErrorSummary: task.ErrorSummary,
		})
	}

	return previews
}

type selectorMatchRule[T any] struct {
	match           func(string, T) bool
	exactDetail     string
	ambiguousDetail func(int) string
}

func resolveSelectorByRules[T any, R any](
	selector string,
	items []T,
	rules []selectorMatchRule[T],
	exact func(string, T, string) R,
	ambiguous func(string, []T, string) R,
	notFound func(string) R,
) R {
	normalizedSelector := strings.TrimSpace(selector)

	for _, rule := range rules {
		matches := make([]T, 0, 1)
		for _, item := range items {
			if rule.match(normalizedSelector, item) {
				matches = append(matches, item)
			}
		}

		if len(matches) == 1 {
			return exact(normalizedSelector, matches[0], rule.exactDetail)
		}
		if len(matches) > 1 {
			return ambiguous(normalizedSelector, matches, rule.ambiguousDetail(len(matches)))
		}
	}

	return notFound(normalizedSelector)
}

func taskSelectorRules() []selectorMatchRule[domain.Task] {
	return []selectorMatchRule[domain.Task]{
		{
			match:       MatchCanonicalID,
			exactDetail: "matched canonical id",
			ambiguousDetail: func(count int) string {
				return fmt.Sprintf("%d tasks matched canonical id; use a full workdir path if available", count)
			},
		},
		{
			match:       MatchWorkdirPath,
			exactDetail: "matched full workdir path",
			ambiguousDetail: func(count int) string {
				return fmt.Sprintf("%d tasks matched full workdir path; use a canonical id if available", count)
			},
		},
		{
			match:       MatchHumanSelector,
			exactDetail: "matched process/name/tag",
			ambiguousDetail: func(count int) string {
				return fmt.Sprintf("%d tasks matched process/name/tag; use a canonical id or full workdir path", count)
			},
		},
	}
}

func selectorDiagnostic(severity domain.DiagnosticSeverity, code string, message string, detail string) []domain.Diagnostic {
	return []domain.Diagnostic{
		{
			Severity: severity,
			Code:     code,
			Message:  message,
			Detail:   detail,
		},
	}
}

func selectorExactDiagnostic(detail string) []domain.Diagnostic {
	return selectorDiagnostic(domain.DiagnosticInfo, "selector_exact", "selector resolved exactly", detail)
}

func selectorAmbiguousDiagnostic(message string, detail string) []domain.Diagnostic {
	return selectorDiagnostic(domain.DiagnosticWarning, "selector_ambiguous", message, detail)
}

func selectorNotFoundDiagnostic(message string, detail string) []domain.Diagnostic {
	return selectorDiagnostic(domain.DiagnosticError, "selector_not_found", message, detail)
}

func taskExactSelectorResolution(selector string, task domain.Task, detail string) domain.SelectorResolution {
	matchedTask := task
	return domain.SelectorResolution{
		Kind:        domain.SelectorResolutionExact,
		Selector:    selector,
		Task:        &matchedTask,
		Matches:     []domain.Task{},
		Diagnostics: selectorExactDiagnostic(detail),
	}
}

func taskAmbiguousSelectorResolution(selector string, matches []domain.Task, detail string) domain.SelectorResolution {
	return domain.SelectorResolution{
		Kind:        domain.SelectorResolutionAmbiguous,
		Selector:    selector,
		Matches:     matches,
		Diagnostics: selectorAmbiguousDiagnostic("selector matched more than one task", detail),
	}
}

func taskNotFoundSelectorResolution(selector string) domain.SelectorResolution {
	return domain.SelectorResolution{
		Kind:     domain.SelectorResolutionNotFound,
		Selector: selector,
		Matches:  []domain.Task{},
		Diagnostics: selectorNotFoundDiagnostic(
			"selector did not match any indexed task",
			"searched canonical id, full workdir path, process, name, and tag",
		),
	}
}

func ResolveSelector(selector string, taskList []domain.Task) (domain.SelectorResolution, error) {
	return resolveSelectorByRules(
		selector,
		taskList,
		taskSelectorRules(),
		taskExactSelectorResolution,
		taskAmbiguousSelectorResolution,
		taskNotFoundSelectorResolution,
	), nil
}

func MatchCanonicalID(selector string, task domain.Task) bool {
	selectorID, ok := trace.NormalizeCanonicalTaskID(selector)
	if !ok {
		return false
	}

	taskID, ok := trace.NormalizeCanonicalTaskID(task.ID)
	if !ok {
		return false
	}

	return selectorID == taskID
}

func MatchWorkdirPath(selector string, task domain.Task) bool {
	normalizePath := func(value string) (string, bool) {
		value = strings.TrimSpace(value)
		if value == "" {
			return "", false
		}

		absolute, err := filepath.Abs(value)
		if err == nil {
			return absolute, true
		}

		return filepath.Clean(value), true
	}

	selectorPath, ok := normalizePath(selector)
	if !ok {
		return false
	}

	workdirPath, ok := normalizePath(task.Workdir)
	if !ok {
		return false
	}

	return selectorPath == workdirPath
}

func MatchHumanSelector(selector string, task domain.Task) bool {
	selector = strings.TrimSpace(selector)
	if selector == "" {
		return false
	}

	needle := strings.ToLower(selector)
	return strings.Contains(strings.ToLower(task.Process), needle) ||
		strings.Contains(strings.ToLower(task.Name), needle) ||
		strings.Contains(strings.ToLower(task.Tag), needle)
}

type selectorDossierState struct {
	kind        domain.SelectorResolutionKind
	selector    string
	diagnostics []domain.Diagnostic
	matchCount  int
}

func newSelectorDossierState(kind domain.SelectorResolutionKind, selector string, diagnostics []domain.Diagnostic, matchCount int) selectorDossierState {
	return selectorDossierState{kind: kind, selector: selector, diagnostics: diagnostics, matchCount: matchCount}
}

func buildSelectorDossier[D any](
	context string,
	resolution selectorDossierState,
	empty func([]domain.Diagnostic) D,
	exact func([]domain.Diagnostic) (D, error),
) (D, error) {
	diagnosticsCopy := append([]domain.Diagnostic(nil), resolution.diagnostics...)
	emptyDossier := empty(diagnosticsCopy)

	switch resolution.kind {
	case domain.SelectorResolutionExact:
		dossier, err := exact(diagnosticsCopy)
		if err != nil {
			return emptyDossier, err
		}
		return dossier, nil
	case domain.SelectorResolutionAmbiguous:
		return emptyDossier, fmt.Errorf("%s: selector %q is ambiguous (%d matches)", context, resolution.selector, resolution.matchCount)
	case domain.SelectorResolutionNotFound:
		return emptyDossier, fmt.Errorf("%s: selector %q not found", context, resolution.selector)
	default:
		return emptyDossier, fmt.Errorf("%s: unsupported selector resolution kind %q", context, resolution.kind)
	}
}

func BuildTaskDossier(resolution domain.SelectorResolution, inventory domain.CommandFileInventory) (domain.TaskDossier, error) {
	return buildSelectorDossier(
		"build task dossier",
		newSelectorDossierState(resolution.Kind, resolution.Selector, resolution.Diagnostics, len(resolution.Matches)),
		func(diagnostics []domain.Diagnostic) domain.TaskDossier {
			return domain.TaskDossier{Diagnostics: diagnostics}
		},
		func(diagnostics []domain.Diagnostic) (domain.TaskDossier, error) {
			if resolution.Task == nil {
				return domain.TaskDossier{}, fmt.Errorf("build task dossier: exact selector %q has no resolved task", resolution.Selector)
			}
			return domain.TaskDossier{
				Task:        *resolution.Task,
				Inventory:   inventory,
				Diagnostics: diagnostics,
			}, nil
		},
	)
}

func logOnlyEvidenceSelectorRules() []selectorMatchRule[domain.LogOnlyTaskEvidence] {
	return []selectorMatchRule[domain.LogOnlyTaskEvidence]{
		{
			match: func(selector string, evidence domain.LogOnlyTaskEvidence) bool {
				return MatchCanonicalID(selector, domain.Task{ID: evidence.ID})
			},
			exactDetail: "matched canonical id",
			ambiguousDetail: func(count int) string {
				return fmt.Sprintf("%d log-only evidence rows matched canonical id; use a full workdir path if available", count)
			},
		},
		{
			match: func(selector string, evidence domain.LogOnlyTaskEvidence) bool {
				return MatchWorkdirPath(selector, domain.Task{Workdir: evidence.Workdir})
			},
			exactDetail: "matched full workdir path",
			ambiguousDetail: func(count int) string {
				return fmt.Sprintf("%d log-only evidence rows matched full workdir path; use a canonical id if available", count)
			},
		},
		{
			match:       MatchLogOnlyEvidenceSelector,
			exactDetail: "matched log-only evidence",
			ambiguousDetail: func(count int) string {
				return fmt.Sprintf("%d log-only evidence rows matched observed log-only fields; use a canonical id or full workdir path if available", count)
			},
		},
	}
}

func logOnlyExactSelectorResolution(selector string, evidence domain.LogOnlyTaskEvidence, detail string) domain.LogOnlySelectorResolution {
	matchedEvidence := evidence
	return domain.LogOnlySelectorResolution{
		Kind:        domain.SelectorResolutionExact,
		Selector:    selector,
		Evidence:    &matchedEvidence,
		Matches:     []domain.LogOnlyTaskEvidence{},
		Diagnostics: selectorExactDiagnostic(detail),
	}
}

func logOnlyAmbiguousSelectorResolution(selector string, matches []domain.LogOnlyTaskEvidence, detail string) domain.LogOnlySelectorResolution {
	return domain.LogOnlySelectorResolution{
		Kind:        domain.SelectorResolutionAmbiguous,
		Selector:    selector,
		Matches:     matches,
		Diagnostics: selectorAmbiguousDiagnostic("selector matched more than one log-only evidence row", detail),
	}
}

func logOnlyNotFoundSelectorResolution(selector string) domain.LogOnlySelectorResolution {
	return domain.LogOnlySelectorResolution{
		Kind:     domain.SelectorResolutionNotFound,
		Selector: selector,
		Matches:  []domain.LogOnlyTaskEvidence{},
		Diagnostics: selectorNotFoundDiagnostic(
			"selector did not match any log-only evidence",
			"searched observed canonical id, full workdir path, process, name, and display label",
		),
	}
}

func ResolveLogOnlySelector(selector string, evidence []domain.LogOnlyTaskEvidence) (domain.LogOnlySelectorResolution, error) {
	return resolveSelectorByRules(
		selector,
		evidence,
		logOnlyEvidenceSelectorRules(),
		logOnlyExactSelectorResolution,
		logOnlyAmbiguousSelectorResolution,
		logOnlyNotFoundSelectorResolution,
	), nil
}

func MatchLogOnlyEvidenceSelector(selector string, evidence domain.LogOnlyTaskEvidence) bool {
	selector = strings.TrimSpace(selector)
	if selector == "" {
		return false
	}

	if MatchCanonicalID(selector, domain.Task{ID: evidence.ID}) {
		return true
	}
	if MatchWorkdirPath(selector, domain.Task{Workdir: evidence.Workdir}) {
		return true
	}

	needle := strings.ToLower(selector)
	containsSelector := func(value string) bool {
		value = strings.TrimSpace(value)
		return value != "" && strings.Contains(strings.ToLower(value), needle)
	}
	if containsSelector(evidence.Process) || containsSelector(evidence.Name) {
		return true
	}

	process := strings.Join(strings.Fields(evidence.Process), " ")
	name := strings.Join(strings.Fields(evidence.Name), " ")
	if process == "" || name == "" {
		return false
	}

	displayStrings := []string{
		process + " (" + name + ")",
		process + "/" + name,
		process + " " + name,
	}
	for _, display := range displayStrings {
		if strings.Contains(strings.ToLower(display), needle) {
			return true
		}
	}

	return false
}

func BuildLogOnlyTaskDossier(resolution domain.LogOnlySelectorResolution, inventory domain.CommandFileInventory) (domain.LogOnlyTaskDossier, error) {
	return buildSelectorDossier(
		"build log-only task dossier",
		newSelectorDossierState(resolution.Kind, resolution.Selector, resolution.Diagnostics, len(resolution.Matches)),
		func(diagnostics []domain.Diagnostic) domain.LogOnlyTaskDossier {
			return domain.LogOnlyTaskDossier{Diagnostics: diagnostics}
		},
		func(diagnostics []domain.Diagnostic) (domain.LogOnlyTaskDossier, error) {
			if resolution.Evidence == nil {
				return domain.LogOnlyTaskDossier{}, fmt.Errorf("build log-only task dossier: exact selector %q has no resolved evidence", resolution.Selector)
			}
			return domain.LogOnlyTaskDossier{
				Evidence:    *resolution.Evidence,
				Inventory:   inventory,
				Diagnostics: diagnostics,
			}, nil
		},
	)
}
