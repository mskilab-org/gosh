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

func ResolveSelector(selector string, taskList []domain.Task) (domain.SelectorResolution, error) {
	normalizedSelector := strings.TrimSpace(selector)

	exactResolution := func(task domain.Task, detail string) domain.SelectorResolution {
		matchedTask := task
		return domain.SelectorResolution{
			Kind:     domain.SelectorResolutionExact,
			Selector: normalizedSelector,
			Task:     &matchedTask,
			Matches:  []domain.Task{},
			Diagnostics: []domain.Diagnostic{
				{
					Severity: domain.DiagnosticInfo,
					Code:     "selector_exact",
					Message:  "selector resolved exactly",
					Detail:   detail,
				},
			},
		}
	}

	ambiguousResolution := func(matches []domain.Task, detail string) domain.SelectorResolution {
		return domain.SelectorResolution{
			Kind:     domain.SelectorResolutionAmbiguous,
			Selector: normalizedSelector,
			Matches:  matches,
			Diagnostics: []domain.Diagnostic{
				{
					Severity: domain.DiagnosticWarning,
					Code:     "selector_ambiguous",
					Message:  "selector matched more than one task",
					Detail:   detail,
				},
			},
		}
	}

	resolveExactMatches := func(matches []domain.Task, exactDetail string, ambiguousDetail string) (domain.SelectorResolution, bool) {
		if len(matches) == 1 {
			return exactResolution(matches[0], exactDetail), true
		}
		if len(matches) > 1 {
			return ambiguousResolution(matches, ambiguousDetail), true
		}
		return domain.SelectorResolution{}, false
	}

	canonicalMatches := make([]domain.Task, 0, 1)
	for _, task := range taskList {
		if MatchCanonicalID(normalizedSelector, task) {
			canonicalMatches = append(canonicalMatches, task)
		}
	}
	if resolution, ok := resolveExactMatches(
		canonicalMatches,
		"matched canonical id",
		fmt.Sprintf("%d tasks matched canonical id; use a full workdir path if available", len(canonicalMatches)),
	); ok {
		return resolution, nil
	}

	workdirMatches := make([]domain.Task, 0, 1)
	for _, task := range taskList {
		if MatchWorkdirPath(normalizedSelector, task) {
			workdirMatches = append(workdirMatches, task)
		}
	}
	if resolution, ok := resolveExactMatches(
		workdirMatches,
		"matched full workdir path",
		fmt.Sprintf("%d tasks matched full workdir path; use a canonical id if available", len(workdirMatches)),
	); ok {
		return resolution, nil
	}

	humanMatches := make([]domain.Task, 0)
	for _, task := range taskList {
		if MatchHumanSelector(normalizedSelector, task) {
			humanMatches = append(humanMatches, task)
		}
	}
	if resolution, ok := resolveExactMatches(
		humanMatches,
		"matched process/name/tag",
		fmt.Sprintf("%d tasks matched process/name/tag; use a canonical id or full workdir path", len(humanMatches)),
	); ok {
		return resolution, nil
	}

	return domain.SelectorResolution{
		Kind:     domain.SelectorResolutionNotFound,
		Selector: normalizedSelector,
		Matches:  []domain.Task{},
		Diagnostics: []domain.Diagnostic{
			{
				Severity: domain.DiagnosticError,
				Code:     "selector_not_found",
				Message:  "selector did not match any indexed task",
				Detail:   "searched canonical id, full workdir path, process, name, and tag",
			},
		},
	}, nil
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

func BuildTaskDossier(resolution domain.SelectorResolution, inventory domain.CommandFileInventory) (domain.TaskDossier, error) {
	diagnostics := append([]domain.Diagnostic(nil), resolution.Diagnostics...)
	emptyDossier := domain.TaskDossier{Diagnostics: diagnostics}

	switch resolution.Kind {
	case domain.SelectorResolutionExact:
		if resolution.Task == nil {
			return emptyDossier, fmt.Errorf("build task dossier: exact selector %q has no resolved task", resolution.Selector)
		}
		return domain.TaskDossier{
			Task:        *resolution.Task,
			Inventory:   inventory,
			Diagnostics: diagnostics,
		}, nil
	case domain.SelectorResolutionAmbiguous:
		return emptyDossier, fmt.Errorf("build task dossier: selector %q is ambiguous (%d matches)", resolution.Selector, len(resolution.Matches))
	case domain.SelectorResolutionNotFound:
		return emptyDossier, fmt.Errorf("build task dossier: selector %q not found", resolution.Selector)
	default:
		return emptyDossier, fmt.Errorf("build task dossier: unsupported selector resolution kind %q", resolution.Kind)
	}
}
