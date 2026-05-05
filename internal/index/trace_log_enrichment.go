package index

import (
	"context"
	"fmt"
	"path/filepath"
	"strings"

	"github.com/mskilab-org/gosh/internal/domain"
	taskmatch "github.com/mskilab-org/gosh/internal/tasks"
	"github.com/mskilab-org/gosh/internal/trace"
)

// EnrichTraceTasksWithLogEvidence returns the trace-backed task rows with any
// deterministic paired-log enrichments applied. It must preserve the trace row
// set and source order; log evidence may only fill missing Workdir and
// ErrorSummary fields for failed/aborted rows when there is exactly one safe
// match.
func EnrichTraceTasksWithLogEvidence(ctx context.Context, runDir domain.RunDir, tasks []domain.Task, evidence []domain.LogOnlyTaskEvidence) (domain.TraceLogEnrichmentResult, error) {
	if ctx == nil {
		return domain.TraceLogEnrichmentResult{}, fmt.Errorf("enrich trace tasks with log evidence: nil context")
	}
	if err := ctx.Err(); err != nil {
		return domain.TraceLogEnrichmentResult{}, fmt.Errorf("enrich trace tasks with log evidence: %w", err)
	}

	enrichedTasks := append([]domain.Task(nil), tasks...)
	diagnostics := make([]domain.Diagnostic, 0)
	for index, task := range tasks {
		if err := ctx.Err(); err != nil {
			return domain.TraceLogEnrichmentResult{}, fmt.Errorf("enrich trace tasks with log evidence: %w", err)
		}

		if !IsTraceTaskEligibleForLogEnrichment(task) {
			continue
		}

		matchedEvidence, diagnostic := FindDeterministicTraceLogEvidenceMatch(task, evidence)
		if diagnostic != nil {
			diagnostics = append(diagnostics, *diagnostic)
		}
		if matchedEvidence == nil {
			continue
		}

		enrichedTasks[index] = MergeTraceTaskLogEvidence(runDir, task, *matchedEvidence)
	}

	return domain.TraceLogEnrichmentResult{Tasks: enrichedTasks, Diagnostics: diagnostics}, nil
}

// IsTraceTaskEligibleForLogEnrichment reports whether a trace-backed task row is
// allowed to receive paired-log enrichment.
func IsTraceTaskEligibleForLogEnrichment(task domain.Task) bool {
	if task.Status != domain.TaskStatusFailed && task.Status != domain.TaskStatusAborted {
		return false
	}

	return strings.TrimSpace(task.Workdir) == "" || strings.TrimSpace(task.ErrorSummary) == ""
}

// FindDeterministicTraceLogEvidenceMatch finds the single log evidence row that
// can safely enrich a trace-backed task, or returns a diagnostic/no-match result
// when evidence is absent or ambiguous.
func FindDeterministicTraceLogEvidenceMatch(task domain.Task, evidence []domain.LogOnlyTaskEvidence) (*domain.LogOnlyTaskEvidence, *domain.Diagnostic) {
	canonicalIDFromValue := func(value string) (string, bool) {
		if canonicalID, ok := trace.NormalizeCanonicalTaskID(value); ok {
			return canonicalID, true
		}
		canonicalID, err := trace.DeriveCanonicalTaskID(value)
		if err != nil {
			return "", false
		}
		return canonicalID, true
	}

	taskCanonicalID, hasTaskCanonicalID := canonicalIDFromValue(task.ID)
	if !hasTaskCanonicalID {
		taskCanonicalID, hasTaskCanonicalID = canonicalIDFromValue(task.Workdir)
	}

	hasTaskWorkdir := strings.TrimSpace(task.Workdir) != ""

	normalizeText := func(value string) string {
		return strings.ToLower(strings.Join(strings.Fields(value), " "))
	}

	taskNameParts := trace.DeriveNameParts(task.Name)
	taskNameKeys := make(map[string]bool)
	for _, value := range []string{task.Name, task.Tag, taskNameParts.FullName, taskNameParts.Tag} {
		if key := normalizeText(value); key != "" {
			taskNameKeys[key] = true
		}
	}
	taskProcessKeys := make(map[string]bool)
	for _, value := range []string{task.Process, taskNameParts.Process} {
		if key := normalizeText(value); key != "" {
			taskProcessKeys[key] = true
		}
	}

	evidenceHasCanonicalID := func(item domain.LogOnlyTaskEvidence) bool {
		if !hasTaskCanonicalID {
			return false
		}
		for _, value := range []string{item.ID, item.Workdir} {
			canonicalID, ok := canonicalIDFromValue(value)
			if ok && canonicalID == taskCanonicalID {
				return true
			}
		}
		return false
	}

	evidenceHasWorkdir := func(item domain.LogOnlyTaskEvidence) bool {
		return hasTaskWorkdir && taskmatch.MatchWorkdirPath(task.Workdir, domain.Task{Workdir: item.Workdir})
	}

	evidenceHasProcessName := func(item domain.LogOnlyTaskEvidence) bool {
		evidenceProcess := normalizeText(item.Process)
		evidenceName := normalizeText(item.Name)
		return evidenceProcess != "" && evidenceName != "" && taskProcessKeys[evidenceProcess] && taskNameKeys[evidenceName]
	}

	evidenceHasReferencedWorkdir := func(item domain.LogOnlyTaskEvidence) bool {
		for _, source := range item.Sources {
			if source.Kind != domain.LogOnlyEvidenceSourceWorkdir {
				continue
			}
			if hasTaskWorkdir && taskmatch.MatchWorkdirPath(task.Workdir, domain.Task{Workdir: source.Path}) {
				return true
			}
			if hasTaskCanonicalID {
				canonicalID, ok := canonicalIDFromValue(source.Path)
				if ok && canonicalID == taskCanonicalID {
					return true
				}
			}
		}
		return false
	}

	taskIdentity := func() string {
		if hasTaskCanonicalID {
			return taskCanonicalID
		}
		if strings.TrimSpace(task.Workdir) != "" {
			return strings.TrimSpace(task.Workdir)
		}
		if strings.TrimSpace(task.Process) != "" || strings.TrimSpace(task.Name) != "" {
			return strings.TrimSpace(strings.TrimSpace(task.Process) + " " + strings.TrimSpace(task.Name))
		}
		return "<unknown>"
	}

	type matchRule struct {
		name  string
		match func(domain.LogOnlyTaskEvidence) bool
	}
	rules := []matchRule{
		{name: "canonical hash", match: evidenceHasCanonicalID},
		{name: "workdir path", match: evidenceHasWorkdir},
		{name: "process/name", match: evidenceHasProcessName},
		{name: "referenced workdir", match: evidenceHasReferencedWorkdir},
	}

	for _, rule := range rules {
		matches := make([]domain.LogOnlyTaskEvidence, 0, 1)
		for _, item := range evidence {
			if rule.match(item) {
				matches = append(matches, item)
			}
		}

		if len(matches) == 1 {
			matched := matches[0]
			return &matched, nil
		}
		if len(matches) > 1 {
			return nil, &domain.Diagnostic{
				Severity: domain.DiagnosticWarning,
				Code:     "trace_log_evidence_ambiguous",
				Message:  "multiple log evidence rows matched trace task",
				Detail:   fmt.Sprintf("%d log evidence rows matched trace task %s by %s; leaving trace row unchanged", len(matches), taskIdentity(), rule.name),
			}
		}
	}

	return nil, &domain.Diagnostic{
		Severity: domain.DiagnosticInfo,
		Code:     "trace_log_evidence_no_match",
		Message:  "no deterministic log evidence matched trace task",
		Detail:   fmt.Sprintf("searched canonical hash, workdir path, process/name, and referenced workdir for trace task %s", taskIdentity()),
	}
}

// MergeTraceTaskLogEvidence fills only missing enrichable fields on a
// trace-backed task row from deterministic paired-log evidence.
func MergeTraceTaskLogEvidence(runDir domain.RunDir, task domain.Task, evidence domain.LogOnlyTaskEvidence) domain.Task {
	merged := task

	if strings.TrimSpace(merged.Workdir) == "" && strings.TrimSpace(evidence.Workdir) != "" {
		canonicalID := strings.TrimSpace(task.ID)
		if canonicalID == "" {
			canonicalID = strings.TrimSpace(evidence.ID)
		}

		workdir, err := trace.ResolveTaskWorkdir(runDir, canonicalID, evidence.Workdir)
		if err == nil && strings.TrimSpace(workdir) != "" {
			merged.Workdir = filepath.Clean(workdir)
		}
	}

	if strings.TrimSpace(merged.ErrorSummary) == "" && strings.TrimSpace(evidence.ErrorSummary) != "" {
		merged.ErrorSummary = evidence.ErrorSummary
	}

	return merged
}
