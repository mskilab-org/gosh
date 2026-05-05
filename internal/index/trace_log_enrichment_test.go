package index

import (
	"context"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/mskilab-org/gosh/internal/domain"
)

func TestEnrichTraceTasksWithLogEvidencePreservesTraceRowsAndEnrichesEligibleMatches(t *testing.T) {
	runPath := t.TempDir()
	failedWorkdir := filepath.Join(runPath, "work", "aa", "111111")
	abortedWorkdir := filepath.Join(runPath, "work", "bb", "222222")
	traceExit := 1

	tasks := []domain.Task{
		{RowOrder: 10, ID: "cc/333333", Status: domain.TaskStatusCompleted, Process: "QC", Name: "QC (normal)"},
		{RowOrder: 20, ID: "aa/111111", Status: domain.TaskStatusFailed, Process: "ALIGN", Name: "ALIGN (tumor)", Tag: "tumor", Exit: &traceExit},
		{RowOrder: 30, ID: "bb/222222", Status: domain.TaskStatusAborted, Process: "CALL", Name: "CALL (tumor)", Workdir: abortedWorkdir},
	}
	evidence := []domain.LogOnlyTaskEvidence{
		{ID: "cc/333333", Workdir: filepath.Join(runPath, "work", "cc", "333333"), ErrorSummary: "completed row should not be enriched"},
		{ID: "AA/111111", Workdir: failedWorkdir, ErrorSummary: "failed after dispatch"},
		{ID: "bb/222222", Workdir: filepath.Join(runPath, "work", "bb", "ignored"), ErrorSummary: "aborted by scheduler"},
	}

	got, err := EnrichTraceTasksWithLogEvidence(context.Background(), domain.RunDir{Path: runPath}, tasks, evidence)
	if err != nil {
		t.Fatalf("EnrichTraceTasksWithLogEvidence() error = %v, want nil", err)
	}
	if len(got.Diagnostics) != 0 {
		t.Fatalf("diagnostics = %#v, want none", got.Diagnostics)
	}
	if len(got.Tasks) != len(tasks) {
		t.Fatalf("tasks length = %d, want %d", len(got.Tasks), len(tasks))
	}
	if got.Tasks[0] != tasks[0] {
		t.Fatalf("completed task = %#v, want unchanged %#v", got.Tasks[0], tasks[0])
	}
	if got.Tasks[1].RowOrder != 20 || got.Tasks[2].RowOrder != 30 {
		t.Fatalf("row order = [%d, %d], want [20, 30]", got.Tasks[1].RowOrder, got.Tasks[2].RowOrder)
	}
	if got.Tasks[1].Workdir != filepath.Clean(failedWorkdir) || got.Tasks[1].ErrorSummary != "failed after dispatch" {
		t.Fatalf("failed task enrichment = %#v, want workdir and error summary from deterministic evidence", got.Tasks[1])
	}
	if got.Tasks[1].Exit != tasks[1].Exit || got.Tasks[1].Process != tasks[1].Process || got.Tasks[1].Name != tasks[1].Name || got.Tasks[1].Tag != tasks[1].Tag {
		t.Fatalf("failed task trace-owned fields = %#v, want preserved from %#v", got.Tasks[1], tasks[1])
	}
	if got.Tasks[2].Workdir != abortedWorkdir || got.Tasks[2].ErrorSummary != "aborted by scheduler" {
		t.Fatalf("aborted task enrichment = %#v, want existing workdir plus evidence summary", got.Tasks[2])
	}
}

func TestEnrichTraceTasksWithLogEvidenceReportsDiagnosticsWithoutHardError(t *testing.T) {
	tasks := []domain.Task{
		{RowOrder: 1, ID: "aa/111111", Status: domain.TaskStatusFailed},
		{RowOrder: 2, ID: "bb/222222", Status: domain.TaskStatusAborted},
	}
	evidence := []domain.LogOnlyTaskEvidence{
		{ID: "cc/333333", ErrorSummary: "wrong task"},
		{ID: "bb/222222", ErrorSummary: "first candidate"},
		{ID: "BB/222222", ErrorSummary: "second candidate"},
	}

	got, err := EnrichTraceTasksWithLogEvidence(context.Background(), domain.RunDir{Path: t.TempDir()}, tasks, evidence)
	if err != nil {
		t.Fatalf("EnrichTraceTasksWithLogEvidence() error = %v, want nil", err)
	}
	if len(got.Tasks) != len(tasks) || got.Tasks[0] != tasks[0] || got.Tasks[1] != tasks[1] {
		t.Fatalf("tasks = %#v, want unchanged %#v", got.Tasks, tasks)
	}
	if len(got.Diagnostics) != 2 {
		t.Fatalf("diagnostics = %#v, want no-match and ambiguous diagnostics", got.Diagnostics)
	}
	if got.Diagnostics[0].Code != "trace_log_evidence_no_match" || got.Diagnostics[0].Severity != domain.DiagnosticInfo {
		t.Fatalf("first diagnostic = %#v, want no-match info", got.Diagnostics[0])
	}
	if got.Diagnostics[1].Code != "trace_log_evidence_ambiguous" || got.Diagnostics[1].Severity != domain.DiagnosticWarning {
		t.Fatalf("second diagnostic = %#v, want ambiguous warning", got.Diagnostics[1])
	}
}

func TestEnrichTraceTasksWithLogEvidenceRejectsNilOrCanceledContext(t *testing.T) {
	if got, err := EnrichTraceTasksWithLogEvidence(nil, domain.RunDir{}, nil, nil); err == nil {
		t.Fatalf("EnrichTraceTasksWithLogEvidence(nil context) returned nil error and result %#v", got)
	} else if !strings.Contains(err.Error(), "nil context") {
		t.Fatalf("error = %q, want nil context message", err.Error())
	}

	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if got, err := EnrichTraceTasksWithLogEvidence(ctx, domain.RunDir{}, []domain.Task{{Status: domain.TaskStatusFailed}}, nil); err == nil {
		t.Fatalf("EnrichTraceTasksWithLogEvidence(canceled context) returned nil error and result %#v", got)
	} else if !strings.Contains(err.Error(), context.Canceled.Error()) {
		t.Fatalf("error = %q, want context canceled message", err.Error())
	}
}

func TestIsTraceTaskEligibleForLogEnrichmentAllowsFailedOrAbortedRowsWithMissingEnrichableFields(t *testing.T) {
	tests := []struct {
		name string
		task domain.Task
	}{
		{
			name: "failed row missing workdir and error summary",
			task: domain.Task{Status: domain.TaskStatusFailed},
		},
		{
			name: "failed row missing workdir only",
			task: domain.Task{Status: domain.TaskStatusFailed, ErrorSummary: "command exited with 137"},
		},
		{
			name: "aborted row missing error summary only",
			task: domain.Task{Status: domain.TaskStatusAborted, Workdir: "/runs/example/work/aa/111111"},
		},
		{
			name: "aborted row missing workdir only",
			task: domain.Task{Status: domain.TaskStatusAborted, ErrorSummary: "workflow aborted"},
		},
		{
			name: "failed row with whitespace-only enrichable fields",
			task: domain.Task{Status: domain.TaskStatusFailed, Workdir: " \t", ErrorSummary: "\n"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if !IsTraceTaskEligibleForLogEnrichment(tt.task) {
				t.Fatalf("IsTraceTaskEligibleForLogEnrichment(%#v) = false, want true", tt.task)
			}
		})
	}
}

func TestFindDeterministicTraceLogEvidenceMatchMatchesBySafeIdentityRules(t *testing.T) {
	root := t.TempDir()
	workdir := filepath.Join(root, "work", "aa", "111111")
	referencedWorkdir := filepath.Join(root, "work", "dd", "444444")

	tests := []struct {
		name     string
		task     domain.Task
		evidence []domain.LogOnlyTaskEvidence
		want     domain.LogOnlyTaskEvidence
	}{
		{
			name: "canonical hash ignores case and surrounding space",
			task: domain.Task{ID: "aa/111111"},
			evidence: []domain.LogOnlyTaskEvidence{
				{ID: "bb/222222", ErrorSummary: "wrong evidence"},
				{ID: " AA/111111\t", ErrorSummary: "matched by id"},
			},
			want: domain.LogOnlyTaskEvidence{ID: " AA/111111\t", ErrorSummary: "matched by id"},
		},
		{
			name: "full workdir path when canonical hash is unavailable",
			task: domain.Task{Workdir: filepath.Join(filepath.Dir(workdir), ".", filepath.Base(workdir))},
			evidence: []domain.LogOnlyTaskEvidence{
				{Workdir: filepath.Join(root, "work", "cc", "333333"), ErrorSummary: "wrong workdir"},
				{Workdir: workdir, ErrorSummary: "matched by workdir"},
			},
			want: domain.LogOnlyTaskEvidence{Workdir: workdir, ErrorSummary: "matched by workdir"},
		},
		{
			name: "process and trace tag tuple",
			task: domain.Task{Process: "ALIGN_STAR", Name: "ALIGN_STAR (tumor-1)", Tag: "tumor-1"},
			evidence: []domain.LogOnlyTaskEvidence{
				{Process: "ALIGN_STAR", Name: "normal-1", ErrorSummary: "wrong sample"},
				{Process: " align_star ", Name: " tumor-1 ", ErrorSummary: "matched by process/name"},
			},
			want: domain.LogOnlyTaskEvidence{Process: " align_star ", Name: " tumor-1 ", ErrorSummary: "matched by process/name"},
		},
		{
			name: "unique referenced workdir source",
			task: domain.Task{ID: "dd/444444"},
			evidence: []domain.LogOnlyTaskEvidence{
				{ID: "", ErrorSummary: "wrong source", Sources: []domain.LogOnlyEvidenceSource{{Kind: domain.LogOnlyEvidenceSourceWorkdir, Path: filepath.Join(root, "work", "ee", "555555")}}},
				{ID: "", ErrorSummary: "matched by referenced workdir", Sources: []domain.LogOnlyEvidenceSource{{Kind: domain.LogOnlyEvidenceSourceWorkdir, Path: referencedWorkdir}}},
			},
			want: domain.LogOnlyTaskEvidence{ID: "", ErrorSummary: "matched by referenced workdir", Sources: []domain.LogOnlyEvidenceSource{{Kind: domain.LogOnlyEvidenceSourceWorkdir, Path: referencedWorkdir}}},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, diagnostic := FindDeterministicTraceLogEvidenceMatch(tt.task, tt.evidence)
			if diagnostic != nil {
				t.Fatalf("FindDeterministicTraceLogEvidenceMatch() diagnostic = %#v, want nil", *diagnostic)
			}
			if got == nil {
				t.Fatalf("FindDeterministicTraceLogEvidenceMatch() evidence = nil, want %#v", tt.want)
			}
			assertLogEvidenceEqual(t, *got, tt.want)
		})
	}
}

func TestFindDeterministicTraceLogEvidenceMatchReportsNoMatch(t *testing.T) {
	tests := []struct {
		name     string
		evidence []domain.LogOnlyTaskEvidence
	}{
		{
			name:     "empty evidence",
			evidence: []domain.LogOnlyTaskEvidence{},
		},
		{
			name:     "unrelated evidence",
			evidence: []domain.LogOnlyTaskEvidence{{ID: "bb/222222", Process: "CALL", Name: "normal"}},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, diagnostic := FindDeterministicTraceLogEvidenceMatch(
				domain.Task{ID: "aa/111111", Workdir: filepath.Join(t.TempDir(), "work", "aa", "111111"), Process: "ALIGN", Name: "ALIGN (tumor)", Tag: "tumor"},
				tt.evidence,
			)

			if got != nil {
				t.Fatalf("FindDeterministicTraceLogEvidenceMatch() evidence = %#v, want nil on no match", *got)
			}
			assertTraceLogEvidenceDiagnostic(t, diagnostic, domain.DiagnosticInfo, "trace_log_evidence_no_match", "no deterministic log evidence matched trace task")
		})
	}
}

func TestFindDeterministicTraceLogEvidenceMatchReportsAmbiguityWithoutGuessing(t *testing.T) {
	got, diagnostic := FindDeterministicTraceLogEvidenceMatch(
		domain.Task{ID: "aa/111111", Process: "ALIGN", Name: "ALIGN (tumor)", Tag: "tumor"},
		[]domain.LogOnlyTaskEvidence{
			{ID: "aa/111111", Process: "ALIGN", Name: "tumor", ErrorSummary: "first candidate"},
			{ID: "AA/111111", Process: "ALIGN", Name: "tumor", ErrorSummary: "second candidate"},
		},
	)

	if got != nil {
		t.Fatalf("FindDeterministicTraceLogEvidenceMatch() evidence = %#v, want nil on ambiguity", *got)
	}
	assertTraceLogEvidenceDiagnostic(t, diagnostic, domain.DiagnosticWarning, "trace_log_evidence_ambiguous", "multiple log evidence rows matched trace task")
}

func TestIsTraceTaskEligibleForLogEnrichmentRejectsRowsThatCannotBeImproved(t *testing.T) {
	completeFailedTask := domain.Task{
		Status:       domain.TaskStatusFailed,
		Workdir:      "/runs/example/work/bb/222222",
		ErrorSummary: "command exited with 137",
	}

	tests := []struct {
		name string
		task domain.Task
	}{
		{
			name: "failed row already has workdir and error summary",
			task: completeFailedTask,
		},
		{
			name: "aborted row already has workdir and error summary",
			task: domain.Task{Status: domain.TaskStatusAborted, Workdir: "/runs/example/work/cc/333333", ErrorSummary: "workflow aborted"},
		},
		{
			name: "completed row with missing fields",
			task: domain.Task{Status: domain.TaskStatusCompleted},
		},
		{
			name: "cached row with missing fields",
			task: domain.Task{Status: domain.TaskStatusCached},
		},
		{
			name: "unknown row with missing fields",
			task: domain.Task{Status: domain.TaskStatusUnknown},
		},
		{
			name: "submitted row with missing fields",
			task: domain.Task{Status: domain.TaskStatusSubmitted},
		},
		{
			name: "running row with missing fields",
			task: domain.Task{Status: domain.TaskStatusRunning},
		},
		{
			name: "empty status row with missing fields",
			task: domain.Task{},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if IsTraceTaskEligibleForLogEnrichment(tt.task) {
				t.Fatalf("IsTraceTaskEligibleForLogEnrichment(%#v) = true, want false", tt.task)
			}
		})
	}
}

func TestMergeTraceTaskLogEvidenceFillsOnlyMissingEnrichableFields(t *testing.T) {
	workdir := filepath.Join(t.TempDir(), "work", "aa", "111111")
	traceExit := 1
	evidenceExit := 137
	task := domain.Task{
		RowOrder:     42,
		ID:           "aa/111111",
		Status:       domain.TaskStatusFailed,
		Process:      "TRACE_PROCESS",
		Name:         "TRACE_PROCESS (tumor-1)",
		Tag:          "tumor-1",
		Workdir:      " \t",
		Exit:         &traceExit,
		Duration:     "2m",
		Realtime:     "120s",
		CPUs:         "4",
		Memory:       "8 GB",
		ErrorSummary: "\n",
	}
	evidence := domain.LogOnlyTaskEvidence{
		ID:             "bb/222222",
		Workdir:        workdir,
		Process:        "LOG_PROCESS",
		Name:           "different sample",
		ObservedStatus: domain.TaskStatusAborted,
		Exit:           &evidenceExit,
		ErrorSummary:   "killed by scheduler",
		ErrorBlock:     "large log block",
	}

	got := MergeTraceTaskLogEvidence(domain.RunDir{Path: filepath.Dir(filepath.Dir(filepath.Dir(workdir)))}, task, evidence)

	if got.Workdir != filepath.Clean(workdir) {
		t.Fatalf("merged Workdir = %q, want %q", got.Workdir, filepath.Clean(workdir))
	}
	if got.ErrorSummary != evidence.ErrorSummary {
		t.Fatalf("merged ErrorSummary = %q, want %q", got.ErrorSummary, evidence.ErrorSummary)
	}
	assertTraceOwnedTaskFieldsPreserved(t, got, task)
}

func TestMergeTraceTaskLogEvidencePreservesExistingEnrichableFields(t *testing.T) {
	task := domain.Task{
		ID:           "aa/111111",
		Status:       domain.TaskStatusFailed,
		Workdir:      filepath.Join(t.TempDir(), "work", "aa", "111111"),
		ErrorSummary: "trace-owned error summary",
	}
	evidence := domain.LogOnlyTaskEvidence{
		Workdir:      filepath.Join(t.TempDir(), "work", "bb", "222222"),
		ErrorSummary: "log evidence error summary",
	}

	got := MergeTraceTaskLogEvidence(domain.RunDir{Path: t.TempDir()}, task, evidence)

	if got.Workdir != task.Workdir {
		t.Fatalf("merged Workdir = %q, want existing trace workdir %q", got.Workdir, task.Workdir)
	}
	if got.ErrorSummary != task.ErrorSummary {
		t.Fatalf("merged ErrorSummary = %q, want existing trace error summary %q", got.ErrorSummary, task.ErrorSummary)
	}
}

func TestMergeTraceTaskLogEvidenceResolvesHashOnlyEvidenceWorkdirFromRunDir(t *testing.T) {
	runPath := t.TempDir()
	resolvedWorkdir := filepath.Join(runPath, "work", "9e", "c300c50150c213c2d44ca9e4624d8c")
	if err := os.MkdirAll(resolvedWorkdir, 0o755); err != nil {
		t.Fatalf("create resolved workdir: %v", err)
	}

	got := MergeTraceTaskLogEvidence(
		domain.RunDir{Path: runPath},
		domain.Task{ID: "9e/c300c5", Status: domain.TaskStatusFailed},
		domain.LogOnlyTaskEvidence{ID: "9e/c300c5", Workdir: "9e/c300c5", ErrorSummary: "failed after dispatch"},
	)

	if got.Workdir != filepath.Clean(resolvedWorkdir) {
		t.Fatalf("merged Workdir = %q, want resolved run workdir %q", got.Workdir, filepath.Clean(resolvedWorkdir))
	}
	if got.ID != "9e/c300c5" {
		t.Fatalf("merged ID = %q, want trace canonical ID preserved", got.ID)
	}
}

func TestMergeTraceTaskLogEvidenceIgnoresBlankEvidenceFields(t *testing.T) {
	task := domain.Task{ID: "aa/111111", Status: domain.TaskStatusFailed}

	got := MergeTraceTaskLogEvidence(
		domain.RunDir{Path: t.TempDir()},
		task,
		domain.LogOnlyTaskEvidence{Workdir: " \t", ErrorSummary: "\n"},
	)

	if got.Workdir != "" {
		t.Fatalf("merged Workdir = %q, want empty when evidence workdir is blank", got.Workdir)
	}
	if got.ErrorSummary != "" {
		t.Fatalf("merged ErrorSummary = %q, want empty when evidence summary is blank", got.ErrorSummary)
	}
}

func assertTraceOwnedTaskFieldsPreserved(t *testing.T, got domain.Task, want domain.Task) {
	t.Helper()
	if got.RowOrder != want.RowOrder || got.ID != want.ID || got.Status != want.Status || got.Process != want.Process || got.Name != want.Name || got.Tag != want.Tag || got.Duration != want.Duration || got.Realtime != want.Realtime || got.CPUs != want.CPUs || got.Memory != want.Memory {
		t.Fatalf("trace-owned fields = %#v, want preserved from %#v", got, want)
	}
	if got.Exit != want.Exit {
		t.Fatalf("merged Exit pointer = %#v, want original trace-owned pointer %#v", got.Exit, want.Exit)
	}
}

func assertLogEvidenceEqual(t *testing.T, got domain.LogOnlyTaskEvidence, want domain.LogOnlyTaskEvidence) {
	t.Helper()
	if got.ID != want.ID || got.Workdir != want.Workdir || got.Process != want.Process || got.Name != want.Name || got.ErrorSummary != want.ErrorSummary {
		t.Fatalf("evidence = %#v, want %#v", got, want)
	}
	if len(got.Sources) != len(want.Sources) {
		t.Fatalf("evidence sources = %#v, want %#v", got.Sources, want.Sources)
	}
	for index := range want.Sources {
		if got.Sources[index] != want.Sources[index] {
			t.Fatalf("evidence source[%d] = %#v, want %#v", index, got.Sources[index], want.Sources[index])
		}
	}
}

func assertTraceLogEvidenceDiagnostic(t *testing.T, got *domain.Diagnostic, wantSeverity domain.DiagnosticSeverity, wantCode string, wantMessage string) {
	t.Helper()
	if got == nil {
		t.Fatalf("diagnostic = nil, want %s/%s", wantSeverity, wantCode)
	}
	if got.Severity != wantSeverity || got.Code != wantCode || got.Message != wantMessage {
		t.Fatalf("diagnostic = %#v, want severity=%q code=%q message=%q", *got, wantSeverity, wantCode, wantMessage)
	}
}
