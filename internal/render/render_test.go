package render

import (
	"bytes"
	"encoding/json"
	"errors"
	"reflect"
	"strings"
	"testing"
	"time"

	"github.com/mskilab-org/gosh/internal/domain"
)

type failingWriter struct {
	err error
}

func (w failingWriter) Write([]byte) (int, error) {
	return 0, w.err
}

func TestBuildDiagnosticBlocksEmptyInput(t *testing.T) {
	if got := BuildDiagnosticBlocks(nil); len(got) != 0 {
		t.Fatalf("BuildDiagnosticBlocks(nil) length = %d, want 0", len(got))
	}
	if got := BuildDiagnosticBlocks([]domain.Diagnostic{}); len(got) != 0 {
		t.Fatalf("BuildDiagnosticBlocks(empty) length = %d, want 0", len(got))
	}
}

func TestBuildDiagnosticBlocksParsesContextDetailsAndHints(t *testing.T) {
	diagnostics := []domain.Diagnostic{
		{
			Severity: domain.DiagnosticWarning,
			Code:     "log_only_degraded",
			Message:  "complete task/resource/status data is unavailable without a Nextflow trace file",
			Detail: strings.Join([]string{
				"Selected log: /runs/example/.nextflow.log",
				"Observed log-only evidence rows: 2",
				"Complete task counts require a trace file.",
				"",
				"Use `nextflow run ... -with-trace` for future runs.",
			}, "\r\n"),
		},
	}

	got := BuildDiagnosticBlocks(diagnostics)
	want := []domain.DiagnosticBlock{
		{
			Severity: domain.DiagnosticWarning,
			Code:     "log_only_degraded",
			Title:    "complete task/resource/status data is unavailable without a Nextflow trace file",
			Context: []domain.DiagnosticContextLine{
				{Label: "Selected log", Value: "/runs/example/.nextflow.log"},
				{Label: "Observed log-only evidence rows", Value: "2"},
			},
			Details: []string{"Complete task counts require a trace file."},
			Hints:   []string{"Use `nextflow run ... -with-trace` for future runs."},
		},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("BuildDiagnosticBlocks() = %#v, want %#v", got, want)
	}
}

func TestBuildDiagnosticBlocksOrdersBySeverityAndGroupsRecommendationHints(t *testing.T) {
	diagnostics := []domain.Diagnostic{
		{
			Severity: domain.DiagnosticInfo,
			Code:     "nextflow_with_trace_recommended",
			Message:  "Run future Nextflow workflows with -with-trace",
			Detail:   "Use `nextflow run ... -with-trace` for future runs so gosh can build a complete trace-backed task index.",
		},
		{
			Severity: domain.DiagnosticWarning,
			Code:     "index_stale",
			Message:  "index is stale",
			Detail:   "Refresh with `gosh index --refresh`.",
		},
		{
			Severity: domain.DiagnosticError,
			Code:     "unsupported_artifacts",
			Message:  "No supported Nextflow trace or log artifacts found",
			Detail:   "Searched trace patterns: trace*.txt\nSearched log patterns: .nextflow.log",
		},
	}

	got := BuildDiagnosticBlocks(diagnostics)
	want := []domain.DiagnosticBlock{
		{
			Severity: domain.DiagnosticError,
			Code:     "unsupported_artifacts",
			Title:    "No supported Nextflow trace or log artifacts found",
			Context: []domain.DiagnosticContextLine{
				{Label: "Searched trace patterns", Value: "trace*.txt"},
				{Label: "Searched log patterns", Value: ".nextflow.log"},
			},
			Hints: []string{"Use `nextflow run ... -with-trace` for future runs so gosh can build a complete trace-backed task index."},
		},
		{
			Severity: domain.DiagnosticWarning,
			Code:     "index_stale",
			Title:    "index is stale",
			Hints:    []string{"Refresh with `gosh index --refresh`."},
		},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("BuildDiagnosticBlocks() = %#v, want %#v", got, want)
	}
}

func TestBuildDiagnosticBlocksFallsBackToCodeTitleAndDoesNotMutateDiagnostics(t *testing.T) {
	diagnostics := []domain.Diagnostic{
		{
			Severity: domain.DiagnosticError,
			Code:     "selector_not_found",
			Detail:   "Try `gosh tasks` to list available task IDs.",
		},
	}
	original := append([]domain.Diagnostic(nil), diagnostics...)

	got := BuildDiagnosticBlocks(diagnostics)
	want := []domain.DiagnosticBlock{
		{
			Severity: domain.DiagnosticError,
			Code:     "selector_not_found",
			Title:    "selector_not_found",
			Hints:    []string{"Try `gosh tasks` to list available task IDs."},
		},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("BuildDiagnosticBlocks() = %#v, want %#v", got, want)
	}
	if !reflect.DeepEqual(diagnostics, original) {
		t.Fatalf("BuildDiagnosticBlocks() mutated diagnostics to %#v, want %#v", diagnostics, original)
	}
}

func TestRenderDiagnosticBlocksHumanRendersSeverityBlocksContextDetailsAndHints(t *testing.T) {
	blocks := []domain.DiagnosticBlock{
		{
			Severity: domain.DiagnosticError,
			Code:     "unsupported_artifacts",
			Title:    "No supported Nextflow trace or log artifacts found",
			Context: []domain.DiagnosticContextLine{
				{Label: "Searched trace patterns", Value: "trace*.txt"},
				{Label: "Searched log patterns", Value: ".nextflow.log"},
			},
			Details: []string{"Complete task counts require a trace file."},
			Hints:   []string{"Use `nextflow run ... -with-trace` for future runs."},
		},
		{
			Severity: domain.DiagnosticWarning,
			Code:     "index_stale",
			Title:    "index is stale",
			Hints:    []string{"Refresh with `gosh index --refresh`."},
		},
		{
			Severity: domain.DiagnosticInfo,
			Code:     "selected_trace",
			Title:    "selected newest trace",
			Context:  []domain.DiagnosticContextLine{{Label: "trace", Value: "/runs/example/trace.txt"}},
		},
	}

	var buf bytes.Buffer
	if err := RenderDiagnosticBlocksHuman(&buf, blocks); err != nil {
		t.Fatalf("RenderDiagnosticBlocksHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"error: No supported Nextflow trace or log artifacts found",
		"  code: unsupported_artifacts",
		"  Searched trace patterns: trace*.txt",
		"  Searched log patterns: .nextflow.log",
		"  Complete task counts require a trace file.",
		"hint: Use `nextflow run ... -with-trace` for future runs.",
		"",
		"warning: index is stale",
		"  code: index_stale",
		"hint: Refresh with `gosh index --refresh`.",
		"",
		"info: selected newest trace",
		"  code: selected_trace",
		"  trace: /runs/example/trace.txt",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderDiagnosticBlocksHuman() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderDiagnosticBlocksHumanEmptyInputWritesNothing(t *testing.T) {
	var buf bytes.Buffer
	if err := RenderDiagnosticBlocksHuman(&buf, nil); err != nil {
		t.Fatalf("RenderDiagnosticBlocksHuman(nil) error = %v, want nil", err)
	}
	if got := buf.String(); got != "" {
		t.Fatalf("RenderDiagnosticBlocksHuman(nil) = %q, want empty output", got)
	}

	if err := RenderDiagnosticBlocksHuman(&buf, []domain.DiagnosticBlock{}); err != nil {
		t.Fatalf("RenderDiagnosticBlocksHuman(empty) error = %v, want nil", err)
	}
	if got := buf.String(); got != "" {
		t.Fatalf("RenderDiagnosticBlocksHuman(empty) = %q, want empty output", got)
	}
}

func TestRenderDiagnosticBlocksHumanFallsBackToCodeTitleWithoutDuplicateCodeLine(t *testing.T) {
	blocks := []domain.DiagnosticBlock{{Severity: domain.DiagnosticInfo, Code: "selector_not_found"}}

	var buf bytes.Buffer
	if err := RenderDiagnosticBlocksHuman(&buf, blocks); err != nil {
		t.Fatalf("RenderDiagnosticBlocksHuman() error = %v, want nil", err)
	}

	want := "info: selector_not_found\n"
	if got := buf.String(); got != want {
		t.Fatalf("RenderDiagnosticBlocksHuman() = %q, want %q", got, want)
	}
}

func TestRenderDiagnosticBlocksHumanNilWriter(t *testing.T) {
	err := RenderDiagnosticBlocksHuman(nil, nil)
	if err == nil {
		t.Fatal("RenderDiagnosticBlocksHuman() error = nil, want nil writer error")
	}
	if !strings.Contains(err.Error(), "nil writer") {
		t.Fatalf("RenderDiagnosticBlocksHuman() error = %v, want nil writer message", err)
	}
}

func TestRenderDiagnosticBlocksHumanReturnsWriterError(t *testing.T) {
	wantErr := errors.New("write failed")
	err := RenderDiagnosticBlocksHuman(failingWriter{err: wantErr}, []domain.DiagnosticBlock{
		{Severity: domain.DiagnosticError, Title: "No supported artifacts"},
	})
	if !errors.Is(err, wantErr) {
		t.Fatalf("RenderDiagnosticBlocksHuman() error = %v, want wrapping %v", err, wantErr)
	}
}

func TestRenderStatusHumanTraceBackedSummaryDeterministic(t *testing.T) {
	builtAt := time.Date(2024, 4, 28, 12, 34, 56, 0, time.UTC)
	traceMod := time.Date(2024, 4, 28, 12, 0, 0, 0, time.UTC)
	logMod := time.Date(2024, 4, 28, 12, 1, 0, 0, time.UTC)
	exit := 137
	view := domain.StatusView{
		Summary: domain.StatusSummary{
			RunDir:    domain.RunDir{Path: "/runs/example"},
			Mode:      domain.IndexModeTraceBacked,
			IndexPath: "/runs/example/.gosh/index.sqlite",
			Freshness: domain.IndexFreshnessFresh,
			BuiltAt:   &builtAt,
			Sources: domain.ArtifactSet{
				Trace: &domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: "/runs/example/trace.txt", ModTime: traceMod, Size: 1200},
				Log:   &domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/example/.nextflow.log", ModTime: logMod, Size: 3400},
			},
			Counts: []domain.StatusCount{
				{Status: domain.TaskStatusFailed, Count: 1},
				{Status: domain.TaskStatusCompleted, Count: 7},
			},
			FailedCount: 1,
			FailedPreview: []domain.FailedTaskPreview{
				{
					ID:           "aa/111111",
					Status:       domain.TaskStatusFailed,
					Process:      "ALIGN",
					Name:         "sample-1",
					Tag:          "tumor",
					Workdir:      "/runs/example/work/aa/111111",
					Exit:         &exit,
					ErrorSummary: "command exited with 137",
				},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderStatusHuman(&buf, view); err != nil {
		t.Fatalf("RenderStatusHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"run_dir: /runs/example",
		"mode: trace-backed",
		"index: /runs/example/.gosh/index.sqlite",
		"freshness: fresh",
		"built_at: 2024-04-28T12:34:56Z",
		"sources:",
		"  trace: /runs/example/trace.txt (mtime=2024-04-28T12:00:00Z size=1200)",
		"  log: /runs/example/.nextflow.log (mtime=2024-04-28T12:01:00Z size=3400)",
		"counts:",
		"  COMPLETED: 7",
		"  FAILED: 1",
		"failed_count: 1",
		"failed_preview:",
		"  - id=aa/111111 status=FAILED process=ALIGN name=sample-1 tag=tumor workdir=/runs/example/work/aa/111111 exit=137 error=command exited with 137",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderStatusHuman() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderStatusHumanTraceBackedNoFailures(t *testing.T) {
	traceMod := time.Date(2024, 4, 28, 13, 0, 0, 0, time.UTC)
	view := domain.StatusView{
		Summary: domain.StatusSummary{
			RunDir:    domain.RunDir{Path: "/runs/ok"},
			Mode:      domain.IndexModeTraceBacked,
			IndexPath: "/runs/ok/.gosh/index.sqlite",
			Freshness: domain.IndexFreshnessFresh,
			Sources: domain.ArtifactSet{
				Trace: &domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: "/runs/ok/trace.tsv", ModTime: traceMod, Size: 100},
			},
			Counts: []domain.StatusCount{
				{Status: domain.TaskStatusCached, Count: 2},
				{Status: domain.TaskStatusCompleted, Count: 5},
			},
			FailedCount: 0,
		},
	}

	var buf bytes.Buffer
	if err := RenderStatusHuman(&buf, view); err != nil {
		t.Fatalf("RenderStatusHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"run_dir: /runs/ok",
		"mode: trace-backed",
		"index: /runs/ok/.gosh/index.sqlite",
		"freshness: fresh",
		"sources:",
		"  trace: /runs/ok/trace.tsv (mtime=2024-04-28T13:00:00Z size=100)",
		"  log: none",
		"counts:",
		"  CACHED: 2",
		"  COMPLETED: 5",
		"failed_count: 0",
		"failed_preview: none",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderStatusHuman() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderStatusHumanTraceBackedDiagnosticsUseGitStyleBlocks(t *testing.T) {
	view := domain.StatusView{
		Summary: domain.StatusSummary{
			RunDir:      domain.RunDir{Path: "/runs/stale"},
			Mode:        domain.IndexModeTraceBacked,
			IndexPath:   "/runs/stale/.gosh/index.sqlite",
			Freshness:   domain.IndexFreshnessStale,
			FailedCount: 0,
			Diagnostics: []domain.Diagnostic{
				{
					Severity: domain.DiagnosticWarning,
					Code:     "index_stale",
					Message:  "index is stale",
					Detail:   "Index path: /runs/stale/.gosh/index.sqlite\nRefresh with `gosh index --refresh`.",
				},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderStatusHuman(&buf, view); err != nil {
		t.Fatalf("RenderStatusHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"run_dir: /runs/stale",
		"mode: trace-backed",
		"index: /runs/stale/.gosh/index.sqlite",
		"freshness: stale",
		"sources:",
		"  trace: none",
		"  log: none",
		"counts: none",
		"failed_count: 0",
		"failed_preview: none",
		"",
		"warning: index is stale",
		"  code: index_stale",
		"  Index path: /runs/stale/.gosh/index.sqlite",
		"hint: Refresh with `gosh index --refresh`.",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderStatusHuman() =\n%s\nwant\n%s", got, want)
	}
	if strings.Contains(buf.String(), "diagnostics:") {
		t.Fatalf("RenderStatusHuman() output = %q, did not want nested diagnostics heading", buf.String())
	}
}

func TestRenderStatusHumanLogOnlyIncludesEvidenceAndGitStyleDiagnostics(t *testing.T) {
	logMod := time.Date(2024, 4, 28, 14, 0, 0, 0, time.UTC)
	exit := 2
	view := domain.StatusView{
		Summary: domain.StatusSummary{
			RunDir:    domain.RunDir{Path: "/runs/log-only"},
			Mode:      domain.IndexModeLogOnly,
			Freshness: domain.IndexFreshnessUnsupported,
			Sources: domain.ArtifactSet{
				Log:              &domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/log-only/.nextflow.log", ModTime: logMod, Size: 900},
				SearchedPatterns: []string{"trace*.txt", "trace*.csv", ".nextflow*.log"},
			},
			Counts: []domain.StatusCount{
				{Status: domain.TaskStatusCompleted, Count: 1},
				{Status: domain.TaskStatusFailed, Count: 1},
			},
			FailedCount: 1,
			LogOnlyEvidence: []domain.LogOnlyTaskEvidence{
				{
					ID:             "aa/111111",
					Process:        "PIPE:QC",
					Name:           "normal-01",
					ObservedStatus: domain.TaskStatusCompleted,
					Completeness:   domain.LogOnlyEvidencePartial,
				},
				{
					ID:                    "bb/222222",
					Process:               "PIPE:CALL",
					Name:                  "tumor-02",
					Workdir:               "/runs/log-only/work/bb/222222",
					ObservedStatus:        domain.TaskStatusFailed,
					Exit:                  &exit,
					ErrorSummary:          "No such file or directory",
					Completeness:          domain.LogOnlyEvidencePartial,
					CommandFilesAvailable: true,
				},
			},
			Diagnostics: []domain.Diagnostic{
				{
					Severity: domain.DiagnosticWarning,
					Code:     "log_only_degraded",
					Message:  "log-only status is degraded; complete task/resource/status data is unavailable",
					Detail:   "Selected log: /runs/log-only/.nextflow.log\nComplete task counts require a trace file.",
				},
				{
					Severity: domain.DiagnosticInfo,
					Code:     "nextflow_with_trace_recommended",
					Message:  "Run future Nextflow workflows with -with-trace",
				},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderStatusHuman(&buf, view); err != nil {
		t.Fatalf("RenderStatusHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"run_dir: /runs/log-only",
		"mode: log-only",
		"freshness: unsupported",
		"sources:",
		"  trace: none",
		"  log: /runs/log-only/.nextflow.log (mtime=2024-04-28T14:00:00Z size=900)",
		"  searched_patterns: trace*.txt, trace*.csv, .nextflow*.log",
		"counts: incomplete (log-only evidence; trace file required)",
		"observed_counts:",
		"  COMPLETED: 1",
		"  FAILED: 1",
		"failed_count: 1",
		"log_only_evidence:",
		"  - id=aa/111111 status=COMPLETED process=PIPE:QC name=normal-01 workdir=- exit=- completeness=partial command_files_available=false error=-",
		"  - id=bb/222222 status=FAILED process=PIPE:CALL name=tumor-02 workdir=/runs/log-only/work/bb/222222 exit=2 completeness=partial command_files_available=true error=No such file or directory",
		"",
		"warning: log-only status is degraded; complete task/resource/status data is unavailable",
		"  code: log_only_degraded",
		"  Selected log: /runs/log-only/.nextflow.log",
		"  Complete task counts require a trace file.",
		"hint: Run future Nextflow workflows with -with-trace",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderStatusHuman() =\n%s\nwant\n%s", got, want)
	}
	if strings.Contains(buf.String(), "diagnostics:") {
		t.Fatalf("RenderStatusHuman() output = %q, did not want nested diagnostics heading", buf.String())
	}
}

func TestRenderStatusJSONTraceBackedSummaryStable(t *testing.T) {
	builtAt := time.Date(2024, 4, 28, 12, 34, 56, 0, time.UTC)
	selectedAt := time.Date(2024, 4, 28, 11, 59, 0, 0, time.UTC)
	traceMod := time.Date(2024, 4, 28, 12, 0, 0, 0, time.UTC)
	logMod := time.Date(2024, 4, 28, 12, 1, 0, 0, time.UTC)
	exit := 137
	view := domain.StatusView{
		Format: domain.OutputFormatJSON,
		Summary: domain.StatusSummary{
			RunDir:    domain.RunDir{Path: "/runs/example"},
			Mode:      domain.IndexModeTraceBacked,
			IndexPath: "/runs/example/.gosh/index.sqlite",
			Freshness: domain.IndexFreshnessFresh,
			BuiltAt:   &builtAt,
			Sources: domain.ArtifactSet{
				RunDir:           domain.RunDir{Path: "/runs/example"},
				Mode:             domain.IndexModeTraceBacked,
				Trace:            &domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: "/runs/example/trace.txt", ModTime: traceMod, Size: 1200},
				Log:              &domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/example/.nextflow.log", ModTime: logMod, Size: 3400},
				SelectedAt:       selectedAt,
				SearchedPatterns: []string{"trace*.txt", ".nextflow*.log"},
				Diagnostics: []domain.Diagnostic{
					{Severity: domain.DiagnosticInfo, Code: "selected_trace", Message: "selected newest trace", Detail: "trace detail"},
				},
			},
			Counts: []domain.StatusCount{
				{Status: domain.TaskStatusFailed, Count: 1},
				{Status: domain.TaskStatusCompleted, Count: 7},
			},
			FailedCount: 1,
			FailedPreview: []domain.FailedTaskPreview{
				{
					ID:           "aa/111111",
					Status:       domain.TaskStatusFailed,
					Process:      "ALIGN",
					Name:         "sample-1",
					Tag:          "tumor",
					Workdir:      "/runs/example/work/aa/111111",
					Exit:         &exit,
					ErrorSummary: "command exited with 137",
				},
			},
			Diagnostics: []domain.Diagnostic{
				{Severity: domain.DiagnosticWarning, Code: "stale_ignored", Message: "index was refreshed", Detail: "old index replaced"},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderStatusJSON(&buf, view); err != nil {
		t.Fatalf("RenderStatusJSON() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"{",
		`  "format": "json",`,
		`  "summary": {`,
		`    "run_dir": "/runs/example",`,
		`    "mode": "trace-backed",`,
		`    "index_path": "/runs/example/.gosh/index.sqlite",`,
		`    "freshness": "fresh",`,
		`    "built_at": "2024-04-28T12:34:56Z",`,
		`    "sources": {`,
		`      "run_dir": "/runs/example",`,
		`      "mode": "trace-backed",`,
		`      "trace": {`,
		`        "kind": "trace",`,
		`        "path": "/runs/example/trace.txt",`,
		`        "mod_time": "2024-04-28T12:00:00Z",`,
		`        "size": 1200`,
		`      },`,
		`      "log": {`,
		`        "kind": "log",`,
		`        "path": "/runs/example/.nextflow.log",`,
		`        "mod_time": "2024-04-28T12:01:00Z",`,
		`        "size": 3400`,
		`      },`,
		`      "selected_at": "2024-04-28T11:59:00Z",`,
		`      "searched_patterns": [`,
		`        "trace*.txt",`,
		`        ".nextflow*.log"`,
		`      ],`,
		`      "diagnostics": [`,
		`        {`,
		`          "severity": "info",`,
		`          "code": "selected_trace",`,
		`          "message": "selected newest trace",`,
		`          "detail": "trace detail"`,
		`        }`,
		`      ]`,
		`    },`,
		`    "counts": [`,
		`      {`,
		`        "status": "COMPLETED",`,
		`        "count": 7`,
		`      },`,
		`      {`,
		`        "status": "FAILED",`,
		`        "count": 1`,
		`      }`,
		`    ],`,
		`    "failed_count": 1,`,
		`    "failed_preview": [`,
		`      {`,
		`        "id": "aa/111111",`,
		`        "status": "FAILED",`,
		`        "process": "ALIGN",`,
		`        "name": "sample-1",`,
		`        "tag": "tumor",`,
		`        "workdir": "/runs/example/work/aa/111111",`,
		`        "exit": 137,`,
		`        "error_summary": "command exited with 137"`,
		`      }`,
		`    ],`,
		`    "log_only_failures": [],`,
		`    "log_only_evidence": [],`,
		`    "diagnostics": [`,
		`      {`,
		`        "severity": "warning",`,
		`        "code": "stale_ignored",`,
		`        "message": "index was refreshed",`,
		`        "detail": "old index replaced"`,
		`      }`,
		`    ]`,
		`  }`,
		"}",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderStatusJSON() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderStatusJSONLogOnlyIncludesFailureEvidenceAndEmptyArrays(t *testing.T) {
	logMod := time.Date(2024, 4, 28, 14, 0, 0, 0, time.UTC)
	exit := 2
	view := domain.StatusView{
		Summary: domain.StatusSummary{
			RunDir:    domain.RunDir{Path: "/runs/log-only"},
			Mode:      domain.IndexModeLogOnly,
			Freshness: domain.IndexFreshnessUnsupported,
			Sources: domain.ArtifactSet{
				RunDir:           domain.RunDir{Path: "/runs/log-only"},
				Mode:             domain.IndexModeLogOnly,
				Log:              &domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/log-only/.nextflow.log", ModTime: logMod, Size: 900},
				SearchedPatterns: []string{"trace*.txt", "trace*.csv", ".nextflow*.log"},
			},
			FailedCount: 1,
			LogOnlyFailures: []domain.LogOnlyFailure{
				{
					ID:           "bb/222222",
					Workdir:      "/runs/log-only/work/bb/222222",
					Process:      "PIPE:CALL",
					Name:         "tumor-02",
					Exit:         &exit,
					ErrorSummary: "No such file or directory",
					ErrorBlock:   "ERROR ~ Process PIPE:CALL (tumor-02) failed",
				},
			},
			Diagnostics: []domain.Diagnostic{
				{
					Severity: domain.DiagnosticWarning,
					Code:     "log_only_degraded",
					Message:  "log-only status is degraded; complete task/resource/status data is unavailable",
					Detail:   "Selected log: /runs/log-only/.nextflow.log\nComplete task counts require a trace file.",
				},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderStatusJSON(&buf, view); err != nil {
		t.Fatalf("RenderStatusJSON() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"{",
		`  "format": "json",`,
		`  "summary": {`,
		`    "run_dir": "/runs/log-only",`,
		`    "mode": "log-only",`,
		`    "index_path": "",`,
		`    "freshness": "unsupported",`,
		`    "built_at": null,`,
		`    "sources": {`,
		`      "run_dir": "/runs/log-only",`,
		`      "mode": "log-only",`,
		`      "trace": null,`,
		`      "log": {`,
		`        "kind": "log",`,
		`        "path": "/runs/log-only/.nextflow.log",`,
		`        "mod_time": "2024-04-28T14:00:00Z",`,
		`        "size": 900`,
		`      },`,
		`      "selected_at": null,`,
		`      "searched_patterns": [`,
		`        "trace*.txt",`,
		`        "trace*.csv",`,
		`        ".nextflow*.log"`,
		`      ],`,
		`      "diagnostics": []`,
		`    },`,
		`    "counts": [],`,
		`    "failed_count": 1,`,
		`    "failed_preview": [],`,
		`    "log_only_failures": [`,
		`      {`,
		`        "id": "bb/222222",`,
		`        "workdir": "/runs/log-only/work/bb/222222",`,
		`        "process": "PIPE:CALL",`,
		`        "name": "tumor-02",`,
		`        "exit": 2,`,
		`        "error_summary": "No such file or directory",`,
		`        "error_block": "ERROR ~ Process PIPE:CALL (tumor-02) failed"`,
		`      }`,
		`    ],`,
		`    "log_only_evidence": [],`,
		`    "diagnostics": [`,
		`      {`,
		`        "severity": "warning",`,
		`        "code": "log_only_degraded",`,
		`        "message": "log-only status is degraded; complete task/resource/status data is unavailable",`,
		`        "detail": "Selected log: /runs/log-only/.nextflow.log\nComplete task counts require a trace file."`,
		`      }`,
		`    ]`,
		`  }`,
		"}",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderStatusJSON() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderStatusJSONLogOnlyIncludesObservedEvidenceAndPreservesDiagnostics(t *testing.T) {
	logMod := time.Date(2024, 4, 28, 14, 0, 0, 0, time.UTC)
	exit := 2
	logPath := "/runs/log-only/.nextflow.log"
	view := domain.StatusView{
		Summary: domain.StatusSummary{
			RunDir:    domain.RunDir{Path: "/runs/log-only"},
			Mode:      domain.IndexModeLogOnly,
			Freshness: domain.IndexFreshnessUnsupported,
			Sources: domain.ArtifactSet{
				RunDir: domain.RunDir{Path: "/runs/log-only"},
				Mode:   domain.IndexModeLogOnly,
				Log:    &domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath, ModTime: logMod, Size: 900},
				Diagnostics: []domain.Diagnostic{
					{Severity: domain.DiagnosticInfo, Code: "selected_log", Message: "selected newest Nextflow log", Detail: "log source detail"},
				},
			},
			Counts: []domain.StatusCount{
				{Status: domain.TaskStatusFailed, Count: 1},
				{Status: domain.TaskStatusCompleted, Count: 1},
			},
			FailedCount: 1,
			LogOnlyEvidence: []domain.LogOnlyTaskEvidence{
				{
					ID:             "aa/111111",
					Process:        "PIPE:QC",
					Name:           "normal-01",
					ObservedStatus: domain.TaskStatusCompleted,
					Sources: []domain.LogOnlyEvidenceSource{
						{Kind: domain.LogOnlyEvidenceSourceLog, Path: logPath, Detail: "lifecycle line"},
					},
					Completeness: domain.LogOnlyEvidencePartial,
				},
				{
					ID:             "bb/222222",
					Workdir:        "/runs/log-only/work/bb/222222",
					Process:        "PIPE:CALL",
					Name:           "tumor-02",
					ObservedStatus: domain.TaskStatusFailed,
					Exit:           &exit,
					ErrorSummary:   "No such file or directory",
					ErrorBlock:     "ERROR ~ Process PIPE:CALL (tumor-02) failed",
					Sources: []domain.LogOnlyEvidenceSource{
						{Kind: domain.LogOnlyEvidenceSourceLog, Path: logPath, Detail: "failure block"},
						{Kind: domain.LogOnlyEvidenceSourceCommand, Path: "/runs/log-only/work/bb/222222/.command.err", Detail: "stderr snippet"},
					},
					Completeness:          domain.LogOnlyEvidencePartial,
					CommandFilesAvailable: true,
				},
			},
			Diagnostics: []domain.Diagnostic{
				{
					Severity: domain.DiagnosticWarning,
					Code:     "log_only_degraded",
					Message:  "log-only status is degraded; complete task/resource/status data is unavailable",
					Detail:   "Selected log: /runs/log-only/.nextflow.log\nComplete task counts require a trace file.",
				},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderStatusJSON(&buf, view); err != nil {
		t.Fatalf("RenderStatusJSON() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"{",
		`  "format": "json",`,
		`  "summary": {`,
		`    "run_dir": "/runs/log-only",`,
		`    "mode": "log-only",`,
		`    "index_path": "",`,
		`    "freshness": "unsupported",`,
		`    "built_at": null,`,
		`    "sources": {`,
		`      "run_dir": "/runs/log-only",`,
		`      "mode": "log-only",`,
		`      "trace": null,`,
		`      "log": {`,
		`        "kind": "log",`,
		`        "path": "/runs/log-only/.nextflow.log",`,
		`        "mod_time": "2024-04-28T14:00:00Z",`,
		`        "size": 900`,
		`      },`,
		`      "selected_at": null,`,
		`      "searched_patterns": [],`,
		`      "diagnostics": [`,
		`        {`,
		`          "severity": "info",`,
		`          "code": "selected_log",`,
		`          "message": "selected newest Nextflow log",`,
		`          "detail": "log source detail"`,
		`        }`,
		`      ]`,
		`    },`,
		`    "counts": [`,
		`      {`,
		`        "status": "COMPLETED",`,
		`        "count": 1`,
		`      },`,
		`      {`,
		`        "status": "FAILED",`,
		`        "count": 1`,
		`      }`,
		`    ],`,
		`    "failed_count": 1,`,
		`    "failed_preview": [],`,
		`    "log_only_failures": [],`,
		`    "log_only_evidence": [`,
		`      {`,
		`        "id": "aa/111111",`,
		`        "workdir": "",`,
		`        "process": "PIPE:QC",`,
		`        "name": "normal-01",`,
		`        "observed_status": "COMPLETED",`,
		`        "exit": null,`,
		`        "error_summary": "",`,
		`        "error_block": "",`,
		`        "sources": [`,
		`          {`,
		`            "kind": "log",`,
		`            "path": "/runs/log-only/.nextflow.log",`,
		`            "detail": "lifecycle line"`,
		`          }`,
		`        ],`,
		`        "completeness": "partial",`,
		`        "command_files_available": false`,
		`      },`,
		`      {`,
		`        "id": "bb/222222",`,
		`        "workdir": "/runs/log-only/work/bb/222222",`,
		`        "process": "PIPE:CALL",`,
		`        "name": "tumor-02",`,
		`        "observed_status": "FAILED",`,
		`        "exit": 2,`,
		`        "error_summary": "No such file or directory",`,
		`        "error_block": "ERROR ~ Process PIPE:CALL (tumor-02) failed",`,
		`        "sources": [`,
		`          {`,
		`            "kind": "log",`,
		`            "path": "/runs/log-only/.nextflow.log",`,
		`            "detail": "failure block"`,
		`          },`,
		`          {`,
		`            "kind": "command-file",`,
		`            "path": "/runs/log-only/work/bb/222222/.command.err",`,
		`            "detail": "stderr snippet"`,
		`          }`,
		`        ],`,
		`        "completeness": "partial",`,
		`        "command_files_available": true`,
		`      }`,
		`    ],`,
		`    "diagnostics": [`,
		`      {`,
		`        "severity": "warning",`,
		`        "code": "log_only_degraded",`,
		`        "message": "log-only status is degraded; complete task/resource/status data is unavailable",`,
		`        "detail": "Selected log: /runs/log-only/.nextflow.log\nComplete task counts require a trace file."`,
		`      }`,
		`    ]`,
		`  }`,
		"}",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderStatusJSON() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderIndexHumanShowsSearchLocationsAndSeveritySortedDiagnostics(t *testing.T) {
	builtAt := time.Date(2024, 4, 28, 12, 34, 56, 0, time.UTC)
	traceMod := time.Date(2024, 4, 28, 12, 0, 0, 0, time.UTC)
	logMod := time.Date(2024, 4, 28, 12, 1, 0, 0, time.UTC)
	view := domain.IndexView{
		Diagnostics: domain.IndexDiagnostics{
			RunDir: domain.RunDir{Path: "/runs/example"},
			Artifacts: domain.ArtifactSet{
				RunDir:           domain.RunDir{Path: "/runs/example"},
				Mode:             domain.IndexModeTraceBacked,
				Trace:            &domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: "/runs/example/trace.txt", ModTime: traceMod, Size: 1200},
				Log:              &domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/example/.nextflow.log", ModTime: logMod, Size: 3400},
				SearchedPatterns: []string{"trace*.txt", ".nextflow.log"},
				SearchLocations: []domain.ArtifactSearchLocation{
					{Kind: domain.SourceKindTrace, BaseDir: "/runs/example", Patterns: []string{"trace*.txt", "trace*.csv"}, Description: "run directory trace files"},
					{Kind: domain.SourceKindTrace, BaseDir: "/runs/example/results/pipeline_info", Patterns: []string{"execution_trace*.txt"}, Description: "pipeline_info execution trace files"},
					{Kind: domain.SourceKindLog, BaseDir: "/runs/example", Patterns: []string{".nextflow.log"}, Description: "run directory log files"},
				},
				Diagnostics: []domain.Diagnostic{
					{Severity: domain.DiagnosticInfo, Code: "source_selected", Message: "selected newest trace"},
				},
			},
			Metadata: &domain.IndexMetadata{
				SchemaVersion: 1,
				RunDir:        "/runs/example",
				IndexPath:     "/runs/example/.gosh/index.sqlite",
				Mode:          domain.IndexModeTraceBacked,
				BuiltAt:       builtAt,
				Freshness:     domain.IndexFreshnessStale,
				StaleReason:   "selected trace changed",
				TaskCount:     42,
			},
			Diagnostics: []domain.Diagnostic{
				{Severity: domain.DiagnosticWarning, Code: "index_stale", Message: "index is stale", Detail: "Refresh with `gosh index --refresh`."},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderIndexHuman(&buf, view); err != nil {
		t.Fatalf("RenderIndexHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"run_dir: /runs/example",
		"mode: trace-backed",
		"index: /runs/example/.gosh/index.sqlite",
		"freshness: stale",
		"stale_reason: selected trace changed",
		"built_at: 2024-04-28T12:34:56Z",
		"task_count: 42",
		"sources:",
		"  trace: /runs/example/trace.txt (mtime=2024-04-28T12:00:00Z size=1200)",
		"  log: /runs/example/.nextflow.log (mtime=2024-04-28T12:01:00Z size=3400)",
		"  searched_locations:",
		"    - trace: /runs/example (run directory trace files; patterns: trace*.txt, trace*.csv)",
		"    - trace: /runs/example/results/pipeline_info (pipeline_info execution trace files; patterns: execution_trace*.txt)",
		"    - log: /runs/example (run directory log files; patterns: .nextflow.log)",
		"",
		"warning: index is stale",
		"  code: index_stale",
		"hint: Refresh with `gosh index --refresh`.",
		"",
		"info: selected newest trace",
		"  code: source_selected",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderIndexHuman() =\n%s\nwant\n%s", got, want)
	}
	if strings.Contains(buf.String(), "diagnostics:") {
		t.Fatalf("RenderIndexHuman() output = %q, did not want nested diagnostics heading", buf.String())
	}
}

func TestRenderIndexHumanUnsupportedNoSourcesShowsZeroTaskCount(t *testing.T) {
	view := domain.IndexView{
		Diagnostics: domain.IndexDiagnostics{
			RunDir: domain.RunDir{Path: "/runs/empty"},
			Artifacts: domain.ArtifactSet{
				RunDir: domain.RunDir{Path: "/runs/empty"},
				Mode:   domain.IndexModeUnsupported,
				SearchLocations: []domain.ArtifactSearchLocation{
					{Kind: domain.SourceKindTrace, BaseDir: "/runs/empty", Patterns: []string{"trace*.txt"}, Description: "run directory trace files"},
					{Kind: domain.SourceKindTrace, BaseDir: "/runs/empty/results/pipeline_info", Patterns: []string{"execution_trace*.txt"}, Description: "pipeline_info execution trace files"},
					{Kind: domain.SourceKindLog, BaseDir: "/runs/empty", Patterns: []string{".nextflow.log"}, Description: "run directory log files"},
				},
				Diagnostics: []domain.Diagnostic{
					{Severity: domain.DiagnosticInfo, Code: "nextflow_with_trace_recommended", Message: "Run future Nextflow workflows with -with-trace", Detail: "Use `nextflow run ... -with-trace` for future runs so gosh can build a complete trace-backed task index."},
				},
			},
			Metadata: &domain.IndexMetadata{
				Mode:        domain.IndexModeUnsupported,
				Freshness:   domain.IndexFreshnessUnsupported,
				TaskCount:   0,
				StaleReason: "no supported artifacts",
			},
			Diagnostics: []domain.Diagnostic{
				{Severity: domain.DiagnosticError, Code: "unsupported_artifacts", Message: "No supported Nextflow trace or log artifacts found", Detail: "Searched trace patterns: trace*.txt\nSearched log patterns: .nextflow.log"},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderIndexHuman(&buf, view); err != nil {
		t.Fatalf("RenderIndexHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"run_dir: /runs/empty",
		"mode: unsupported",
		"freshness: unsupported",
		"stale_reason: no supported artifacts",
		"task_count: 0",
		"sources:",
		"  trace: none",
		"  log: none",
		"  searched_locations:",
		"    - trace: /runs/empty (run directory trace files; patterns: trace*.txt)",
		"    - trace: /runs/empty/results/pipeline_info (pipeline_info execution trace files; patterns: execution_trace*.txt)",
		"    - log: /runs/empty (run directory log files; patterns: .nextflow.log)",
		"",
		"error: No supported Nextflow trace or log artifacts found",
		"  code: unsupported_artifacts",
		"  Searched trace patterns: trace*.txt",
		"  Searched log patterns: .nextflow.log",
		"hint: Use `nextflow run ... -with-trace` for future runs so gosh can build a complete trace-backed task index.",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderIndexHuman() =\n%s\nwant\n%s", got, want)
	}
	if strings.Contains(buf.String(), "diagnostics:") {
		t.Fatalf("RenderIndexHuman() output = %q, did not want nested diagnostics heading", buf.String())
	}
}

func TestRenderIndexJSONTraceBackedStableIncludesArtifactsMetadataAndDiagnostics(t *testing.T) {
	builtAt := time.Date(2024, 4, 28, 12, 34, 56, 0, time.UTC)
	selectedAt := time.Date(2024, 4, 28, 11, 59, 0, 0, time.UTC)
	traceMod := time.Date(2024, 4, 28, 12, 0, 0, 0, time.UTC)
	logMod := time.Date(2024, 4, 28, 12, 1, 0, 0, time.UTC)
	view := domain.IndexView{
		Format: domain.OutputFormatJSON,
		Diagnostics: domain.IndexDiagnostics{
			RunDir: domain.RunDir{Path: "/runs/example"},
			Artifacts: domain.ArtifactSet{
				RunDir:           domain.RunDir{Path: "/runs/example"},
				Mode:             domain.IndexModeTraceBacked,
				Trace:            &domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: "/runs/example/trace.txt", ModTime: traceMod, Size: 1200},
				Log:              &domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/example/.nextflow.log", ModTime: logMod, Size: 3400},
				SelectedAt:       selectedAt,
				SearchedPatterns: []string{"trace*.txt", ".nextflow.log"},
				Diagnostics: []domain.Diagnostic{
					{Severity: domain.DiagnosticInfo, Code: "source_selected", Message: "selected newest trace", Detail: "trace detail"},
				},
			},
			Metadata: &domain.IndexMetadata{
				SchemaVersion: 1,
				RunDir:        "/runs/example",
				IndexPath:     "/runs/example/.gosh/index.sqlite",
				Mode:          domain.IndexModeTraceBacked,
				Trace:         &domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: "/runs/example/trace.txt", ModTime: traceMod, Size: 1200},
				Log:           &domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/example/.nextflow.log", ModTime: logMod, Size: 3400},
				BuiltAt:       builtAt,
				Freshness:     domain.IndexFreshnessStale,
				StaleReason:   "selected trace changed",
				TaskCount:     42,
			},
			Diagnostics: []domain.Diagnostic{
				{Severity: domain.DiagnosticWarning, Code: "index_stale", Message: "index is stale", Detail: "Refresh with `gosh index --refresh`."},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderIndexJSON(&buf, view); err != nil {
		t.Fatalf("RenderIndexJSON() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"{",
		`  "format": "json",`,
		`  "run_dir": "/runs/example",`,
		`  "artifacts": {`,
		`    "run_dir": "/runs/example",`,
		`    "mode": "trace-backed",`,
		`    "trace": {`,
		`      "kind": "trace",`,
		`      "path": "/runs/example/trace.txt",`,
		`      "mod_time": "2024-04-28T12:00:00Z",`,
		`      "size": 1200`,
		`    },`,
		`    "log": {`,
		`      "kind": "log",`,
		`      "path": "/runs/example/.nextflow.log",`,
		`      "mod_time": "2024-04-28T12:01:00Z",`,
		`      "size": 3400`,
		`    },`,
		`    "selected_at": "2024-04-28T11:59:00Z",`,
		`    "searched_patterns": [`,
		`      "trace*.txt",`,
		`      ".nextflow.log"`,
		`    ],`,
		`    "diagnostics": [`,
		`      {`,
		`        "severity": "info",`,
		`        "code": "source_selected",`,
		`        "message": "selected newest trace",`,
		`        "detail": "trace detail"`,
		`      }`,
		`    ]`,
		`  },`,
		`  "metadata": {`,
		`    "schema_version": 1,`,
		`    "run_dir": "/runs/example",`,
		`    "index_path": "/runs/example/.gosh/index.sqlite",`,
		`    "mode": "trace-backed",`,
		`    "trace": {`,
		`      "kind": "trace",`,
		`      "path": "/runs/example/trace.txt",`,
		`      "mod_time": "2024-04-28T12:00:00Z",`,
		`      "size": 1200`,
		`    },`,
		`    "log": {`,
		`      "kind": "log",`,
		`      "path": "/runs/example/.nextflow.log",`,
		`      "mod_time": "2024-04-28T12:01:00Z",`,
		`      "size": 3400`,
		`    },`,
		`    "built_at": "2024-04-28T12:34:56Z",`,
		`    "freshness": "stale",`,
		`    "stale_reason": "selected trace changed",`,
		`    "task_count": 42`,
		`  },`,
		`  "diagnostics": [`,
		`    {`,
		`      "severity": "warning",`,
		`      "code": "index_stale",`,
		`      "message": "index is stale",`,
		`      "detail": "Refresh with ` + "`" + `gosh index --refresh` + "`" + `."`,
		`    }`,
		`  ]`,
		"}",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderIndexJSON() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderIndexJSONUnsupportedWithoutMetadataIncludesNullsAndEmptyArrays(t *testing.T) {
	view := domain.IndexView{
		Diagnostics: domain.IndexDiagnostics{
			RunDir: domain.RunDir{Path: "/runs/empty"},
			Artifacts: domain.ArtifactSet{
				RunDir: domain.RunDir{Path: "/runs/empty"},
				Mode:   domain.IndexModeUnsupported,
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderIndexJSON(&buf, view); err != nil {
		t.Fatalf("RenderIndexJSON() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"{",
		`  "format": "json",`,
		`  "run_dir": "/runs/empty",`,
		`  "artifacts": {`,
		`    "run_dir": "/runs/empty",`,
		`    "mode": "unsupported",`,
		`    "trace": null,`,
		`    "log": null,`,
		`    "selected_at": null,`,
		`    "searched_patterns": [],`,
		`    "diagnostics": []`,
		`  },`,
		`  "metadata": null,`,
		`  "diagnostics": []`,
		"}",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderIndexJSON() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderUnsupportedDiagnosticsHumanShowsGitStyleBlocks(t *testing.T) {
	diagnostics := []domain.Diagnostic{
		{
			Severity: domain.DiagnosticInfo,
			Code:     "nextflow_with_trace_recommended",
			Message:  "Re-run with `-with-trace` to enable trace-backed summaries.",
		},
		{
			Severity: domain.DiagnosticError,
			Code:     "unsupported_artifacts",
			Message:  "No supported Nextflow trace or log artifacts found",
			Detail:   "Searched trace patterns: trace*.txt\nSearched log patterns: .nextflow.log",
		},
	}

	var buf bytes.Buffer
	if err := RenderUnsupportedDiagnostics(&buf, diagnostics, domain.OutputFormatHuman); err != nil {
		t.Fatalf("RenderUnsupportedDiagnostics() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"error: No supported Nextflow trace or log artifacts found",
		"  code: unsupported_artifacts",
		"  Searched trace patterns: trace*.txt",
		"  Searched log patterns: .nextflow.log",
		"hint: Re-run with `-with-trace` to enable trace-backed summaries.",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderUnsupportedDiagnostics() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderUnsupportedDiagnosticsJSONStable(t *testing.T) {
	diagnostics := []domain.Diagnostic{
		{Severity: domain.DiagnosticWarning, Code: "log_only_degraded", Message: "log-only status is degraded", Detail: "trace file required"},
		{Severity: domain.DiagnosticError, Code: "unsupported_artifacts", Message: "No supported Nextflow trace or log artifacts found"},
	}

	var buf bytes.Buffer
	if err := RenderUnsupportedDiagnostics(&buf, diagnostics, domain.OutputFormatJSON); err != nil {
		t.Fatalf("RenderUnsupportedDiagnostics() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"{",
		`  "format": "json",`,
		`  "diagnostics": [`,
		`    {`,
		`      "severity": "warning",`,
		`      "code": "log_only_degraded",`,
		`      "message": "log-only status is degraded",`,
		`      "detail": "trace file required"`,
		`    },`,
		`    {`,
		`      "severity": "error",`,
		`      "code": "unsupported_artifacts",`,
		`      "message": "No supported Nextflow trace or log artifacts found",`,
		`      "detail": ""`,
		`    }`,
		`  ]`,
		"}",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderUnsupportedDiagnostics() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderUnsupportedDiagnosticsHumanEmptyDiagnosticsWritesNothing(t *testing.T) {
	var buf bytes.Buffer
	if err := RenderUnsupportedDiagnostics(&buf, nil, domain.OutputFormatHuman); err != nil {
		t.Fatalf("RenderUnsupportedDiagnostics() error = %v, want nil", err)
	}

	if got := buf.String(); got != "" {
		t.Fatalf("RenderUnsupportedDiagnostics() = %q, want empty output", got)
	}
}

func TestRenderUnsupportedDiagnosticsNilWriter(t *testing.T) {
	err := RenderUnsupportedDiagnostics(nil, nil, domain.OutputFormatHuman)
	if err == nil {
		t.Fatal("RenderUnsupportedDiagnostics() error = nil, want nil writer error")
	}
	if !strings.Contains(err.Error(), "nil writer") {
		t.Fatalf("RenderUnsupportedDiagnostics() error = %v, want nil writer message", err)
	}
}

func TestRenderUnsupportedDiagnosticsHumanReturnsWriterError(t *testing.T) {
	wantErr := errors.New("write failed")
	err := RenderUnsupportedDiagnostics(failingWriter{err: wantErr}, []domain.Diagnostic{
		{Severity: domain.DiagnosticError, Code: "unsupported_artifacts", Message: "No supported artifacts"},
	}, domain.OutputFormatHuman)
	if !errors.Is(err, wantErr) {
		t.Fatalf("RenderUnsupportedDiagnostics() error = %v, want wrapping %v", err, wantErr)
	}
}

func TestRenderUnsupportedDiagnosticsJSONReturnsWriterError(t *testing.T) {
	wantErr := errors.New("write failed")
	err := RenderUnsupportedDiagnostics(failingWriter{err: wantErr}, []domain.Diagnostic{
		{Severity: domain.DiagnosticError, Code: "unsupported_artifacts", Message: "No supported artifacts"},
	}, domain.OutputFormatJSON)
	if !errors.Is(err, wantErr) {
		t.Fatalf("RenderUnsupportedDiagnostics() error = %v, want wrapping %v", err, wantErr)
	}
}

func TestRenderTasksHumanRendersDefaultColumnsInProvidedOrder(t *testing.T) {
	exit := 137
	view := domain.TasksView{
		Tasks: []domain.Task{
			{
				RowOrder: 2,
				ID:       "bb/222222",
				Status:   domain.TaskStatusFailed,
				Process:  "CALL_VARIANTS",
				Name:     "sample-2",
				Tag:      "tumor",
				Workdir:  "/runs/example/work/bb/222222",
				Exit:     &exit,
				Duration: "1h 2m",
				Realtime: "62m",
				CPUs:     "8",
				Memory:   "16 GB",
			},
			{
				RowOrder: 1,
				ID:       "aa/111111",
				Status:   domain.TaskStatusCompleted,
				Process:  "ALIGN",
				Name:     "sample-1",
				Tag:      "normal",
				Workdir:  "/runs/example/work/aa/111111",
				Duration: "5m",
				Realtime: "300s",
				CPUs:     "2",
				Memory:   "4 GB",
			},
			{ID: "cc/333333"},
		},
	}

	var buf bytes.Buffer
	if err := RenderTasksHuman(&buf, view); err != nil {
		t.Fatalf("RenderTasksHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"id\tstatus\tprocess\tname/tag\tworkdir\texit\tduration\trealtime\tcpus\tmemory",
		"bb/222222\tFAILED\tCALL_VARIANTS\tsample-2/tumor\t/runs/example/work/bb/222222\t137\t1h 2m\t62m\t8\t16 GB",
		"aa/111111\tCOMPLETED\tALIGN\tsample-1/normal\t/runs/example/work/aa/111111\t-\t5m\t300s\t2\t4 GB",
		"cc/333333\t-\t-\t-\t-\t-\t-\t-\t-\t-",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderTasksHuman() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderTasksHumanEmptyResultsClear(t *testing.T) {
	var buf bytes.Buffer
	if err := RenderTasksHuman(&buf, domain.TasksView{}); err != nil {
		t.Fatalf("RenderTasksHuman() error = %v, want nil", err)
	}

	want := "tasks: none\n"
	if got := buf.String(); got != want {
		t.Fatalf("RenderTasksHuman() = %q, want %q", got, want)
	}
}

func TestRenderTasksHumanDiagnosticsUseGitStyleBlocks(t *testing.T) {
	view := domain.TasksView{
		Tasks: []domain.Task{{ID: "aa/111111", Status: domain.TaskStatusCompleted, Process: "ALIGN", Name: "sample-1", Workdir: "/runs/example/work/aa/111111"}},
		Diagnostics: []domain.Diagnostic{
			{
				Severity: domain.DiagnosticWarning,
				Code:     "index_stale",
				Message:  "index is stale",
				Detail:   "Index path: /runs/example/.gosh/index.sqlite\nRefresh with `gosh index --refresh`.",
			},
			{
				Severity: domain.DiagnosticInfo,
				Code:     "nextflow_with_trace_recommended",
				Message:  "Run future Nextflow workflows with -with-trace",
				Detail:   "Use `nextflow run ... -with-trace` for future runs so gosh can build a complete trace-backed task index.",
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderTasksHuman(&buf, view); err != nil {
		t.Fatalf("RenderTasksHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"id\tstatus\tprocess\tname/tag\tworkdir\texit\tduration\trealtime\tcpus\tmemory",
		"aa/111111\tCOMPLETED\tALIGN\tsample-1\t/runs/example/work/aa/111111\t-\t-\t-\t-\t-",
		"",
		"warning: index is stale",
		"  code: index_stale",
		"  Index path: /runs/example/.gosh/index.sqlite",
		"hint: Refresh with `gosh index --refresh`.",
		"hint: Use `nextflow run ... -with-trace` for future runs so gosh can build a complete trace-backed task index.",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderTasksHuman() =\n%s\nwant\n%s", got, want)
	}
	if strings.Contains(buf.String(), "diagnostics:") {
		t.Fatalf("RenderTasksHuman() output = %q, did not want nested diagnostics heading", buf.String())
	}
}

func TestRenderTasksHumanErrorDiagnosticsDoNotPretendEmptyTasks(t *testing.T) {
	view := domain.TasksView{
		Diagnostics: []domain.Diagnostic{
			{
				Severity: domain.DiagnosticError,
				Code:     "tasks_unavailable_log_only",
				Message:  "gosh tasks requires a trace-backed task index; complete task/resource/status data is unavailable in log-only mode",
				Detail: strings.Join([]string{
					"Mode: log-only",
					"Run dir: /runs/log-only",
					"Selected log: /runs/log-only/.nextflow.log",
					"Only deterministic log-only failure evidence may be available; complete task rows require a Nextflow trace file.",
				}, "\n"),
			},
			{
				Severity: domain.DiagnosticInfo,
				Code:     "nextflow_with_trace_recommended",
				Message:  "Run future Nextflow workflows with -with-trace",
				Detail:   "Use `nextflow run ... -with-trace` for future runs so gosh can build a complete trace-backed task index.",
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderTasksHuman(&buf, view); err != nil {
		t.Fatalf("RenderTasksHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"error: gosh tasks requires a trace-backed task index; complete task/resource/status data is unavailable in log-only mode",
		"  code: tasks_unavailable_log_only",
		"  Mode: log-only",
		"  Run dir: /runs/log-only",
		"  Selected log: /runs/log-only/.nextflow.log",
		"  Only deterministic log-only failure evidence may be available; complete task rows require a Nextflow trace file.",
		"hint: Use `nextflow run ... -with-trace` for future runs so gosh can build a complete trace-backed task index.",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderTasksHuman() =\n%s\nwant\n%s", got, want)
	}
	for _, notWant := range []string{"tasks: none", "id\tstatus\tprocess", "diagnostics:"} {
		if strings.Contains(buf.String(), notWant) {
			t.Fatalf("RenderTasksHuman() output = %q, did not want %q", buf.String(), notWant)
		}
	}
}

func TestRenderTasksJSONPreservesProvidedOrderWithQueryMetadataAndDiagnostics(t *testing.T) {
	exit := 137
	builtAt := time.Date(2024, 4, 28, 12, 34, 56, 0, time.UTC)
	traceMod := time.Date(2024, 4, 28, 12, 0, 0, 0, time.UTC)
	logMod := time.Date(2024, 4, 28, 12, 1, 0, 0, time.UTC)
	view := domain.TasksView{
		Format: domain.OutputFormatJSON,
		Query: domain.TaskQuery{
			ProcessSubstring: "CALL",
			NameSubstring:    "sample",
			SampleSubstring:  "tumor",
			Status:           domain.TaskStatusFailed,
			StatusRaw:        "failed",
		},
		Metadata: &domain.IndexMetadata{
			SchemaVersion: 1,
			RunDir:        "/runs/example",
			IndexPath:     "/runs/example/.gosh/index.sqlite",
			Mode:          domain.IndexModeTraceBacked,
			Trace:         &domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: "/runs/example/trace.txt", ModTime: traceMod, Size: 1200},
			Log:           &domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/example/.nextflow.log", ModTime: logMod, Size: 3400},
			BuiltAt:       builtAt,
			Freshness:     domain.IndexFreshnessFresh,
			TaskCount:     3,
		},
		Tasks: []domain.Task{
			{
				RowOrder:     2,
				ID:           "bb/222222",
				Status:       domain.TaskStatusFailed,
				Process:      "CALL_VARIANTS",
				Name:         "sample-2",
				Tag:          "tumor",
				Workdir:      "/runs/example/work/bb/222222",
				Exit:         &exit,
				Duration:     "1h 2m",
				Realtime:     "62m",
				CPUs:         "8",
				Memory:       "16 GB",
				ErrorSummary: "command exited with 137",
			},
			{
				RowOrder: 1,
				ID:       "aa/111111",
				Status:   domain.TaskStatusCompleted,
				Process:  "ALIGN",
				Name:     "sample-1",
				Tag:      "normal",
				Workdir:  "/runs/example/work/aa/111111",
				Duration: "5m",
				Realtime: "300s",
				CPUs:     "2",
				Memory:   "4 GB",
			},
		},
		Diagnostics: []domain.Diagnostic{
			{Severity: domain.DiagnosticWarning, Code: "filtered", Message: "status filter applied", Detail: "2 rows matched"},
		},
	}

	var buf bytes.Buffer
	if err := RenderTasksJSON(&buf, view); err != nil {
		t.Fatalf("RenderTasksJSON() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"{",
		`  "format": "json",`,
		`  "query": {`,
		`    "process_substring": "CALL",`,
		`    "name_substring": "sample",`,
		`    "sample_substring": "tumor",`,
		`    "status": "FAILED",`,
		`    "status_raw": "failed"`,
		`  },`,
		`  "metadata": {`,
		`    "schema_version": 1,`,
		`    "run_dir": "/runs/example",`,
		`    "index_path": "/runs/example/.gosh/index.sqlite",`,
		`    "mode": "trace-backed",`,
		`    "trace": {`,
		`      "kind": "trace",`,
		`      "path": "/runs/example/trace.txt",`,
		`      "mod_time": "2024-04-28T12:00:00Z",`,
		`      "size": 1200`,
		`    },`,
		`    "log": {`,
		`      "kind": "log",`,
		`      "path": "/runs/example/.nextflow.log",`,
		`      "mod_time": "2024-04-28T12:01:00Z",`,
		`      "size": 3400`,
		`    },`,
		`    "built_at": "2024-04-28T12:34:56Z",`,
		`    "freshness": "fresh",`,
		`    "stale_reason": "",`,
		`    "task_count": 3`,
		`  },`,
		`  "tasks": [`,
		`    {`,
		`      "id": "bb/222222",`,
		`      "row_order": 2,`,
		`      "status": "FAILED",`,
		`      "process": "CALL_VARIANTS",`,
		`      "name": "sample-2",`,
		`      "tag": "tumor",`,
		`      "workdir": "/runs/example/work/bb/222222",`,
		`      "exit": 137,`,
		`      "duration": "1h 2m",`,
		`      "realtime": "62m",`,
		`      "cpus": "8",`,
		`      "memory": "16 GB",`,
		`      "error_summary": "command exited with 137"`,
		`    },`,
		`    {`,
		`      "id": "aa/111111",`,
		`      "row_order": 1,`,
		`      "status": "COMPLETED",`,
		`      "process": "ALIGN",`,
		`      "name": "sample-1",`,
		`      "tag": "normal",`,
		`      "workdir": "/runs/example/work/aa/111111",`,
		`      "exit": null,`,
		`      "duration": "5m",`,
		`      "realtime": "300s",`,
		`      "cpus": "2",`,
		`      "memory": "4 GB",`,
		`      "error_summary": ""`,
		`    }`,
		`  ],`,
		`  "diagnostics": [`,
		`    {`,
		`      "severity": "warning",`,
		`      "code": "filtered",`,
		`      "message": "status filter applied",`,
		`      "detail": "2 rows matched"`,
		`    }`,
		`  ]`,
		"}",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderTasksJSON() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderTasksJSONEmptyResultsStable(t *testing.T) {
	var buf bytes.Buffer
	if err := RenderTasksJSON(&buf, domain.TasksView{}); err != nil {
		t.Fatalf("RenderTasksJSON() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"{",
		`  "format": "json",`,
		`  "query": {`,
		`    "process_substring": "",`,
		`    "name_substring": "",`,
		`    "sample_substring": "",`,
		`    "status": "",`,
		`    "status_raw": ""`,
		`  },`,
		`  "metadata": null,`,
		`  "tasks": [],`,
		`  "diagnostics": []`,
		"}",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderTasksJSON() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderInspectJSONExactDossierStable(t *testing.T) {
	exit := 137
	view := domain.InspectView{
		Format: domain.OutputFormatJSON,
		Resolution: domain.SelectorResolution{
			Kind:     domain.SelectorResolutionExact,
			Selector: "ab/c123def",
			Diagnostics: []domain.Diagnostic{
				{Severity: domain.DiagnosticInfo, Code: "selector_exact", Message: "selector resolved exactly", Detail: "matched canonical id"},
			},
		},
		Dossier: &domain.TaskDossier{
			Task: domain.Task{
				RowOrder:     42,
				ID:           "ab/c123def",
				Status:       domain.TaskStatusFailed,
				Process:      "ALIGN_STAR",
				Name:         "tumor-sample",
				Tag:          "tumor replicate 1",
				Workdir:      "/runs/example/work/ab/c123def",
				Exit:         &exit,
				Duration:     "1h",
				Realtime:     "58m",
				CPUs:         "8",
				Memory:       "32 GB",
				ErrorSummary: "process failed: command exited with 137",
			},
			Inventory: domain.CommandFileInventory{
				Workdir: "/runs/example/work/ab/c123def",
				Files: []domain.CommandFile{
					{
						Kind:   domain.CommandFileErr,
						Path:   "/runs/example/work/ab/c123def/.command.err",
						Exists: true,
						Size:   512,
						Snippet: &domain.Snippet{
							Path:      "/runs/example/work/ab/c123def/.command.err",
							Strategy:  domain.SnippetStrategyError,
							StartLine: 10,
							EndLine:   12,
							Content:   "before\nERROR <failed>\nafter",
							MaxBytes:  4096,
						},
					},
					{
						Kind:   domain.CommandFileShell,
						Path:   "/runs/example/work/ab/c123def/.command.sh",
						Exists: true,
						Size:   42,
						Snippet: &domain.Snippet{
							Path:      "/runs/example/work/ab/c123def/.command.sh",
							Strategy:  domain.SnippetStrategyHead,
							StartLine: 1,
							EndLine:   1,
							Content:   "echo <hello>",
							Truncated: true,
							MaxBytes:  1024,
						},
					},
				},
			},
			Diagnostics: []domain.Diagnostic{
				{Severity: domain.DiagnosticWarning, Code: "command_missing", Message: ".command.log was missing"},
			},
		},
		Diagnostics: []domain.Diagnostic{
			{Severity: domain.DiagnosticWarning, Code: "inspect_partial", Message: "one command file was missing"},
		},
	}

	var buf bytes.Buffer
	if err := RenderInspectJSON(&buf, view); err != nil {
		t.Fatalf("RenderInspectJSON() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"{",
		`  "format": "json",`,
		`  "resolution": {`,
		`    "kind": "exact",`,
		`    "selector": "ab/c123def",`,
		`    "task": null,`,
		`    "matches": [],`,
		`    "diagnostics": [`,
		`      {`,
		`        "severity": "info",`,
		`        "code": "selector_exact",`,
		`        "message": "selector resolved exactly",`,
		`        "detail": "matched canonical id"`,
		`      }`,
		`    ]`,
		`  },`,
		`  "dossier": {`,
		`    "task": {`,
		`      "id": "ab/c123def",`,
		`      "row_order": 42,`,
		`      "status": "FAILED",`,
		`      "process": "ALIGN_STAR",`,
		`      "name": "tumor-sample",`,
		`      "tag": "tumor replicate 1",`,
		`      "workdir": "/runs/example/work/ab/c123def",`,
		`      "exit": 137,`,
		`      "duration": "1h",`,
		`      "realtime": "58m",`,
		`      "cpus": "8",`,
		`      "memory": "32 GB",`,
		`      "error_summary": "process failed: command exited with 137"`,
		`    },`,
		`    "inventory": {`,
		`      "workdir": "/runs/example/work/ab/c123def",`,
		`      "files": [`,
		`        {`,
		`          "kind": ".command.sh",`,
		`          "path": "/runs/example/work/ab/c123def/.command.sh",`,
		`          "exists": true,`,
		`          "size": 42,`,
		`          "snippet": {`,
		`            "path": "/runs/example/work/ab/c123def/.command.sh",`,
		`            "strategy": "head",`,
		`            "start_line": 1,`,
		`            "end_line": 1,`,
		`            "content": "echo <hello>",`,
		`            "truncated": true,`,
		`            "max_bytes": 1024`,
		`          }`,
		`        },`,
		`        {`,
		`          "kind": ".command.err",`,
		`          "path": "/runs/example/work/ab/c123def/.command.err",`,
		`          "exists": true,`,
		`          "size": 512,`,
		`          "snippet": {`,
		`            "path": "/runs/example/work/ab/c123def/.command.err",`,
		`            "strategy": "error-focused",`,
		`            "start_line": 10,`,
		`            "end_line": 12,`,
		`            "content": "before\nERROR <failed>\nafter",`,
		`            "truncated": false,`,
		`            "max_bytes": 4096`,
		`          }`,
		`        }`,
		`      ]`,
		`    },`,
		`    "diagnostics": [`,
		`      {`,
		`        "severity": "warning",`,
		`        "code": "command_missing",`,
		`        "message": ".command.log was missing",`,
		`        "detail": ""`,
		`      }`,
		`    ]`,
		`  },`,
		`  "diagnostics": [`,
		`    {`,
		`      "severity": "warning",`,
		`      "code": "inspect_partial",`,
		`      "message": "one command file was missing",`,
		`      "detail": ""`,
		`    }`,
		`  ]`,
		"}",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderInspectJSON() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderInspectJSONAmbiguousSelectorStable(t *testing.T) {
	view := domain.InspectView{
		Resolution: domain.SelectorResolution{
			Kind:     domain.SelectorResolutionAmbiguous,
			Selector: "ALIGN",
			Matches: []domain.Task{
				{RowOrder: 2, ID: "bb/222222", Status: domain.TaskStatusFailed, Process: "ALIGN_STAR", Name: "sample-2", Tag: "tumor", Workdir: "/runs/example/work/bb/222222"},
				{RowOrder: 1, ID: "aa/111111", Status: domain.TaskStatusCompleted, Process: "ALIGN_STAR", Name: "sample-1", Workdir: "/runs/example/work/aa/111111"},
			},
			Diagnostics: []domain.Diagnostic{
				{Severity: domain.DiagnosticWarning, Code: "selector_ambiguous", Message: "selector matched more than one task"},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderInspectJSON(&buf, view); err != nil {
		t.Fatalf("RenderInspectJSON() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"{",
		`  "format": "json",`,
		`  "resolution": {`,
		`    "kind": "ambiguous",`,
		`    "selector": "ALIGN",`,
		`    "task": null,`,
		`    "matches": [`,
		`      {`,
		`        "id": "bb/222222",`,
		`        "row_order": 2,`,
		`        "status": "FAILED",`,
		`        "process": "ALIGN_STAR",`,
		`        "name": "sample-2",`,
		`        "tag": "tumor",`,
		`        "workdir": "/runs/example/work/bb/222222",`,
		`        "exit": null,`,
		`        "duration": "",`,
		`        "realtime": "",`,
		`        "cpus": "",`,
		`        "memory": "",`,
		`        "error_summary": ""`,
		`      },`,
		`      {`,
		`        "id": "aa/111111",`,
		`        "row_order": 1,`,
		`        "status": "COMPLETED",`,
		`        "process": "ALIGN_STAR",`,
		`        "name": "sample-1",`,
		`        "tag": "",`,
		`        "workdir": "/runs/example/work/aa/111111",`,
		`        "exit": null,`,
		`        "duration": "",`,
		`        "realtime": "",`,
		`        "cpus": "",`,
		`        "memory": "",`,
		`        "error_summary": ""`,
		`      }`,
		`    ],`,
		`    "diagnostics": [`,
		`      {`,
		`        "severity": "warning",`,
		`        "code": "selector_ambiguous",`,
		`        "message": "selector matched more than one task",`,
		`        "detail": ""`,
		`      }`,
		`    ]`,
		`  },`,
		`  "dossier": null,`,
		`  "diagnostics": []`,
		"}",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderInspectJSON() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderInspectJSONNotFoundStable(t *testing.T) {
	view := domain.InspectView{
		Resolution: domain.SelectorResolution{
			Kind:     domain.SelectorResolutionNotFound,
			Selector: "missing-task",
			Diagnostics: []domain.Diagnostic{
				{
					Severity: domain.DiagnosticError,
					Code:     "selector_not_found",
					Message:  "selector did not match any indexed task",
					Detail:   "Try `gosh tasks` to list available task IDs.",
				},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderInspectJSON(&buf, view); err != nil {
		t.Fatalf("RenderInspectJSON() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"{",
		`  "format": "json",`,
		`  "resolution": {`,
		`    "kind": "not-found",`,
		`    "selector": "missing-task",`,
		`    "task": null,`,
		`    "matches": [],`,
		`    "diagnostics": [`,
		`      {`,
		`        "severity": "error",`,
		`        "code": "selector_not_found",`,
		`        "message": "selector did not match any indexed task",`,
		"        \"detail\": \"Try `gosh tasks` to list available task IDs.\"",
		`      }`,
		`    ]`,
		`  },`,
		`  "dossier": null,`,
		`  "diagnostics": []`,
		"}",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderInspectJSON() =\n%s\nwant\n%s", got, want)
	}
}

func decodeRenderedInspectJSON(t *testing.T, view domain.InspectView) map[string]any {
	t.Helper()

	var buf bytes.Buffer
	if err := RenderInspectJSON(&buf, view); err != nil {
		t.Fatalf("RenderInspectJSON() error = %v, want nil", err)
	}

	var payload map[string]any
	if err := json.Unmarshal(buf.Bytes(), &payload); err != nil {
		t.Fatalf("RenderInspectJSON() produced invalid JSON: %v\noutput:\n%s", err, buf.String())
	}
	return payload
}

func jsonObject(t *testing.T, value any, name string) map[string]any {
	t.Helper()
	object, ok := value.(map[string]any)
	if !ok {
		t.Fatalf("%s = %#v, want JSON object", name, value)
	}
	return object
}

func jsonArray(t *testing.T, value any, name string) []any {
	t.Helper()
	array, ok := value.([]any)
	if !ok {
		t.Fatalf("%s = %#v, want JSON array", name, value)
	}
	return array
}

func TestRenderInspectJSONLogOnlyExactDossierIncludesEvidenceAndPreservesTraceBackedFields(t *testing.T) {
	exit := 2
	view := domain.InspectView{
		EvidenceKind: domain.InspectEvidenceLogOnly,
		LogOnlyResolution: &domain.LogOnlySelectorResolution{
			Kind:     domain.SelectorResolutionExact,
			Selector: "bb/222222",
			Evidence: &domain.LogOnlyTaskEvidence{ID: "bb/222222", ObservedStatus: domain.TaskStatusFailed},
			Diagnostics: []domain.Diagnostic{
				{Severity: domain.DiagnosticInfo, Code: "selector_exact", Message: "selector resolved from log-only evidence"},
			},
		},
		LogOnlyDossier: &domain.LogOnlyTaskDossier{
			Evidence: domain.LogOnlyTaskEvidence{
				ID:             "bb/222222",
				Workdir:        "/runs/log-only/work/bb/222222",
				Process:        "PIPE:CALL",
				Name:           "tumor-02",
				ObservedStatus: domain.TaskStatusFailed,
				Exit:           &exit,
				ErrorSummary:   "No such file or directory",
				ErrorBlock:     "Command exit status: 2\nCommand error:\nNo such file or directory",
				Sources: []domain.LogOnlyEvidenceSource{
					{Kind: domain.LogOnlyEvidenceSourceLog, Path: "/runs/log-only/.nextflow.log", Detail: "failure block"},
					{Kind: domain.LogOnlyEvidenceSourceCommand, Path: "/runs/log-only/work/bb/222222/.command.err", Detail: "stderr snippet"},
				},
				Completeness:          domain.LogOnlyEvidencePartial,
				CommandFilesAvailable: true,
			},
			Inventory: domain.CommandFileInventory{
				Workdir: "/runs/log-only/work/bb/222222",
				Files: []domain.CommandFile{
					{Kind: domain.CommandFileErr, Path: "/runs/log-only/work/bb/222222/.command.err", Exists: true, Size: 128},
					{Kind: domain.CommandFileShell, Path: "/runs/log-only/work/bb/222222/.command.sh", Exists: true, Size: 33},
				},
			},
			Diagnostics: []domain.Diagnostic{
				{Severity: domain.DiagnosticWarning, Code: "log_only_partial", Message: "log-only inspect is partial"},
			},
		},
		Diagnostics: []domain.Diagnostic{
			{Severity: domain.DiagnosticWarning, Code: "inspect_partial", Message: "complete task rows require a trace file"},
		},
	}

	payload := decodeRenderedInspectJSON(t, view)
	if got := payload["format"]; got != "json" {
		t.Fatalf("format = %#v, want %q", got, "json")
	}
	if got := payload["evidence_kind"]; got != string(domain.InspectEvidenceLogOnly) {
		t.Fatalf("evidence_kind = %#v, want %q", got, domain.InspectEvidenceLogOnly)
	}

	traceResolution := jsonObject(t, payload["resolution"], "resolution")
	if got := traceResolution["task"]; got != nil {
		t.Fatalf("resolution.task = %#v, want nil for log-only output", got)
	}
	if matches := jsonArray(t, traceResolution["matches"], "resolution.matches"); len(matches) != 0 {
		t.Fatalf("resolution.matches length = %d, want 0", len(matches))
	}
	if diagnostics := jsonArray(t, traceResolution["diagnostics"], "resolution.diagnostics"); len(diagnostics) != 0 {
		t.Fatalf("resolution.diagnostics length = %d, want 0", len(diagnostics))
	}
	if got := payload["dossier"]; got != nil {
		t.Fatalf("dossier = %#v, want nil for log-only output", got)
	}

	logOnlyResolution := jsonObject(t, payload["log_only_resolution"], "log_only_resolution")
	if got := logOnlyResolution["kind"]; got != string(domain.SelectorResolutionExact) {
		t.Fatalf("log_only_resolution.kind = %#v, want %q", got, domain.SelectorResolutionExact)
	}
	if got := logOnlyResolution["selector"]; got != "bb/222222" {
		t.Fatalf("log_only_resolution.selector = %#v, want %q", got, "bb/222222")
	}
	resolutionEvidence := jsonObject(t, logOnlyResolution["evidence"], "log_only_resolution.evidence")
	if got := resolutionEvidence["id"]; got != "bb/222222" {
		t.Fatalf("log_only_resolution.evidence.id = %#v, want %q", got, "bb/222222")
	}
	if matches := jsonArray(t, logOnlyResolution["matches"], "log_only_resolution.matches"); len(matches) != 0 {
		t.Fatalf("log_only_resolution.matches length = %d, want 0", len(matches))
	}
	if diagnostics := jsonArray(t, logOnlyResolution["diagnostics"], "log_only_resolution.diagnostics"); len(diagnostics) != 1 {
		t.Fatalf("log_only_resolution.diagnostics length = %d, want 1", len(diagnostics))
	}

	logOnlyDossier := jsonObject(t, payload["log_only_dossier"], "log_only_dossier")
	dossierEvidence := jsonObject(t, logOnlyDossier["evidence"], "log_only_dossier.evidence")
	if got := dossierEvidence["exit"]; got != float64(2) {
		t.Fatalf("log_only_dossier.evidence.exit = %#v, want 2", got)
	}
	if got := dossierEvidence["error_block"]; got != "Command exit status: 2\nCommand error:\nNo such file or directory" {
		t.Fatalf("log_only_dossier.evidence.error_block = %#v", got)
	}
	if sources := jsonArray(t, dossierEvidence["sources"], "log_only_dossier.evidence.sources"); len(sources) != 2 {
		t.Fatalf("log_only_dossier.evidence.sources length = %d, want 2", len(sources))
	}
	if got := dossierEvidence["command_files_available"]; got != true {
		t.Fatalf("log_only_dossier.evidence.command_files_available = %#v, want true", got)
	}

	inventory := jsonObject(t, logOnlyDossier["inventory"], "log_only_dossier.inventory")
	files := jsonArray(t, inventory["files"], "log_only_dossier.inventory.files")
	if len(files) != 2 {
		t.Fatalf("log_only_dossier.inventory.files length = %d, want 2", len(files))
	}
	firstFile := jsonObject(t, files[0], "log_only_dossier.inventory.files[0]")
	if got := firstFile["kind"]; got != string(domain.CommandFileShell) {
		t.Fatalf("first command file kind = %#v, want shell sorted first", got)
	}
	if diagnostics := jsonArray(t, logOnlyDossier["diagnostics"], "log_only_dossier.diagnostics"); len(diagnostics) != 1 {
		t.Fatalf("log_only_dossier.diagnostics length = %d, want 1", len(diagnostics))
	}
	if diagnostics := jsonArray(t, payload["diagnostics"], "diagnostics"); len(diagnostics) != 1 {
		t.Fatalf("diagnostics length = %d, want 1", len(diagnostics))
	}
}

func TestRenderInspectJSONLogOnlyAmbiguousSelectorIncludesMatchesAndEmptyArrays(t *testing.T) {
	view := domain.InspectView{
		LogOnlyResolution: &domain.LogOnlySelectorResolution{
			Kind:     domain.SelectorResolutionAmbiguous,
			Selector: "ALIGN",
			Matches: []domain.LogOnlyTaskEvidence{
				{ID: "bb/222222", ObservedStatus: domain.TaskStatusFailed, Process: "ALIGN_STAR", Name: "tumor-02", Workdir: "/runs/log-only/work/bb/222222", CommandFilesAvailable: true},
				{ObservedStatus: domain.TaskStatusFailed, Process: "ALIGN_STAR", Name: "tumor-03"},
			},
		},
	}

	payload := decodeRenderedInspectJSON(t, view)
	if got := payload["evidence_kind"]; got != string(domain.InspectEvidenceLogOnly) {
		t.Fatalf("evidence_kind = %#v, want default %q", got, domain.InspectEvidenceLogOnly)
	}
	if got := payload["log_only_dossier"]; got != nil {
		t.Fatalf("log_only_dossier = %#v, want nil for ambiguous selector", got)
	}

	logOnlyResolution := jsonObject(t, payload["log_only_resolution"], "log_only_resolution")
	if got := logOnlyResolution["kind"]; got != string(domain.SelectorResolutionAmbiguous) {
		t.Fatalf("log_only_resolution.kind = %#v, want %q", got, domain.SelectorResolutionAmbiguous)
	}
	if got := logOnlyResolution["evidence"]; got != nil {
		t.Fatalf("log_only_resolution.evidence = %#v, want nil for ambiguous selector", got)
	}
	matches := jsonArray(t, logOnlyResolution["matches"], "log_only_resolution.matches")
	if len(matches) != 2 {
		t.Fatalf("log_only_resolution.matches length = %d, want 2", len(matches))
	}
	firstMatch := jsonObject(t, matches[0], "log_only_resolution.matches[0]")
	if got := firstMatch["command_files_available"]; got != true {
		t.Fatalf("first match command_files_available = %#v, want true", got)
	}
	secondMatch := jsonObject(t, matches[1], "log_only_resolution.matches[1]")
	if got := secondMatch["id"]; got != "" {
		t.Fatalf("second match id = %#v, want empty string", got)
	}
	if sources := jsonArray(t, secondMatch["sources"], "log_only_resolution.matches[1].sources"); len(sources) != 0 {
		t.Fatalf("second match sources length = %d, want 0", len(sources))
	}
	if diagnostics := jsonArray(t, logOnlyResolution["diagnostics"], "log_only_resolution.diagnostics"); len(diagnostics) != 0 {
		t.Fatalf("log_only_resolution.diagnostics length = %d, want 0", len(diagnostics))
	}
	if diagnostics := jsonArray(t, payload["diagnostics"], "diagnostics"); len(diagnostics) != 0 {
		t.Fatalf("diagnostics length = %d, want 0", len(diagnostics))
	}
}

func TestRenderLogOnlyInspectHumanExactDossierIncludesEvidenceSourcesWorkdirCommandSnippetsAndHints(t *testing.T) {
	exit := 2
	view := domain.InspectView{
		EvidenceKind: domain.InspectEvidenceLogOnly,
		LogOnlyResolution: &domain.LogOnlySelectorResolution{
			Kind:     domain.SelectorResolutionExact,
			Selector: "bb/222222",
		},
		LogOnlyDossier: &domain.LogOnlyTaskDossier{
			Evidence: domain.LogOnlyTaskEvidence{
				ID:             "bb/222222",
				Workdir:        "/runs/log-only/work/bb/222222",
				Process:        "PIPE:CALL",
				Name:           "tumor-02",
				ObservedStatus: domain.TaskStatusFailed,
				Exit:           &exit,
				ErrorSummary:   "No such file or directory",
				ErrorBlock:     "Command exit status: 2\nCommand error:\nNo such file or directory",
				Sources: []domain.LogOnlyEvidenceSource{
					{Kind: domain.LogOnlyEvidenceSourceLog, Path: "/runs/log-only/.nextflow.log", Detail: "failure block"},
					{Kind: domain.LogOnlyEvidenceSourceCommand, Path: "/runs/log-only/work/bb/222222/.command.err", Detail: "stderr snippet"},
				},
				Completeness:          domain.LogOnlyEvidencePartial,
				CommandFilesAvailable: true,
			},
			Inventory: domain.CommandFileInventory{
				Workdir: "/runs/log-only/work/bb/222222",
				Files: []domain.CommandFile{
					{
						Kind:   domain.CommandFileErr,
						Path:   "/runs/log-only/work/bb/222222/.command.err",
						Exists: true,
						Size:   128,
						Snippet: &domain.Snippet{
							Path:      "/runs/log-only/work/bb/222222/.command.err",
							Strategy:  domain.SnippetStrategyError,
							StartLine: 4,
							EndLine:   6,
							Content:   "before\nNo such file or directory\nafter",
							Truncated: true,
							MaxBytes:  4096,
						},
					},
					{
						Kind:   domain.CommandFileShell,
						Path:   "/runs/log-only/work/bb/222222/.command.sh",
						Exists: true,
						Size:   33,
						Snippet: &domain.Snippet{
							Path:      "/runs/log-only/work/bb/222222/.command.sh",
							Strategy:  domain.SnippetStrategyHead,
							StartLine: 1,
							EndLine:   2,
							Content:   "#!/usr/bin/env bash\nmissing-tool --input tumor",
							MaxBytes:  4096,
						},
					},
				},
			},
			Diagnostics: []domain.Diagnostic{
				{
					Severity: domain.DiagnosticWarning,
					Code:     "log_only_partial",
					Message:  "log-only inspect is partial",
					Detail:   "Selected log: /runs/log-only/.nextflow.log\nComplete task rows require a trace file.\nHint: Run future Nextflow workflows with -with-trace to produce complete task/resource/status data.",
				},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderLogOnlyInspectHuman(&buf, view); err != nil {
		t.Fatalf("RenderLogOnlyInspectHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"selector: bb/222222",
		"resolution: exact",
		"evidence_kind: log-only-partial",
		"evidence:",
		"  id: bb/222222",
		"  observed_status: FAILED",
		"  process: PIPE:CALL",
		"  name: tumor-02",
		"  workdir: /runs/log-only/work/bb/222222",
		"  exit: 2",
		"  completeness: partial",
		"  command_files_available: true",
		"  error_summary: No such file or directory",
		"  error_block:",
		"    Command exit status: 2",
		"    Command error:",
		"    No such file or directory",
		"sources:",
		"  - kind=log path=/runs/log-only/.nextflow.log detail=failure block",
		"  - kind=command-file path=/runs/log-only/work/bb/222222/.command.err detail=stderr snippet",
		"workdir:",
		"  path: /runs/log-only/work/bb/222222",
		"  available: true",
		"command_files:",
		"  - kind=.command.sh path=/runs/log-only/work/bb/222222/.command.sh exists=true size=33",
		"    snippet: strategy=head lines=1-2 truncated=false max_bytes=4096",
		"    content:",
		"      #!/usr/bin/env bash",
		"      missing-tool --input tumor",
		"  - kind=.command.err path=/runs/log-only/work/bb/222222/.command.err exists=true size=128",
		"    snippet: strategy=error-focused lines=4-6 truncated=true max_bytes=4096",
		"    content:",
		"      before",
		"      No such file or directory",
		"      after",
		"",
		"warning: log-only inspect is partial",
		"  code: log_only_partial",
		"  Selected log: /runs/log-only/.nextflow.log",
		"  Complete task rows require a trace file.",
		"hint: Run future Nextflow workflows with -with-trace to produce complete task/resource/status data.",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderLogOnlyInspectHuman() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderLogOnlyInspectHumanExactDossierWithoutWorkdirShowsUnavailableCommandFilesAndHints(t *testing.T) {
	view := domain.InspectView{
		EvidenceKind: domain.InspectEvidenceLogOnly,
		LogOnlyResolution: &domain.LogOnlySelectorResolution{
			Kind:     domain.SelectorResolutionExact,
			Selector: "NO_WORKDIR",
		},
		LogOnlyDossier: &domain.LogOnlyTaskDossier{
			Evidence: domain.LogOnlyTaskEvidence{
				Process:        "NFCORE_RNA:NO_WORKDIR",
				Name:           "sample-with-log-only-error",
				ObservedStatus: domain.TaskStatusFailed,
				ErrorSummary:   "process failed before workdir was observed",
				ErrorBlock:     "No workdir line was parseable in the selected log",
				Sources: []domain.LogOnlyEvidenceSource{
					{Kind: domain.LogOnlyEvidenceSourceLog, Path: "/runs/log-only/.nextflow.log", Detail: "failure block"},
				},
				Completeness: domain.LogOnlyEvidencePartial,
			},
			Diagnostics: []domain.Diagnostic{
				{
					Severity: domain.DiagnosticWarning,
					Code:     "inspect_workdir_unknown",
					Message:  "command-file inventory unavailable",
					Detail:   "The selected log-only evidence did not include a resolvable workdir.\nHint: Use the selected log error block or rerun with -with-trace.",
				},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderLogOnlyInspectHuman(&buf, view); err != nil {
		t.Fatalf("RenderLogOnlyInspectHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"selector: NO_WORKDIR",
		"resolution: exact",
		"evidence_kind: log-only-partial",
		"evidence:",
		"  id: -",
		"  observed_status: FAILED",
		"  process: NFCORE_RNA:NO_WORKDIR",
		"  name: sample-with-log-only-error",
		"  workdir: -",
		"  exit: -",
		"  completeness: partial",
		"  command_files_available: false",
		"  error_summary: process failed before workdir was observed",
		"  error_block:",
		"    No workdir line was parseable in the selected log",
		"sources:",
		"  - kind=log path=/runs/log-only/.nextflow.log detail=failure block",
		"workdir:",
		"  path: -",
		"  available: false",
		"command_files: none",
		"",
		"warning: command-file inventory unavailable",
		"  code: inspect_workdir_unknown",
		"  The selected log-only evidence did not include a resolvable workdir.",
		"hint: Use the selected log error block or rerun with -with-trace.",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderLogOnlyInspectHuman() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderLogOnlyInspectHumanAmbiguousSelectorShowsLogOnlyMatchesAndDiagnosticHint(t *testing.T) {
	view := domain.InspectView{
		EvidenceKind: domain.InspectEvidenceLogOnly,
		LogOnlyResolution: &domain.LogOnlySelectorResolution{
			Kind:     domain.SelectorResolutionAmbiguous,
			Selector: "ALIGN",
			Matches: []domain.LogOnlyTaskEvidence{
				{ID: "bb/222222", ObservedStatus: domain.TaskStatusFailed, Process: "ALIGN_STAR", Name: "tumor-02", Workdir: "/runs/log-only/work/bb/222222", CommandFilesAvailable: true},
				{ObservedStatus: domain.TaskStatusFailed, Process: "ALIGN_STAR", Name: "tumor-03"},
			},
			Diagnostics: []domain.Diagnostic{
				{
					Severity: domain.DiagnosticWarning,
					Code:     "selector_ambiguous",
					Message:  "selector matched more than one log-only evidence row",
					Detail:   "2 log-only evidence rows matched observed log-only fields; use a canonical id or full workdir path if available",
				},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderLogOnlyInspectHuman(&buf, view); err != nil {
		t.Fatalf("RenderLogOnlyInspectHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"selector: ALIGN",
		"resolution: ambiguous",
		"evidence_kind: log-only-partial",
		"matches:",
		"id\tstatus\tprocess\tname\tworkdir\tcommand_files_available",
		"bb/222222\tFAILED\tALIGN_STAR\ttumor-02\t/runs/log-only/work/bb/222222\ttrue",
		"-\tFAILED\tALIGN_STAR\ttumor-03\t-\tfalse",
		"",
		"warning: selector matched more than one log-only evidence row",
		"  code: selector_ambiguous",
		"  2 log-only evidence rows matched observed log-only fields; use a canonical id or full workdir path if available",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderLogOnlyInspectHuman() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderInspectHumanRoutesLogOnlyViewToLogOnlyHumanRendering(t *testing.T) {
	exit := 1
	view := domain.InspectView{
		EvidenceKind: domain.InspectEvidenceLogOnly,
		LogOnlyResolution: &domain.LogOnlySelectorResolution{
			Kind:     domain.SelectorResolutionExact,
			Selector: "bb/222222",
		},
		LogOnlyDossier: &domain.LogOnlyTaskDossier{
			Evidence: domain.LogOnlyTaskEvidence{
				ID:             "bb/222222",
				Workdir:        "/runs/log-only/work/bb/222222",
				Process:        "ALIGN_STAR",
				Name:           "tumor-02",
				ObservedStatus: domain.TaskStatusFailed,
				Exit:           &exit,
				ErrorSummary:   "No such file or directory",
				Completeness:   domain.LogOnlyEvidencePartial,
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderInspectHuman(&buf, view); err != nil {
		t.Fatalf("RenderInspectHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"selector: bb/222222",
		"resolution: exact",
		"evidence_kind: log-only-partial",
		"evidence:",
		"  id: bb/222222",
		"  observed_status: FAILED",
		"  process: ALIGN_STAR",
		"  name: tumor-02",
		"  workdir: /runs/log-only/work/bb/222222",
		"  exit: 1",
		"  completeness: partial",
		"  command_files_available: false",
		"  error_summary: No such file or directory",
		"  error_block:",
		"    -",
		"sources: none",
		"workdir:",
		"  path: /runs/log-only/work/bb/222222",
		"  available: true",
		"command_files: none",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderInspectHuman() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderInspectHumanExactDossierIncludesTaskFilesAndSnippets(t *testing.T) {
	exit := 137
	view := domain.InspectView{
		Resolution: domain.SelectorResolution{Kind: domain.SelectorResolutionExact, Selector: "ab/c123def"},
		Dossier: &domain.TaskDossier{
			Task: domain.Task{
				ID:           "ab/c123def",
				Status:       domain.TaskStatusFailed,
				Process:      "ALIGN_STAR",
				Name:         "tumor-sample",
				Tag:          "tumor replicate 1",
				Workdir:      "/runs/example/work/ab/c123def",
				Exit:         &exit,
				Duration:     "1h",
				Realtime:     "58m",
				CPUs:         "8",
				Memory:       "32 GB",
				ErrorSummary: "process failed: command exited with 137",
			},
			Inventory: domain.CommandFileInventory{
				Workdir: "/runs/example/work/ab/c123def",
				Files: []domain.CommandFile{
					{
						Kind:   domain.CommandFileShell,
						Path:   "/runs/example/work/ab/c123def/.command.sh",
						Exists: true,
						Size:   42,
						Snippet: &domain.Snippet{
							Path:      "/runs/example/work/ab/c123def/.command.sh",
							Strategy:  domain.SnippetStrategyHead,
							StartLine: 1,
							EndLine:   2,
							Content:   "echo hello\necho done",
							Truncated: true,
							MaxBytes:  4096,
						},
					},
					{
						Kind:   domain.CommandFileLog,
						Path:   "/runs/example/work/ab/c123def/.command.log",
						Exists: false,
					},
					{
						Kind:   domain.CommandFileErr,
						Path:   "/runs/example/work/ab/c123def/.command.err",
						Exists: true,
						Size:   512,
						Snippet: &domain.Snippet{
							Path:      "/runs/example/work/ab/c123def/.command.err",
							Strategy:  domain.SnippetStrategyError,
							StartLine: 10,
							EndLine:   12,
							Content:   "before\nERROR failed\nafter",
							MaxBytes:  4096,
						},
					},
				},
			},
			Diagnostics: []domain.Diagnostic{
				{Severity: domain.DiagnosticInfo, Code: "selector_exact", Message: "selector resolved exactly"},
			},
		},
		Diagnostics: []domain.Diagnostic{
			{Severity: domain.DiagnosticWarning, Code: "inspect_partial", Message: "one command file was missing"},
		},
	}

	var buf bytes.Buffer
	if err := RenderInspectHuman(&buf, view); err != nil {
		t.Fatalf("RenderInspectHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"selector: ab/c123def",
		"resolution: exact",
		"task:",
		"  id: ab/c123def",
		"  status: FAILED",
		"  process: ALIGN_STAR",
		"  name: tumor-sample",
		"  tag: tumor replicate 1",
		"  workdir: /runs/example/work/ab/c123def",
		"  exit: 137",
		"  duration: 1h",
		"  realtime: 58m",
		"  cpus: 8",
		"  memory: 32 GB",
		"  error_summary: process failed: command exited with 137",
		"command_files:",
		"  - kind=.command.sh path=/runs/example/work/ab/c123def/.command.sh exists=true size=42",
		"    snippet: strategy=head lines=1-2 truncated=true max_bytes=4096",
		"    content:",
		"      echo hello",
		"      echo done",
		"  - kind=.command.log path=/runs/example/work/ab/c123def/.command.log exists=false size=0",
		"  - kind=.command.err path=/runs/example/work/ab/c123def/.command.err exists=true size=512",
		"    snippet: strategy=error-focused lines=10-12 truncated=false max_bytes=4096",
		"    content:",
		"      before",
		"      ERROR failed",
		"      after",
		"",
		"warning: one command file was missing",
		"  code: inspect_partial",
		"",
		"info: selector resolved exactly",
		"  code: selector_exact",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderInspectHuman() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderInspectHumanAmbiguousSelectorShowsDisambiguationTable(t *testing.T) {
	view := domain.InspectView{
		Resolution: domain.SelectorResolution{
			Kind:     domain.SelectorResolutionAmbiguous,
			Selector: "ALIGN",
			Matches: []domain.Task{
				{ID: "bb/222222", Status: domain.TaskStatusFailed, Process: "ALIGN_STAR", Name: "sample-2", Tag: "tumor", Workdir: "/runs/example/work/bb/222222"},
				{ID: "aa/111111", Status: domain.TaskStatusCompleted, Process: "ALIGN_STAR", Name: "sample-1", Workdir: "/runs/example/work/aa/111111"},
			},
			Diagnostics: []domain.Diagnostic{
				{Severity: domain.DiagnosticWarning, Code: "selector_ambiguous", Message: "selector matched more than one task"},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderInspectHuman(&buf, view); err != nil {
		t.Fatalf("RenderInspectHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"selector: ALIGN",
		"resolution: ambiguous",
		"matches:",
		"id\tstatus\tprocess\tname/tag\tworkdir",
		"bb/222222\tFAILED\tALIGN_STAR\tsample-2/tumor\t/runs/example/work/bb/222222",
		"aa/111111\tCOMPLETED\tALIGN_STAR\tsample-1\t/runs/example/work/aa/111111",
		"",
		"warning: selector matched more than one task",
		"  code: selector_ambiguous",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderInspectHuman() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderInspectHumanNotFoundShowsDiagnostics(t *testing.T) {
	view := domain.InspectView{
		Resolution: domain.SelectorResolution{
			Kind:     domain.SelectorResolutionNotFound,
			Selector: "missing-task",
			Diagnostics: []domain.Diagnostic{
				{
					Severity: domain.DiagnosticError,
					Code:     "selector_not_found",
					Message:  "selector did not match any indexed task",
					Detail:   "Try `gosh tasks` to list available task IDs.",
				},
			},
		},
	}

	var buf bytes.Buffer
	if err := RenderInspectHuman(&buf, view); err != nil {
		t.Fatalf("RenderInspectHuman() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"selector: missing-task",
		"resolution: not-found",
		"matches: none",
		"",
		"error: selector did not match any indexed task",
		"  code: selector_not_found",
		"hint: Try `gosh tasks` to list available task IDs.",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderInspectHuman() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderIndexHumanReturnsWriterError(t *testing.T) {
	wantErr := errors.New("write failed")
	err := RenderIndexHuman(failingWriter{err: wantErr}, domain.IndexView{
		Diagnostics: domain.IndexDiagnostics{
			RunDir: domain.RunDir{Path: "/runs/example"},
			Metadata: &domain.IndexMetadata{
				Mode:      domain.IndexModeTraceBacked,
				Freshness: domain.IndexFreshnessFresh,
			},
		},
	})
	if !errors.Is(err, wantErr) {
		t.Fatalf("RenderIndexHuman() error = %v, want wrapping %v", err, wantErr)
	}
}

func TestRenderLogOnlyInspectHumanReturnsWriterError(t *testing.T) {
	wantErr := errors.New("write failed")
	err := RenderLogOnlyInspectHuman(failingWriter{err: wantErr}, domain.InspectView{
		EvidenceKind: domain.InspectEvidenceLogOnly,
		LogOnlyResolution: &domain.LogOnlySelectorResolution{
			Kind:     domain.SelectorResolutionExact,
			Selector: "aa/111111",
		},
		LogOnlyDossier: &domain.LogOnlyTaskDossier{
			Evidence: domain.LogOnlyTaskEvidence{ID: "aa/111111", ObservedStatus: domain.TaskStatusFailed},
		},
	})
	if !errors.Is(err, wantErr) {
		t.Fatalf("RenderLogOnlyInspectHuman() error = %v, want wrapping %v", err, wantErr)
	}
}

func TestRenderInspectHumanReturnsWriterError(t *testing.T) {
	wantErr := errors.New("write failed")
	err := RenderInspectHuman(failingWriter{err: wantErr}, domain.InspectView{
		Resolution: domain.SelectorResolution{Kind: domain.SelectorResolutionExact, Selector: "aa/111111"},
		Dossier:    &domain.TaskDossier{Task: domain.Task{ID: "aa/111111"}},
	})
	if !errors.Is(err, wantErr) {
		t.Fatalf("RenderInspectHuman() error = %v, want wrapping %v", err, wantErr)
	}
}

func TestRenderInspectJSONReturnsWriterError(t *testing.T) {
	wantErr := errors.New("write failed")
	err := RenderInspectJSON(failingWriter{err: wantErr}, domain.InspectView{
		Resolution: domain.SelectorResolution{Kind: domain.SelectorResolutionExact, Selector: "aa/111111"},
		Dossier:    &domain.TaskDossier{Task: domain.Task{ID: "aa/111111"}},
	})
	if !errors.Is(err, wantErr) {
		t.Fatalf("RenderInspectJSON() error = %v, want wrapping %v", err, wantErr)
	}
}

func TestRenderIndexJSONReturnsWriterError(t *testing.T) {
	wantErr := errors.New("write failed")
	err := RenderIndexJSON(failingWriter{err: wantErr}, domain.IndexView{})
	if !errors.Is(err, wantErr) {
		t.Fatalf("RenderIndexJSON() error = %v, want wrapping %v", err, wantErr)
	}
}

func TestRenderTasksJSONReturnsWriterError(t *testing.T) {
	wantErr := errors.New("write failed")
	err := RenderTasksJSON(failingWriter{err: wantErr}, domain.TasksView{Tasks: []domain.Task{{ID: "aa/111111"}}})
	if !errors.Is(err, wantErr) {
		t.Fatalf("RenderTasksJSON() error = %v, want wrapping %v", err, wantErr)
	}
}

func TestRenderTasksHumanReturnsWriterError(t *testing.T) {
	wantErr := errors.New("write failed")
	err := RenderTasksHuman(failingWriter{err: wantErr}, domain.TasksView{Tasks: []domain.Task{{ID: "aa/111111"}}})
	if !errors.Is(err, wantErr) {
		t.Fatalf("RenderTasksHuman() error = %v, want wrapping %v", err, wantErr)
	}
}

func TestRenderStatusHumanReturnsWriterError(t *testing.T) {
	wantErr := errors.New("write failed")
	err := RenderStatusHuman(failingWriter{err: wantErr}, domain.StatusView{})
	if !errors.Is(err, wantErr) {
		t.Fatalf("RenderStatusHuman() error = %v, want wrapping %v", err, wantErr)
	}
}

func TestRenderStatusJSONReturnsWriterError(t *testing.T) {
	wantErr := errors.New("write failed")
	err := RenderStatusJSON(failingWriter{err: wantErr}, domain.StatusView{})
	if !errors.Is(err, wantErr) {
		t.Fatalf("RenderStatusJSON() error = %v, want wrapping %v", err, wantErr)
	}
}
