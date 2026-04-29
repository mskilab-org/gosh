package render

import (
	"bytes"
	"errors"
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

func TestRenderStatusHumanLogOnlyIncludesEvidenceAndDiagnostics(t *testing.T) {
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
			FailedCount: 1,
			LogOnlyFailures: []domain.LogOnlyFailure{
				{
					ID:           "bb/222222",
					Process:      "PIPE:CALL",
					Name:         "tumor-02",
					Workdir:      "/runs/log-only/work/bb/222222",
					Exit:         &exit,
					ErrorSummary: "No such file or directory",
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
		"counts: unavailable (log-only mode; trace file required)",
		"failed_count: 1",
		"log_only_failures:",
		"  - id=bb/222222 process=PIPE:CALL name=tumor-02 workdir=/runs/log-only/work/bb/222222 exit=2 error=No such file or directory",
		"diagnostics:",
		"  - warning log_only_degraded: log-only status is degraded; complete task/resource/status data is unavailable",
		"    detail: Selected log: /runs/log-only/.nextflow.log",
		"    detail: Complete task counts require a trace file.",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderStatusHuman() =\n%s\nwant\n%s", got, want)
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

func TestRenderIndexHumanTraceBackedDiagnosticsDeterministic(t *testing.T) {
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
		"  searched_patterns: trace*.txt, .nextflow.log",
		"diagnostics:",
		"  - info source_selected: selected newest trace",
		"  - warning index_stale: index is stale",
		"    detail: Refresh with `gosh index --refresh`.",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderIndexHuman() =\n%s\nwant\n%s", got, want)
	}
}

func TestRenderIndexHumanUnsupportedNoSourcesShowsZeroTaskCount(t *testing.T) {
	view := domain.IndexView{
		Diagnostics: domain.IndexDiagnostics{
			RunDir: domain.RunDir{Path: "/runs/empty"},
			Artifacts: domain.ArtifactSet{
				RunDir: domain.RunDir{Path: "/runs/empty"},
				Mode:   domain.IndexModeUnsupported,
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
		"diagnostics:",
		"  - error unsupported_artifacts: No supported Nextflow trace or log artifacts found",
		"    detail: Searched trace patterns: trace*.txt",
		"    detail: Searched log patterns: .nextflow.log",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Fatalf("RenderIndexHuman() =\n%s\nwant\n%s", got, want)
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

func TestRenderUnsupportedDiagnosticsHumanShowsDiagnostics(t *testing.T) {
	diagnostics := []domain.Diagnostic{
		{
			Severity: domain.DiagnosticError,
			Code:     "unsupported_artifacts",
			Message:  "No supported Nextflow trace or log artifacts found",
			Detail:   "Searched trace patterns: trace*.txt\nSearched log patterns: .nextflow.log",
		},
		{
			Severity: domain.DiagnosticInfo,
			Code:     "nextflow_with_trace_recommended",
			Message:  "Re-run with `-with-trace` to enable trace-backed summaries.",
		},
	}

	var buf bytes.Buffer
	if err := RenderUnsupportedDiagnostics(&buf, diagnostics, domain.OutputFormatHuman); err != nil {
		t.Fatalf("RenderUnsupportedDiagnostics() error = %v, want nil", err)
	}

	want := strings.Join([]string{
		"diagnostics:",
		"  - error unsupported_artifacts: No supported Nextflow trace or log artifacts found",
		"    detail: Searched trace patterns: trace*.txt",
		"    detail: Searched log patterns: .nextflow.log",
		"  - info nextflow_with_trace_recommended: Re-run with `-with-trace` to enable trace-backed summaries.",
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

func TestRenderUnsupportedDiagnosticsHumanEmptyDiagnosticsClear(t *testing.T) {
	var buf bytes.Buffer
	if err := RenderUnsupportedDiagnostics(&buf, nil, domain.OutputFormatHuman); err != nil {
		t.Fatalf("RenderUnsupportedDiagnostics() error = %v, want nil", err)
	}

	want := "diagnostics: none\n"
	if got := buf.String(); got != want {
		t.Fatalf("RenderUnsupportedDiagnostics() = %q, want %q", got, want)
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
		"diagnostics:",
		"  - info selector_exact: selector resolved exactly",
		"  - warning inspect_partial: one command file was missing",
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
		"diagnostics:",
		"  - warning selector_ambiguous: selector matched more than one task",
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
		"diagnostics:",
		"  - error selector_not_found: selector did not match any indexed task",
		"    detail: Try `gosh tasks` to list available task IDs.",
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
