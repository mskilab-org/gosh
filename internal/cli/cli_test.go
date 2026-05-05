package cli

import (
	"context"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/mskilab-org/gosh/internal/domain"
	indexdb "github.com/mskilab-org/gosh/internal/index"
	"github.com/mskilab-org/gosh/internal/run"
)

func TestExecuteWritesGlobalHelpAndVersion(t *testing.T) {
	for _, args := range [][]string{{"--help"}, {"-h"}} {
		t.Run(strings.Join(args, " "), func(t *testing.T) {
			var stdout strings.Builder
			var stderr strings.Builder

			err := Execute(context.Background(), args, &stdout, &stderr, "v1.2.3-test")
			if err != nil {
				t.Fatalf("Execute(%v) returned error: %v", args, err)
			}

			got := stdout.String()
			for _, want := range []string{
				"Usage:",
				"gosh [global options] <command> [command options]",
				"gosh status [--run-dir DIR] [--results-dir DIR] [--json]",
				"gosh tasks [--run-dir DIR] [--results-dir DIR] [filters] [--json]",
				"gosh inspect <task> [--run-dir DIR] [--results-dir DIR] [--json]",
				"gosh index [--run-dir DIR] [--results-dir DIR] [--refresh] [--json]",
			} {
				if !strings.Contains(got, want) {
					t.Fatalf("help output = %q, want it to contain %q", got, want)
				}
			}
			for _, legacy := range []string{"gosh run", "gosh debug", "gosh help"} {
				if strings.Contains(got, legacy) {
					t.Fatalf("help output = %q, want clean v1 surface without legacy %q", got, legacy)
				}
			}
			if stderr.String() != "" {
				t.Fatalf("stderr = %q, want no error output for help", stderr.String())
			}
		})
	}

	t.Run("version", func(t *testing.T) {
		var stdout strings.Builder
		var stderr strings.Builder

		err := Execute(context.Background(), []string{"--version"}, &stdout, &stderr, "v1.2.3-test")
		if err != nil {
			t.Fatalf("Execute(--version) returned error: %v", err)
		}
		if got, want := stdout.String(), "gosh v1.2.3-test\n"; got != want {
			t.Fatalf("version output = %q, want %q", got, want)
		}
		if stderr.String() != "" {
			t.Fatalf("stderr = %q, want no error output for version", stderr.String())
		}
	})
}

func TestExecuteWritesCommandHelp(t *testing.T) {
	tests := []struct {
		name      string
		args      []string
		wantUsage string
	}{
		{name: "status long help", args: []string{"status", "--help"}, wantUsage: "gosh status [--run-dir DIR] [--results-dir DIR] [--json]"},
		{name: "tasks short help", args: []string{"tasks", "-h"}, wantUsage: "gosh tasks [--run-dir DIR] [--results-dir DIR] [filters] [--json]"},
		{name: "inspect long help", args: []string{"inspect", "--help"}, wantUsage: "gosh inspect <task> [--run-dir DIR] [--results-dir DIR] [--json]"},
		{name: "index long help", args: []string{"index", "--help"}, wantUsage: "gosh index [--run-dir DIR] [--results-dir DIR] [--refresh] [--json]"},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			var stdout strings.Builder
			var stderr strings.Builder

			err := Execute(context.Background(), test.args, &stdout, &stderr, "v1.2.3-test")
			if err != nil {
				t.Fatalf("Execute(%v) returned error: %v", test.args, err)
			}

			got := stdout.String()
			for _, want := range []string{test.wantUsage, "Options:"} {
				if !strings.Contains(got, want) {
					t.Fatalf("command help output = %q, want it to contain %q", got, want)
				}
			}
			if stderr.String() != "" {
				t.Fatalf("stderr = %q, want no error output for command help", stderr.String())
			}
		})
	}
}

func TestNewRootCommandDocumentsResultsDirWithoutRenamingCommands(t *testing.T) {
	root := NewRootCommand("v1.2.3-test")

	if got, want := root.Version, "v1.2.3-test"; got != want {
		t.Fatalf("Version = %q, want %q", got, want)
	}

	wantNames := []CommandName{CommandStatus, CommandTasks, CommandInspect, CommandIndex}
	if len(root.Commands) != len(wantNames) {
		t.Fatalf("len(Commands) = %d, want %d", len(root.Commands), len(wantNames))
	}
	for i, wantName := range wantNames {
		command := root.Commands[i]
		if command.Name != wantName {
			t.Fatalf("Commands[%d].Name = %q, want %q", i, command.Name, wantName)
		}
		for _, wantUsagePart := range []string{"--run-dir DIR", "--results-dir DIR"} {
			if !strings.Contains(command.Usage, wantUsagePart) {
				t.Fatalf("Commands[%d].Usage = %q, want it to contain %q", i, command.Usage, wantUsagePart)
			}
		}
	}
}

func TestExecuteDispatchesV1TopLevelCommands(t *testing.T) {
	tests := []struct {
		name string
		args func(runDir string, failedID string) []string
		want []string
	}{
		{
			name: "status",
			args: func(runDir string, failedID string) []string {
				return []string{"status", "--run-dir", runDir, "--json"}
			},
			want: []string{`"format": "json"`, `"mode": "trace-backed"`, `"failed_count": 1`},
		},
		{
			name: "tasks",
			args: func(runDir string, failedID string) []string {
				return []string{"tasks", "--run-dir", runDir, "--json", "--status", "FAILED"}
			},
			want: []string{`"format": "json"`, `"id": "bb/222222"`, `"process": "QUANT"`},
		},
		{
			name: "inspect",
			args: func(runDir string, failedID string) []string {
				return []string{"inspect", failedID, "--run-dir", runDir, "--json"}
			},
			want: []string{`"format": "json"`, `"kind": "exact"`, `"selector": "bb/222222"`, `.command.sh`},
		},
		{
			name: "index",
			args: func(runDir string, failedID string) []string {
				return []string{"index", "--run-dir", runDir, "--json"}
			},
			want: []string{`"format": "json"`, `"mode": "trace-backed"`},
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			runDir, failedID := writeExecuteTraceFixture(t)
			var stdout strings.Builder
			var stderr strings.Builder

			err := Execute(context.Background(), test.args(runDir, failedID), &stdout, &stderr, "v1.2.3-test")
			if err != nil {
				t.Fatalf("Execute(%s) returned error: %v", test.name, err)
			}

			got := stdout.String()
			for _, want := range test.want {
				if !strings.Contains(got, want) {
					t.Fatalf("Execute(%s) output = %q, want it to contain %q", test.name, got, want)
				}
			}
			if stderr.String() != "" {
				t.Fatalf("stderr = %q, want Execute to return errors instead of writing them", stderr.String())
			}
		})
	}
}

func TestExecuteParsesGlobalOptionsBeforeCommand(t *testing.T) {
	runDir, _ := writeExecuteTraceFixture(t)
	var stdout strings.Builder
	var stderr strings.Builder

	err := Execute(context.Background(), []string{"--run-dir", runDir, "--json", "tasks", "--status", "FAILED"}, &stdout, &stderr, "v1.2.3-test")
	if err != nil {
		t.Fatalf("Execute(global options before command) returned error: %v", err)
	}

	got := stdout.String()
	for _, want := range []string{`"format": "json"`, `"run_dir": "` + runDir + `"`, `"status_raw": "FAILED"`, `"id": "bb/222222"`} {
		if !strings.Contains(got, want) {
			t.Fatalf("Execute(global options before command) output = %q, want it to contain %q", got, want)
		}
	}
	if stderr.String() != "" {
		t.Fatalf("stderr = %q, want Execute to return errors instead of writing them", stderr.String())
	}
}

func TestExecuteDispatchesWithResultsDirInGlobalAndCommandLocalPositions(t *testing.T) {
	runDir, resultsDir, _, tracePath := writeExecutePipelineInfoTraceFixture(t)

	tests := []struct {
		name string
		args []string
		want []string
	}{
		{
			name: "global results dir before command selects pipeline-info trace",
			args: []string{"--run-dir", runDir, "--results-dir", resultsDir, "--json", "status"},
			want: []string{`"format": "json"`, `"mode": "trace-backed"`, `"failed_count": 1`, tracePath, "TRACE:DISPATCH"},
		},
		{
			name: "command-local results dir after command selects relative pipeline-info trace",
			args: []string{"tasks", "--run-dir", runDir, "--results-dir", "custom-results", "--json", "--status", "FAILED"},
			want: []string{`"format": "json"`, `"status_raw": "FAILED"`, `"id": "cc/333333"`, `"process": "TRACE:DISPATCH"`},
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			var stdout strings.Builder
			var stderr strings.Builder

			err := Execute(context.Background(), test.args, &stdout, &stderr, "v1.2.3-test")
			if err != nil {
				t.Fatalf("Execute(%v) returned error: %v", test.args, err)
			}

			got := stdout.String()
			for _, want := range test.want {
				if !strings.Contains(got, want) {
					t.Fatalf("Execute(%s) output = %q, want it to contain %q", test.name, got, want)
				}
			}
			if strings.Contains(got, "should-not-drive-dispatch") {
				t.Fatalf("Execute(%s) output = %q, did not want log-only fallback evidence", test.name, got)
			}
			if stderr.String() != "" {
				t.Fatalf("stderr = %q, want Execute to return errors instead of writing them", stderr.String())
			}
		})
	}
}

func TestExecuteTreatsResultsDirAsGlobalWhenScanningForRootHelpAndVersion(t *testing.T) {
	tests := []struct {
		name       string
		args       []string
		wantOutput string
	}{
		{
			name:       "root help after separate results dir value",
			args:       []string{"--results-dir", "custom-results", "--help"},
			wantOutput: "Usage:",
		},
		{
			name:       "version after equals results dir value",
			args:       []string{"--results-dir=custom-results", "--version"},
			wantOutput: "gosh v1.2.3-test\n",
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			var stdout strings.Builder
			var stderr strings.Builder

			err := Execute(context.Background(), test.args, &stdout, &stderr, "v1.2.3-test")
			if err != nil {
				t.Fatalf("Execute(%v) returned error: %v", test.args, err)
			}
			if got := stdout.String(); !strings.Contains(got, test.wantOutput) {
				t.Fatalf("stdout = %q, want it to contain %q", got, test.wantOutput)
			}
			if stderr.String() != "" {
				t.Fatalf("stderr = %q, want no error output", stderr.String())
			}
		})
	}
}

func TestExecuteRejectsUnknownAndLegacyCommandsWithoutWritingOutput(t *testing.T) {
	for _, command := range []string{"bogus", "run", "debug", "help"} {
		t.Run(command, func(t *testing.T) {
			var stdout strings.Builder
			var stderr strings.Builder

			err := Execute(context.Background(), []string{command}, &stdout, &stderr, "v1.2.3-test")
			if err == nil {
				t.Fatalf("Execute(%q) returned nil error", command)
			}
			want := "unknown command \"" + command + "\""
			if !strings.Contains(err.Error(), want) {
				t.Fatalf("error = %q, want it to contain %q", err.Error(), want)
			}
			if stdout.String() != "" {
				t.Fatalf("stdout = %q, want no output for unknown command", stdout.String())
			}
			if stderr.String() != "" {
				t.Fatalf("stderr = %q, want Execute to return errors instead of writing them", stderr.String())
			}
		})
	}
}

func TestParseGlobalOptionsParsesGlobalFlagsBeforeCommand(t *testing.T) {
	tests := []struct {
		name          string
		args          []string
		wantGlobal    GlobalOptions
		wantRemaining []string
	}{
		{
			name:          "defaults with no arguments",
			args:          nil,
			wantGlobal:    GlobalOptions{RunDir: ".", Format: domain.OutputFormatHuman},
			wantRemaining: nil,
		},
		{
			name:          "long run dir separate value and json shorthand",
			args:          []string{"--run-dir", "/runs/nf", "--json", "status"},
			wantGlobal:    GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatJSON},
			wantRemaining: []string{"status"},
		},
		{
			name:          "long run dir equals value and format flag",
			args:          []string{"--run-dir=relative/run", "--format", "json", "tasks", "--status", "FAILED"},
			wantGlobal:    GlobalOptions{RunDir: "relative/run", Format: domain.OutputFormatJSON},
			wantRemaining: []string{"tasks", "--status", "FAILED"},
		},
		{
			name:          "short run dir leaves default human format",
			args:          []string{"-d", "nf-run", "index"},
			wantGlobal:    GlobalOptions{RunDir: "nf-run", Format: domain.OutputFormatHuman},
			wantRemaining: []string{"index"},
		},
		{
			name:          "long results dir separate value",
			args:          []string{"--results-dir", "/runs/nf/results", "status"},
			wantGlobal:    GlobalOptions{RunDir: ".", ResultsDir: "/runs/nf/results", Format: domain.OutputFormatHuman},
			wantRemaining: []string{"status"},
		},
		{
			name:          "long results dir equals value with other globals",
			args:          []string{"--run-dir=relative/run", "--results-dir=relative/results", "--format=json", "tasks"},
			wantGlobal:    GlobalOptions{RunDir: "relative/run", ResultsDir: "relative/results", Format: domain.OutputFormatJSON},
			wantRemaining: []string{"tasks"},
		},
		{
			name:          "format equals value",
			args:          []string{"--format=json", "inspect", "ab/c123def"},
			wantGlobal:    GlobalOptions{RunDir: ".", Format: domain.OutputFormatJSON},
			wantRemaining: []string{"inspect", "ab/c123def"},
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			gotGlobal, gotRemaining, err := ParseGlobalOptions(test.args)
			if err != nil {
				t.Fatalf("ParseGlobalOptions(%v) returned error: %v", test.args, err)
			}
			if gotGlobal != test.wantGlobal {
				t.Fatalf("global options = %#v, want %#v", gotGlobal, test.wantGlobal)
			}
			if !reflect.DeepEqual(gotRemaining, test.wantRemaining) {
				t.Fatalf("remaining args = %#v, want %#v", gotRemaining, test.wantRemaining)
			}
		})
	}
}

func TestParseGlobalOptionsStopsBeforeCommandSpecificArgs(t *testing.T) {
	tests := []struct {
		name          string
		args          []string
		wantGlobal    GlobalOptions
		wantRemaining []string
	}{
		{
			name:          "first command stops parsing before command flags",
			args:          []string{"tasks", "--run-dir", "command-run", "--json"},
			wantGlobal:    GlobalOptions{RunDir: ".", Format: domain.OutputFormatHuman},
			wantRemaining: []string{"tasks", "--run-dir", "command-run", "--json"},
		},
		{
			name:          "separator stops parsing and is not returned",
			args:          []string{"--run-dir", "global-run", "--", "tasks", "--run-dir", "command-run", "--json"},
			wantGlobal:    GlobalOptions{RunDir: "global-run", Format: domain.OutputFormatHuman},
			wantRemaining: []string{"tasks", "--run-dir", "command-run", "--json"},
		},
		{
			name:          "command format flag remains command-specific",
			args:          []string{"--json", "inspect", "--format", "human", "ab/c123def"},
			wantGlobal:    GlobalOptions{RunDir: ".", Format: domain.OutputFormatJSON},
			wantRemaining: []string{"inspect", "--format", "human", "ab/c123def"},
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			gotGlobal, gotRemaining, err := ParseGlobalOptions(test.args)
			if err != nil {
				t.Fatalf("ParseGlobalOptions(%v) returned error: %v", test.args, err)
			}
			if gotGlobal != test.wantGlobal {
				t.Fatalf("global options = %#v, want %#v", gotGlobal, test.wantGlobal)
			}
			if !reflect.DeepEqual(gotRemaining, test.wantRemaining) {
				t.Fatalf("remaining args = %#v, want %#v", gotRemaining, test.wantRemaining)
			}
		})
	}
}

func TestParseStatusOptionsParsesInheritedAndCommandLocalFormat(t *testing.T) {
	tests := []struct {
		name   string
		global GlobalOptions
		args   []string
		want   StatusOptions
	}{
		{
			name:   "inherits run dir and human format with no args",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   nil,
			want:   StatusOptions{Global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman}},
		},
		{
			name:   "command local json shorthand overrides inherited human format",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   []string{"--json"},
			want:   StatusOptions{Global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatJSON}},
		},
		{
			name:   "command local format flag parses json value",
			global: GlobalOptions{RunDir: "relative/run", Format: domain.OutputFormatHuman},
			args:   []string{"--format", "json"},
			want:   StatusOptions{Global: GlobalOptions{RunDir: "relative/run", Format: domain.OutputFormatJSON}},
		},
		{
			name:   "command local format equals parses json value",
			global: GlobalOptions{RunDir: ".", Format: domain.OutputFormatHuman},
			args:   []string{"--format=json"},
			want:   StatusOptions{Global: GlobalOptions{RunDir: ".", Format: domain.OutputFormatJSON}},
		},
		{
			name:   "command local results dir overrides inherited empty value",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   []string{"--results-dir", "/runs/nf/results"},
			want:   StatusOptions{Global: GlobalOptions{RunDir: "/runs/nf", ResultsDir: "/runs/nf/results", Format: domain.OutputFormatHuman}},
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			got, err := ParseStatusOptions(test.global, test.args)
			if err != nil {
				t.Fatalf("ParseStatusOptions(%#v, %v) returned error: %v", test.global, test.args, err)
			}
			if got != test.want {
				t.Fatalf("status options = %#v, want %#v", got, test.want)
			}
		})
	}
}

func TestParseStatusOptionsRejectsInvalidStatusArgs(t *testing.T) {
	global := GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman}
	tests := []struct {
		name        string
		args        []string
		wantMessage string
	}{
		{name: "unexpected positional arg", args: []string{"extra"}, wantMessage: "unexpected status argument"},
		{name: "unknown flag", args: []string{"--verbose"}, wantMessage: "unknown status option"},
		{name: "missing format value", args: []string{"--format"}, wantMessage: "--format requires FORMAT"},
		{name: "empty format equals", args: []string{"--format="}, wantMessage: "--format requires FORMAT"},
		{name: "unsupported format value", args: []string{"--format", "xml"}, wantMessage: "unsupported output format"},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			got, err := ParseStatusOptions(global, test.args)
			if err == nil {
				t.Fatalf("ParseStatusOptions(%#v, %v) returned nil error", global, test.args)
			}
			if !strings.Contains(err.Error(), test.wantMessage) {
				t.Fatalf("error = %q, want it to contain %q", err.Error(), test.wantMessage)
			}
			if got != (StatusOptions{}) {
				t.Fatalf("status options on error = %#v, want zero value", got)
			}
		})
	}
}

func TestLoadCommandContextWithResultsDirResolvesRelativeResultsDirAndDiscoversPipelineInfoTrace(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "runs", "nf-run")
	cwd := filepath.Join(workspace, "cwd")
	pipelineInfoRoot := filepath.Join(runRoot, "custom-results", run.PipelineInfoDirName)
	for _, dir := range []string{runRoot, cwd, pipelineInfoRoot} {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			t.Fatalf("mkdir fixture directory %q: %v", dir, err)
		}
	}
	t.Chdir(cwd)

	tracePath := filepath.Join(pipelineInfoRoot, "execution_trace_2026-02-23.txt")
	if err := os.WriteFile(tracePath, []byte("task_id\n"), 0o644); err != nil {
		t.Fatalf("write pipeline-info trace fixture: %v", err)
	}

	got, err := loadCommandContextWithResultsDir(context.Background(), GlobalOptions{
		RunDir:     filepath.Join("..", "runs", "nf-run"),
		ResultsDir: "custom-results",
		Format:     domain.OutputFormatJSON,
	})
	if err != nil {
		t.Fatalf("loadCommandContextWithResultsDir(relative results dir) returned error: %v", err)
	}

	wantRunDir := filepath.Clean(runRoot)
	wantResultsDir := filepath.Join(wantRunDir, "custom-results")
	if got.RunDir.Path != wantRunDir {
		t.Fatalf("RunDir.Path = %q, want %q", got.RunDir.Path, wantRunDir)
	}
	if got.ResultsDir.Path != wantResultsDir {
		t.Fatalf("ResultsDir.Path = %q, want %q", got.ResultsDir.Path, wantResultsDir)
	}
	if got.Format != domain.OutputFormatJSON {
		t.Fatalf("Format = %q, want %q", got.Format, domain.OutputFormatJSON)
	}
	if got.Artifacts.RunDir != got.RunDir || got.Artifacts.ResultsDir != got.ResultsDir {
		t.Fatalf("Artifacts run/results dirs = %#v/%#v, want %#v/%#v", got.Artifacts.RunDir, got.Artifacts.ResultsDir, got.RunDir, got.ResultsDir)
	}
	if got.Artifacts.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("Artifacts.Mode = %q, want %q", got.Artifacts.Mode, domain.IndexModeTraceBacked)
	}
	if got.Artifacts.Trace == nil || got.Artifacts.Trace.Path != tracePath {
		t.Fatalf("Artifacts.Trace = %+v, want selected pipeline-info trace %q", got.Artifacts.Trace, tracePath)
	}
	if len(got.Artifacts.SearchLocations) != 3 {
		t.Fatalf("len(SearchLocations) = %d, want 3: %+v", len(got.Artifacts.SearchLocations), got.Artifacts.SearchLocations)
	}
	wantPipelineInfo := filepath.Join(wantResultsDir, run.PipelineInfoDirName)
	if got.Artifacts.SearchLocations[1].BaseDir != wantPipelineInfo {
		t.Fatalf("pipeline-info search base dir = %q, want %q", got.Artifacts.SearchLocations[1].BaseDir, wantPipelineInfo)
	}
	assertNoIndexCacheDir(t, wantRunDir)
}

func TestLoadCommandContextWithResultsDirDefaultsResultsDirAndFormat(t *testing.T) {
	runRoot := t.TempDir()
	pipelineInfoRoot := filepath.Join(runRoot, run.DefaultResultsDirName, run.PipelineInfoDirName)
	if err := os.MkdirAll(pipelineInfoRoot, 0o755); err != nil {
		t.Fatalf("mkdir default pipeline-info fixture: %v", err)
	}
	tracePath := filepath.Join(pipelineInfoRoot, "execution_trace_2026-02-23.csv")
	if err := os.WriteFile(tracePath, []byte("task_id\n"), 0o644); err != nil {
		t.Fatalf("write default pipeline-info trace fixture: %v", err)
	}

	got, err := loadCommandContextWithResultsDir(context.Background(), GlobalOptions{RunDir: runRoot})
	if err != nil {
		t.Fatalf("loadCommandContextWithResultsDir(defaults) returned error: %v", err)
	}

	wantResultsDir := filepath.Join(runRoot, run.DefaultResultsDirName)
	if got.ResultsDir.Path != wantResultsDir {
		t.Fatalf("ResultsDir.Path = %q, want default %q", got.ResultsDir.Path, wantResultsDir)
	}
	if got.Format != domain.OutputFormatHuman {
		t.Fatalf("Format = %q, want default %q", got.Format, domain.OutputFormatHuman)
	}
	if got.Artifacts.Trace == nil || got.Artifacts.Trace.Path != tracePath {
		t.Fatalf("Artifacts.Trace = %+v, want selected default pipeline-info trace %q", got.Artifacts.Trace, tracePath)
	}
	assertNoIndexCacheDir(t, runRoot)
}

func TestLoadCommandContextWithResultsDirRejectsUnsupportedFormat(t *testing.T) {
	runRoot := t.TempDir()

	got, err := loadCommandContextWithResultsDir(context.Background(), GlobalOptions{
		RunDir: runRoot,
		Format: domain.OutputFormat("xml"),
	})
	if err == nil {
		t.Fatalf("loadCommandContextWithResultsDir(unsupported format) returned nil error")
	}
	if !strings.Contains(err.Error(), "unsupported output format") {
		t.Fatalf("error = %q, want it to mention unsupported output format", err.Error())
	}
	if !reflect.DeepEqual(got, commandContext{}) {
		t.Fatalf("context on error = %#v, want zero value", got)
	}
	assertNoIndexCacheDir(t, runRoot)
}

func TestLoadCommandContextWithResultsDirUsesConfiguredResultsDirForPipelineInfoTrace(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "runs", "nf-run")
	resultsRoot := filepath.Join(workspace, "external-results")
	pipelineInfoRoot := filepath.Join(resultsRoot, run.PipelineInfoDirName)
	for _, dir := range []string{runRoot, pipelineInfoRoot} {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			t.Fatalf("mkdir fixture directory %q: %v", dir, err)
		}
	}

	logPath := filepath.Join(runRoot, ".nextflow.log")
	if err := os.WriteFile(logPath, []byte("ERROR ~ Error executing process > 'LOG:ONLY (should-not-drive-inspect)'\n"), 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}
	tracePath := filepath.Join(pipelineInfoRoot, "execution_trace_2026-05-05.csv")
	if err := os.WriteFile(tracePath, []byte("task_id\n"), 0o644); err != nil {
		t.Fatalf("write pipeline-info trace fixture: %v", err)
	}

	got, err := loadCommandContextWithResultsDir(context.Background(), GlobalOptions{
		RunDir:     runRoot,
		ResultsDir: resultsRoot,
		Format:     domain.OutputFormatJSON,
	})
	if err != nil {
		t.Fatalf("loadCommandContextWithResultsDir(configured results dir) returned error: %v", err)
	}

	if got.RunDir.Path != runRoot {
		t.Fatalf("RunDir.Path = %q, want %q", got.RunDir.Path, runRoot)
	}
	if got.ResultsDir.Path != resultsRoot {
		t.Fatalf("ResultsDir.Path = %q, want configured %q", got.ResultsDir.Path, resultsRoot)
	}
	if got.Format != domain.OutputFormatJSON {
		t.Fatalf("Format = %q, want %q", got.Format, domain.OutputFormatJSON)
	}
	if got.Artifacts.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("Artifacts.Mode = %q, want %q", got.Artifacts.Mode, domain.IndexModeTraceBacked)
	}
	if got.Artifacts.Trace == nil || got.Artifacts.Trace.Path != tracePath {
		t.Fatalf("Artifacts.Trace = %+v, want selected pipeline-info trace %q", got.Artifacts.Trace, tracePath)
	}
	if got.Artifacts.Log == nil || got.Artifacts.Log.Path != logPath {
		t.Fatalf("Artifacts.Log = %+v, want selected log %q", got.Artifacts.Log, logPath)
	}
	assertNoIndexCacheDir(t, runRoot)
}

func TestLoadCommandContextWithResultsDirResolvesRelativeResultsDirAgainstRunDir(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "runs", "nf-run")
	cwd := filepath.Join(workspace, "cwd")
	pipelineInfoRoot := filepath.Join(runRoot, "relative-results", run.PipelineInfoDirName)
	for _, dir := range []string{runRoot, cwd, pipelineInfoRoot} {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			t.Fatalf("mkdir fixture directory %q: %v", dir, err)
		}
	}
	t.Chdir(cwd)

	tracePath := filepath.Join(pipelineInfoRoot, "execution_trace_2026-05-05.tsv")
	if err := os.WriteFile(tracePath, []byte("task_id\n"), 0o644); err != nil {
		t.Fatalf("write pipeline-info trace fixture: %v", err)
	}

	got, err := loadCommandContextWithResultsDir(context.Background(), GlobalOptions{
		RunDir:     filepath.Join("..", "runs", "nf-run"),
		ResultsDir: "relative-results",
	})
	if err != nil {
		t.Fatalf("loadCommandContextWithResultsDir(relative results dir) returned error: %v", err)
	}

	wantRunDir := filepath.Clean(runRoot)
	wantResultsDir := filepath.Join(wantRunDir, "relative-results")
	if got.RunDir.Path != wantRunDir {
		t.Fatalf("RunDir.Path = %q, want %q", got.RunDir.Path, wantRunDir)
	}
	if got.ResultsDir.Path != wantResultsDir {
		t.Fatalf("ResultsDir.Path = %q, want relative to run dir %q", got.ResultsDir.Path, wantResultsDir)
	}
	if got.Format != domain.OutputFormatHuman {
		t.Fatalf("Format = %q, want default %q", got.Format, domain.OutputFormatHuman)
	}
	if got.Artifacts.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("Artifacts.Mode = %q, want %q", got.Artifacts.Mode, domain.IndexModeTraceBacked)
	}
	if got.Artifacts.Trace == nil || got.Artifacts.Trace.Path != tracePath {
		t.Fatalf("Artifacts.Trace = %+v, want selected pipeline-info trace %q", got.Artifacts.Trace, tracePath)
	}
	wantPipelineInfo := filepath.Join(wantResultsDir, run.PipelineInfoDirName)
	if len(got.Artifacts.SearchLocations) < 2 || got.Artifacts.SearchLocations[1].BaseDir != wantPipelineInfo {
		t.Fatalf("pipeline-info search location = %+v, want base dir %q", got.Artifacts.SearchLocations, wantPipelineInfo)
	}
	assertNoIndexCacheDir(t, wantRunDir)
}

func TestUnsupportedArtifactDiagnosticsForCommandUsesArtifactDiagnostics(t *testing.T) {
	runDir := domain.RunDir{Path: "/tmp/nf-run"}
	artifactDiagnostics := []domain.Diagnostic{
		{Severity: domain.DiagnosticError, Code: "custom_missing", Message: "custom missing artifacts", Detail: "custom detail"},
		{Severity: domain.DiagnosticInfo, Code: "custom_hint", Message: "custom hint", Detail: "custom hint detail"},
	}
	artifacts := domain.ArtifactSet{
		Diagnostics:     artifactDiagnostics,
		SearchLocations: run.BuildArtifactSearchLocations(runDir, domain.ResultsDir{Path: "/tmp/custom-results"}),
	}

	got := unsupportedArtifactDiagnosticsForCommand(runDir, artifacts)

	if !reflect.DeepEqual(got, artifactDiagnostics) {
		t.Fatalf("unsupportedArtifactDiagnosticsForCommand(custom diagnostics) = %+v, want %+v", got, artifactDiagnostics)
	}
}

func TestUnsupportedArtifactDiagnosticsForCommandDerivesFromSearchLocations(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "runs", "nf-run")
	resultsRoot := filepath.Join(workspace, "custom-results")
	runDir := domain.RunDir{Path: runRoot}
	locations := run.BuildArtifactSearchLocations(runDir, domain.ResultsDir{Path: resultsRoot})
	artifacts := domain.ArtifactSet{SearchLocations: locations}

	got := unsupportedArtifactDiagnosticsForCommand(runDir, artifacts)

	want := run.UnsupportedArtifactDiagnosticsFromSearchLocations(runDir, locations)
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("unsupportedArtifactDiagnosticsForCommand(search locations) = %+v, want %+v", got, want)
	}
	if !strings.Contains(got[0].Detail, filepath.Join(resultsRoot, run.PipelineInfoDirName)) {
		t.Fatalf("diagnostic detail = %q, want custom pipeline_info search location", got[0].Detail)
	}
}

func TestUnsupportedArtifactDiagnosticsForCommandFallsBackToDefaults(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	artifacts := domain.ArtifactSet{}

	got := unsupportedArtifactDiagnosticsForCommand(runDir, artifacts)

	want := run.UnsupportedArtifactDiagnostics(runDir)
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("unsupportedArtifactDiagnosticsForCommand(empty artifacts) = %+v, want %+v", got, want)
	}
	if !strings.Contains(got[0].Detail, filepath.Join(runDir.Path, run.DefaultResultsDirName, run.PipelineInfoDirName)) {
		t.Fatalf("diagnostic detail = %q, want default pipeline_info search location", got[0].Detail)
	}
}

func TestRunStatusTraceBackedRebuildsIndexAndSummarizesFailedLikeTasks(t *testing.T) {
	runDir := t.TempDir()
	failedWorkdir := filepath.Join(runDir, "work", "bb", "222222")
	abortedWorkdir := filepath.Join(runDir, "work", "cc", "333333")
	tracePath := filepath.Join(runDir, "trace.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AA/111111,COMPLETED,ALIGN,ALIGN (sample-ok),sample-ok," + filepath.Join(runDir, "work", "aa", "111111") + ",0,1m,60s,2,4 GB",
		"BB/222222,FAILED,QUANT,QUANT (sample-failed),sample-failed," + failedWorkdir + ",1,2m,120s,4,8 GB",
		"CC/333333,ABORTED,CALL,CALL (sample-aborted),sample-aborted," + abortedWorkdir + ",143,3m,180s,1,2 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write trace fixture: %v", err)
	}

	var output strings.Builder
	err := RunStatus(context.Background(), StatusOptions{
		Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatJSON},
	}, &output)
	if err != nil {
		t.Fatalf("RunStatus(trace-backed) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		`"format": "json"`,
		`"mode": "trace-backed"`,
		`"freshness": "fresh"`,
		`"failed_count": 2`,
		`"status": "FAILED"`,
		`"status": "ABORTED"`,
		failedWorkdir,
		abortedWorkdir,
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunStatus(trace-backed) output = %q, want it to contain %q", got, want)
		}
	}

	store, err := indexdb.OpenStore(context.Background(), domain.RunDir{Path: runDir})
	if err != nil {
		t.Fatalf("open store after RunStatus: %v", err)
	}
	defer store.Close()

	metadata, err := indexdb.ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata after RunStatus: %v", err)
	}
	if metadata == nil {
		t.Fatalf("ReadMetadata after RunStatus = nil, want metadata")
	}
	if metadata.Mode != domain.IndexModeTraceBacked || metadata.Freshness != domain.IndexFreshnessFresh || metadata.TaskCount != 3 {
		t.Fatalf("metadata mode/freshness/task_count = %q/%q/%d, want %q/%q/3", metadata.Mode, metadata.Freshness, metadata.TaskCount, domain.IndexModeTraceBacked, domain.IndexFreshnessFresh)
	}
}

func TestRunStatusPrefersCustomResultsDirPipelineInfoTraceOverSelectedLog(t *testing.T) {
	runDir := t.TempDir()
	resultsDir := filepath.Join(runDir, "custom-results")
	pipelineInfoRoot := filepath.Join(resultsDir, run.PipelineInfoDirName)
	if err := os.MkdirAll(pipelineInfoRoot, 0o755); err != nil {
		t.Fatalf("mkdir pipeline_info fixture: %v", err)
	}

	logPath := filepath.Join(runDir, ".nextflow.log")
	logContent := strings.Join([]string{
		"ERROR ~ Error executing process > 'LOG:ONLY (should-not-drive-status)'",
		"",
		"Command error:",
		"  log-only fallback should not be used when a pipeline-info trace exists",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}

	failedWorkdir := filepath.Join(runDir, "work", "dd", "444444")
	tracePath := filepath.Join(pipelineInfoRoot, "execution_trace_2026-04-29.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"DD/444444,FAILED,TRACE:CALL,TRACE:CALL (tumor-01),tumor-01," + failedWorkdir + ",1,5m,300s,2,6 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write pipeline-info trace fixture: %v", err)
	}

	var output strings.Builder
	err := RunStatus(context.Background(), StatusOptions{
		Global: GlobalOptions{RunDir: runDir, ResultsDir: "custom-results", Format: domain.OutputFormatJSON},
	}, &output)
	if err != nil {
		t.Fatalf("RunStatus(custom pipeline-info trace) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		`"mode": "trace-backed"`,
		`"freshness": "fresh"`,
		`"failed_count": 1`,
		tracePath,
		failedWorkdir,
		"TRACE:CALL",
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunStatus(custom pipeline-info trace) output = %q, want it to contain %q", got, want)
		}
	}
	for _, notWant := range []string{"log_only_degraded", "should-not-drive-status"} {
		if strings.Contains(got, notWant) {
			t.Fatalf("RunStatus(custom pipeline-info trace) output = %q, did not want log-only marker %q", got, notWant)
		}
	}
}

func TestLogOnlyNoParseableEvidenceStatusError(t *testing.T) {
	tests := []struct {
		name        string
		summary     domain.StatusSummary
		wantErr     bool
		wantErrText string
	}{
		{
			name:    "nil diagnostics succeeds",
			summary: domain.StatusSummary{RunDir: domain.RunDir{Path: "/runs/log-only"}},
		},
		{
			name: "log-only diagnostics with parseable evidence succeed",
			summary: domain.StatusSummary{
				RunDir: domain.RunDir{Path: "/runs/log-only"},
				Diagnostics: []domain.Diagnostic{
					{Severity: domain.DiagnosticWarning, Code: "log_only_degraded", Message: "log-only status is degraded"},
					{Severity: domain.DiagnosticInfo, Code: "nextflow_with_trace_recommended", Message: "run future workflows with -with-trace"},
				},
			},
		},
		{
			name: "no parseable evidence diagnostic returns command error",
			summary: domain.StatusSummary{
				RunDir: domain.RunDir{Path: "/runs/log-only"},
				Diagnostics: []domain.Diagnostic{
					{Severity: domain.DiagnosticWarning, Code: "log_only_degraded", Message: "log-only status is degraded"},
					{Severity: domain.DiagnosticError, Code: "log_only_no_parseable_evidence", Message: "No parseable task evidence found in selected Nextflow log"},
				},
			},
			wantErr:     true,
			wantErrText: "no parseable evidence",
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			err := logOnlyNoParseableEvidenceStatusError(test.summary)
			if test.wantErr {
				if err == nil {
					t.Fatalf("logOnlyNoParseableEvidenceStatusError(%#v) returned nil error", test.summary)
				}
				if !strings.Contains(err.Error(), test.wantErrText) {
					t.Fatalf("logOnlyNoParseableEvidenceStatusError(%#v) error = %q, want it to contain %q", test.summary, err.Error(), test.wantErrText)
				}
				return
			}
			if err != nil {
				t.Fatalf("logOnlyNoParseableEvidenceStatusError(%#v) error = %v, want nil", test.summary, err)
			}
		})
	}
}

func TestRunStatusLogOnlyRendersDegradedSummaryFromSelectedLog(t *testing.T) {
	runDir := t.TempDir()
	workdir := filepath.Join(runDir, "work", "ab", "c123def")
	if err := os.MkdirAll(workdir, 0o755); err != nil {
		t.Fatalf("create workdir fixture: %v", err)
	}
	logPath := filepath.Join(runDir, ".nextflow.log")
	logContent := strings.Join([]string{
		"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
		"ERROR ~ Error executing process > 'PIPE:ALIGN (sample-01)'",
		"",
		"Command exit status:",
		"  137",
		"",
		"Command error:",
		"  killed by scheduler",
		"",
		"Work dir:",
		"  AB/C123DEF",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}

	var output strings.Builder
	err := RunStatus(context.Background(), StatusOptions{
		Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
	}, &output)
	if err != nil {
		t.Fatalf("RunStatus(log-only) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		"mode: log-only",
		"counts: incomplete (log-only evidence; trace file required)",
		"observed_counts:",
		"  FAILED: 1",
		"failed_count: 1",
		"log_only_evidence:",
		"status=FAILED",
		"completeness=partial",
		"log_only_degraded",
		"PIPE:ALIGN",
		"sample-01",
		"killed by scheduler",
		logPath,
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunStatus(log-only) output = %q, want it to contain %q", got, want)
		}
	}
	assertNoIndexCacheDir(t, runDir)
}

func TestRunStatusLogOnlyWithoutParseableEvidenceRendersSummaryAndReturnsCommandError(t *testing.T) {
	runDir := t.TempDir()
	logPath := filepath.Join(runDir, ".nextflow.log")
	logContent := strings.Join([]string{
		"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
		"Apr-28 12:02:00.000 [main] ERROR nextflow.Session - Pipeline aborted before any process failure was reported",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}

	var output strings.Builder
	err := RunStatus(context.Background(), StatusOptions{
		Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
	}, &output)
	if err == nil {
		t.Fatalf("RunStatus(log-only no evidence) returned nil error")
	}
	if !strings.Contains(err.Error(), "no parseable evidence") {
		t.Fatalf("RunStatus(log-only no evidence) error = %q, want it to mention no parseable evidence", err.Error())
	}

	got := output.String()
	for _, want := range []string{
		"mode: log-only",
		"counts: unavailable",
		"failed_count: 0",
		"log_only_evidence: none",
		"log_only_no_parseable_evidence",
		"hint: Use `nextflow run ... -with-trace`",
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunStatus(log-only no evidence) output = %q, want it to contain %q", got, want)
		}
	}
	assertNoIndexCacheDir(t, runDir)
}

func TestRunStatusLogOnlyWithoutParseableEvidenceJSONReturnsCommandError(t *testing.T) {
	runDir := t.TempDir()
	logPath := filepath.Join(runDir, ".nextflow.log")
	logContent := strings.Join([]string{
		"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
		"Apr-28 12:02:00.000 [main] ERROR nextflow.Session - Pipeline aborted before any process failure was reported",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}

	var output strings.Builder
	err := RunStatus(context.Background(), StatusOptions{
		Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatJSON},
	}, &output)
	if err == nil {
		t.Fatalf("RunStatus(log-only no evidence JSON) returned nil error")
	}
	if !strings.Contains(err.Error(), "no parseable evidence") {
		t.Fatalf("RunStatus(log-only no evidence JSON) error = %q, want it to mention no parseable evidence", err.Error())
	}

	got := output.String()
	for _, want := range []string{
		`"format": "json"`,
		`"mode": "log-only"`,
		`"failed_count": 0`,
		`"log_only_evidence": []`,
		`"code": "log_only_no_parseable_evidence"`,
		`"code": "nextflow_with_trace_recommended"`,
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunStatus(log-only no evidence JSON) output = %q, want it to contain %q", got, want)
		}
	}
	for _, notWant := range []string{`"observed_status":`, `"command_files_available": true`} {
		if strings.Contains(got, notWant) {
			t.Fatalf("RunStatus(log-only no evidence JSON) output = %q, did not want fabricated task evidence marker %q", got, notWant)
		}
	}
	assertNoIndexCacheDir(t, runDir)
}

func TestRunStatusUnsupportedRendersDiagnosticsAndReturnsCommandError(t *testing.T) {
	runDir := t.TempDir()

	var output strings.Builder
	err := RunStatus(context.Background(), StatusOptions{
		Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
	}, &output)
	if err == nil {
		t.Fatalf("RunStatus(unsupported) returned nil error")
	}
	if !strings.Contains(err.Error(), "no usable artifacts") {
		t.Fatalf("RunStatus(unsupported) error = %q, want it to mention no usable artifacts", err.Error())
	}

	got := output.String()
	for _, want := range []string{"unsupported_artifacts", "No supported Nextflow trace or log artifacts found", "Use `nextflow run ... -with-trace`"} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunStatus(unsupported) output = %q, want it to contain %q", got, want)
		}
	}
	assertNoIndexCacheDir(t, runDir)
}

func TestRunTasksTraceBackedRebuildsIndexAndListsSourceOrder(t *testing.T) {
	runDir := t.TempDir()
	cachedWorkdir := filepath.Join(runDir, "work", "cc", "333333")
	completedWorkdir := filepath.Join(runDir, "work", "aa", "111111")
	failedWorkdir := filepath.Join(runDir, "work", "bb", "222222")
	tracePath := filepath.Join(runDir, "trace.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"CC/333333,CACHED,CACHE_STEP,CACHE_STEP (sample-cached),sample-cached," + cachedWorkdir + ",0,30s,30s,1,1 GB",
		"AA/111111,COMPLETED,ALIGN,ALIGN (sample-ok),sample-ok," + completedWorkdir + ",0,1m,60s,2,4 GB",
		"BB/222222,FAILED,QUANT,QUANT (sample-failed),sample-failed," + failedWorkdir + ",1,2m,120s,4,8 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write trace fixture: %v", err)
	}

	var output strings.Builder
	err := RunTasks(context.Background(), TasksOptions{
		Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
	}, &output)
	if err != nil {
		t.Fatalf("RunTasks(trace-backed) returned error: %v", err)
	}

	got := output.String()
	if !strings.HasPrefix(got, "id\tstatus\tprocess\tname/tag\tworkdir\texit\tduration\trealtime\tcpus\tmemory\n") {
		t.Fatalf("RunTasks(trace-backed) output = %q, want task table header first", got)
	}
	cachedIndex := strings.Index(got, "cc/333333\tCACHED")
	completedIndex := strings.Index(got, "aa/111111\tCOMPLETED")
	failedIndex := strings.Index(got, "bb/222222\tFAILED")
	if cachedIndex < 0 || completedIndex < 0 || failedIndex < 0 {
		t.Fatalf("RunTasks(trace-backed) output = %q, want all task rows", got)
	}
	if !(cachedIndex < completedIndex && completedIndex < failedIndex) {
		t.Fatalf("RunTasks(trace-backed) row order indexes = cached %d completed %d failed %d, want source order", cachedIndex, completedIndex, failedIndex)
	}
	for _, want := range []string{cachedWorkdir, completedWorkdir, failedWorkdir} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunTasks(trace-backed) output = %q, want it to contain workdir %q", got, want)
		}
	}

	store, err := indexdb.OpenStore(context.Background(), domain.RunDir{Path: runDir})
	if err != nil {
		t.Fatalf("open store after RunTasks: %v", err)
	}
	defer store.Close()

	metadata, err := indexdb.ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata after RunTasks: %v", err)
	}
	if metadata == nil {
		t.Fatalf("ReadMetadata after RunTasks = nil, want metadata")
	}
	if metadata.Mode != domain.IndexModeTraceBacked || metadata.Freshness != domain.IndexFreshnessFresh || metadata.TaskCount != 3 {
		t.Fatalf("metadata mode/freshness/task_count = %q/%q/%d, want %q/%q/3", metadata.Mode, metadata.Freshness, metadata.TaskCount, domain.IndexModeTraceBacked, domain.IndexFreshnessFresh)
	}
}

func TestRunTasksPrefersCustomResultsDirPipelineInfoTraceOverSelectedLog(t *testing.T) {
	runDir := t.TempDir()
	resultsDir := filepath.Join(runDir, "custom-results")
	pipelineInfoRoot := filepath.Join(resultsDir, run.PipelineInfoDirName)
	if err := os.MkdirAll(pipelineInfoRoot, 0o755); err != nil {
		t.Fatalf("mkdir pipeline_info fixture: %v", err)
	}

	logPath := filepath.Join(runDir, ".nextflow.log")
	logContent := strings.Join([]string{
		"ERROR ~ Error executing process > 'LOG:ONLY (should-not-drive-tasks)'",
		"",
		"Command error:",
		"  log-only fallback should not be used when a custom pipeline-info trace exists",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}

	traceWorkdir := filepath.Join(runDir, "work", "ee", "555555")
	tracePath := filepath.Join(pipelineInfoRoot, "execution_trace_2026-04-29.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"EE/555555,FAILED,TRACE:CALL,TRACE:CALL (tumor-01),tumor-01," + traceWorkdir + ",1,5m,300s,2,6 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write pipeline-info trace fixture: %v", err)
	}

	var output strings.Builder
	err := RunTasks(context.Background(), TasksOptions{
		Global: GlobalOptions{RunDir: runDir, ResultsDir: "custom-results", Format: domain.OutputFormatHuman},
	}, &output)
	if err != nil {
		t.Fatalf("RunTasks(custom pipeline-info trace) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		"id\tstatus\tprocess\tname/tag\tworkdir\texit\tduration\trealtime\tcpus\tmemory",
		"ee/555555\tFAILED\tTRACE:CALL",
		traceWorkdir,
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunTasks(custom pipeline-info trace) output = %q, want it to contain %q", got, want)
		}
	}
	for _, notWant := range []string{"tasks_unavailable_log_only", "task table unavailable", "LOG:ONLY", "should-not-drive-tasks"} {
		if strings.Contains(got, notWant) {
			t.Fatalf("RunTasks(custom pipeline-info trace) output = %q, did not want log-only marker %q", got, notWant)
		}
	}

	store, err := indexdb.OpenStore(context.Background(), domain.RunDir{Path: runDir})
	if err != nil {
		t.Fatalf("open store after RunTasks(custom pipeline-info trace): %v", err)
	}
	defer store.Close()

	metadata, err := indexdb.ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata after RunTasks(custom pipeline-info trace): %v", err)
	}
	if metadata == nil || metadata.Trace == nil {
		t.Fatalf("metadata after RunTasks(custom pipeline-info trace) = %+v, want trace-backed metadata", metadata)
	}
	if metadata.Mode != domain.IndexModeTraceBacked || metadata.TaskCount != 1 || metadata.Trace.Path != tracePath {
		t.Fatalf("metadata mode/task_count/trace = %q/%d/%q, want %q/1/%q", metadata.Mode, metadata.TaskCount, metadata.Trace.Path, domain.IndexModeTraceBacked, tracePath)
	}
}

func TestRunTasksAppliesV1FiltersWithoutFabricatingRows(t *testing.T) {
	runDir := t.TempDir()
	tracePath := filepath.Join(runDir, "trace.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AA/111111,COMPLETED,ALIGN,ALIGN (tumor-01),tumor-01," + filepath.Join(runDir, "work", "aa", "111111") + ",0,1m,60s,2,4 GB",
		"BB/222222,FAILED,ALIGN,ALIGN (tumor-02),tumor-02," + filepath.Join(runDir, "work", "bb", "222222") + ",1,2m,120s,4,8 GB",
		"CC/333333,FAILED,CALL,CALL (tumor-02),tumor-02," + filepath.Join(runDir, "work", "cc", "333333") + ",1,3m,180s,8,16 GB",
		"DD/444444,FAILED,ALIGN,ALIGN (normal-01),normal-01," + filepath.Join(runDir, "work", "dd", "444444") + ",1,4m,240s,2,4 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write trace fixture: %v", err)
	}

	var output strings.Builder
	err := RunTasks(context.Background(), TasksOptions{
		Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
		Query: domain.TaskQuery{
			ProcessSubstring: " align ",
			SampleSubstring:  "TUMOR-02",
			StatusRaw:        "failed",
		},
	}, &output)
	if err != nil {
		t.Fatalf("RunTasks(filtered) returned error: %v", err)
	}

	got := output.String()
	if !strings.Contains(got, "bb/222222\tFAILED\tALIGN") {
		t.Fatalf("RunTasks(filtered) output = %q, want matching failed ALIGN tumor row", got)
	}
	for _, notWant := range []string{"aa/111111", "cc/333333", "dd/444444", "tasks: none"} {
		if strings.Contains(got, notWant) {
			t.Fatalf("RunTasks(filtered) output = %q, did not want %q", got, notWant)
		}
	}
}

func TestRunTasksJSONRendersFilteredRowsAndMetadata(t *testing.T) {
	runDir := t.TempDir()
	completedWorkdir := filepath.Join(runDir, "work", "aa", "111111")
	tracePath := filepath.Join(runDir, "trace.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AA/111111,COMPLETED,ALIGN,ALIGN (sample-ok),sample-ok," + completedWorkdir + ",0,1m,60s,2,4 GB",
		"BB/222222,FAILED,QUANT,QUANT (sample-failed),sample-failed," + filepath.Join(runDir, "work", "bb", "222222") + ",1,2m,120s,4,8 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write trace fixture: %v", err)
	}

	var output strings.Builder
	err := RunTasks(context.Background(), TasksOptions{
		Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatJSON},
		Query:  domain.TaskQuery{StatusRaw: "completed"},
	}, &output)
	if err != nil {
		t.Fatalf("RunTasks(json) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		`"format": "json"`,
		`"status_raw": "completed"`,
		`"metadata": {`,
		`"mode": "trace-backed"`,
		`"freshness": "fresh"`,
		`"task_count": 2`,
		`"tasks": [`,
		`"id": "aa/111111"`,
		`"status": "COMPLETED"`,
		completedWorkdir,
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunTasks(json) output = %q, want it to contain %q", got, want)
		}
	}
	for _, notWant := range []string{`"id": "bb/222222"`, `"status": "FAILED"`} {
		if strings.Contains(got, notWant) {
			t.Fatalf("RunTasks(json) output = %q, did not want filtered row marker %q", got, notWant)
		}
	}
}

func TestRunTasksLogOnlyHumanDiagnosticsAreGitStyle(t *testing.T) {
	runDir := t.TempDir()
	logPath := filepath.Join(runDir, ".nextflow.log")
	logContent := strings.Join([]string{
		"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
		"ERROR ~ Error executing process > 'PIPE:ALIGN (sample-01)'",
		"Command exit status:",
		"  137",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}

	var output strings.Builder
	err := RunTasks(context.Background(), TasksOptions{
		Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
	}, &output)
	if err == nil {
		t.Fatalf("RunTasks(log-only human diagnostics) returned nil error")
	}
	if !strings.Contains(err.Error(), "log-only") {
		t.Fatalf("RunTasks(log-only human diagnostics) error = %q, want it to mention log-only", err.Error())
	}

	got := output.String()
	if !strings.HasPrefix(got, "error: task table unavailable without a trace file\n") {
		t.Fatalf("RunTasks(log-only human diagnostics) output = %q, want Git-style primary error first", got)
	}
	for _, want := range []string{
		"  code: tasks_unavailable_log_only",
		"  Mode: log-only",
		"  Run dir: " + runDir,
		"  Selected log: " + logPath,
		"complete task/resource/status data is unavailable",
		"hint: Use `nextflow run ... -with-trace`",
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunTasks(log-only human diagnostics) output = %q, want it to contain %q", got, want)
		}
	}
	for _, notWant := range []string{"diagnostics:", "id\tstatus\tprocess", "tasks: none", "PIPE:ALIGN"} {
		if strings.Contains(got, notWant) {
			t.Fatalf("RunTasks(log-only human diagnostics) output = %q, did not want fabricated task-table marker %q", got, notWant)
		}
	}
	assertNoIndexCacheDir(t, runDir)
}

func TestRunTasksLogOnlyAndUnsupportedStateTasksUnavailable(t *testing.T) {
	t.Run("log-only", func(t *testing.T) {
		runDir := t.TempDir()
		logPath := filepath.Join(runDir, ".nextflow.log")
		logContent := strings.Join([]string{
			"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
			"ERROR ~ Error executing process > 'PIPE:ALIGN (sample-01)'",
			"Command exit status:",
			"  137",
			"",
		}, "\n")
		if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
			t.Fatalf("write log fixture: %v", err)
		}

		var output strings.Builder
		err := RunTasks(context.Background(), TasksOptions{
			Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
		}, &output)
		if err == nil {
			t.Fatalf("RunTasks(log-only) returned nil error")
		}
		if !strings.Contains(err.Error(), "log-only") {
			t.Fatalf("RunTasks(log-only) error = %q, want it to mention log-only", err.Error())
		}

		got := output.String()
		for _, want := range []string{"tasks_unavailable_log_only", "complete task/resource/status data is unavailable", "Use `nextflow run ... -with-trace`", logPath} {
			if !strings.Contains(got, want) {
				t.Fatalf("RunTasks(log-only) output = %q, want it to contain %q", got, want)
			}
		}
		for _, notWant := range []string{"id\tstatus\tprocess", "tasks: none", "PIPE:ALIGN"} {
			if strings.Contains(got, notWant) {
				t.Fatalf("RunTasks(log-only) output = %q, did not want fabricated task-table marker %q", got, notWant)
			}
		}
		assertNoIndexCacheDir(t, runDir)
	})

	t.Run("unsupported", func(t *testing.T) {
		runDir := t.TempDir()

		var output strings.Builder
		err := RunTasks(context.Background(), TasksOptions{
			Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatJSON},
		}, &output)
		if err == nil {
			t.Fatalf("RunTasks(unsupported) returned nil error")
		}
		if !strings.Contains(err.Error(), "no usable artifacts") {
			t.Fatalf("RunTasks(unsupported) error = %q, want it to mention no usable artifacts", err.Error())
		}

		got := output.String()
		for _, want := range []string{`"format": "json"`, "tasks_unavailable_unsupported", "unsupported_artifacts", "complete task/resource/status data is unavailable", "Run future Nextflow workflows with -with-trace"} {
			if !strings.Contains(got, want) {
				t.Fatalf("RunTasks(unsupported) output = %q, want it to contain %q", got, want)
			}
		}
		if strings.Contains(got, `"tasks"`) {
			t.Fatalf("RunTasks(unsupported) output = %q, did not want fabricated tasks array", got)
		}
		assertNoIndexCacheDir(t, runDir)
	})
}

func TestRunInspectTraceBackedExactInventoriesCommandFilesWithoutExecuting(t *testing.T) {
	runDir := t.TempDir()
	workdir := filepath.Join(runDir, "work", "ab", "c123def")
	if err := os.MkdirAll(workdir, 0o755); err != nil {
		t.Fatalf("create workdir fixture: %v", err)
	}

	tracePath := filepath.Join(runDir, "trace.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AB/C123DEF,FAILED,ALIGN,ALIGN (sample-01),sample-01," + workdir + ",137,2m,120s,4,8 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write trace fixture: %v", err)
	}

	sentinel := filepath.Join(runDir, "command-file-was-executed")
	commandFiles := map[string]string{
		".command.sh":  "echo inspect\ntouch " + sentinel + "\n",
		".command.err": "before\nERROR failed deterministically\nafter\n",
		".command.run": "#!/usr/bin/env bash\ntouch " + sentinel + "\n",
	}
	for name, content := range commandFiles {
		if err := os.WriteFile(filepath.Join(workdir, name), []byte(content), 0o755); err != nil {
			t.Fatalf("write %s fixture: %v", name, err)
		}
	}

	var output strings.Builder
	err := RunInspect(context.Background(), InspectOptions{
		Global:   GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
		Selector: "ab/c123def",
	}, &output)
	if err != nil {
		t.Fatalf("RunInspect(trace-backed exact) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		"selector: ab/c123def",
		"resolution: exact",
		"task:",
		"id: ab/c123def",
		"status: FAILED",
		"process: ALIGN",
		workdir,
		"command_files:",
		"kind=.command.sh path=" + filepath.Join(workdir, ".command.sh") + " exists=true",
		"echo inspect",
		"kind=.command.log path=" + filepath.Join(workdir, ".command.log") + " exists=false",
		"kind=.command.err path=" + filepath.Join(workdir, ".command.err") + " exists=true",
		"ERROR failed deterministically",
		"kind=.command.run path=" + filepath.Join(workdir, ".command.run") + " exists=true",
		"selector_exact",
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunInspect(trace-backed exact) output = %q, want it to contain %q", got, want)
		}
	}
	if _, err := os.Stat(sentinel); !os.IsNotExist(err) {
		t.Fatalf("command file sentinel stat error = %v, want command files to be inventoried but not executed", err)
	}
}

func TestRunInspectPrefersCustomResultsDirPipelineInfoTraceOverSelectedLog(t *testing.T) {
	runDir := t.TempDir()
	resultsDir := filepath.Join(runDir, "custom-results")
	pipelineInfoRoot := filepath.Join(resultsDir, run.PipelineInfoDirName)
	traceWorkdir := filepath.Join(runDir, "work", "ff", "666666")
	for _, dir := range []string{pipelineInfoRoot, traceWorkdir} {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			t.Fatalf("create inspect fixture dir %s: %v", dir, err)
		}
	}

	if err := os.WriteFile(filepath.Join(traceWorkdir, ".command.sh"), []byte("#!/usr/bin/env bash\necho custom inspect\n"), 0o644); err != nil {
		t.Fatalf("write command fixture: %v", err)
	}

	logPath := filepath.Join(runDir, ".nextflow.log")
	logContent := strings.Join([]string{
		"ERROR ~ Error executing process > 'LOG:ONLY (should-not-drive-inspect)'",
		"",
		"Command error:",
		"  log-only fallback should not be used when a custom pipeline-info trace exists",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}

	tracePath := filepath.Join(pipelineInfoRoot, "execution_trace_inspect.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"FF/666666,FAILED,TRACE:INSPECT,TRACE:INSPECT (tumor-01),tumor-01," + traceWorkdir + ",1,5m,300s,2,6 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write pipeline-info trace fixture: %v", err)
	}

	var output strings.Builder
	err := RunInspect(context.Background(), InspectOptions{
		Global:   GlobalOptions{RunDir: runDir, ResultsDir: "custom-results", Format: domain.OutputFormatHuman},
		Selector: "ff/666666",
	}, &output)
	if err != nil {
		t.Fatalf("RunInspect(custom pipeline-info trace) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		"selector: ff/666666",
		"resolution: exact",
		"task:",
		"id: ff/666666",
		"status: FAILED",
		"process: TRACE:INSPECT",
		traceWorkdir,
		"command_files:",
		"kind=.command.sh path=" + filepath.Join(traceWorkdir, ".command.sh") + " exists=true",
		"echo custom inspect",
		"selector_exact",
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunInspect(custom pipeline-info trace) output = %q, want it to contain %q", got, want)
		}
	}
	for _, notWant := range []string{"log_only_partial", "log-only inspect is partial", "LOG:ONLY", "should-not-drive-inspect"} {
		if strings.Contains(got, notWant) {
			t.Fatalf("RunInspect(custom pipeline-info trace) output = %q, did not want log-only marker %q", got, notWant)
		}
	}
}

func TestRunInspectTraceBackedAmbiguousAndNotFoundRenderResolutionWithoutDossier(t *testing.T) {
	runDir := t.TempDir()
	firstWorkdir := filepath.Join(runDir, "work", "aa", "111111")
	secondWorkdir := filepath.Join(runDir, "work", "bb", "222222")
	tracePath := filepath.Join(runDir, "trace.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AA/111111,COMPLETED,ALIGN,ALIGN (sample-01),sample-01," + firstWorkdir + ",0,1m,60s,2,4 GB",
		"BB/222222,FAILED,ALIGN,ALIGN (sample-02),sample-02," + secondWorkdir + ",1,2m,120s,4,8 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write trace fixture: %v", err)
	}

	t.Run("ambiguous", func(t *testing.T) {
		var output strings.Builder
		err := RunInspect(context.Background(), InspectOptions{
			Global:   GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
			Selector: "ALIGN",
		}, &output)
		if err != nil {
			t.Fatalf("RunInspect(ambiguous) returned error: %v", err)
		}

		got := output.String()
		for _, want := range []string{"resolution: ambiguous", "matches:", "aa/111111", "bb/222222", "selector_ambiguous"} {
			if !strings.Contains(got, want) {
				t.Fatalf("RunInspect(ambiguous) output = %q, want it to contain %q", got, want)
			}
		}
		for _, notWant := range []string{"task:\n", "command_files:"} {
			if strings.Contains(got, notWant) {
				t.Fatalf("RunInspect(ambiguous) output = %q, did not want dossier marker %q", got, notWant)
			}
		}
	})

	t.Run("not found", func(t *testing.T) {
		var output strings.Builder
		err := RunInspect(context.Background(), InspectOptions{
			Global:   GlobalOptions{RunDir: runDir, Format: domain.OutputFormatJSON},
			Selector: "missing-task",
		}, &output)
		if err != nil {
			t.Fatalf("RunInspect(not found) returned error: %v", err)
		}

		got := output.String()
		for _, want := range []string{`"format": "json"`, `"kind": "not-found"`, `"selector": "missing-task"`, `"dossier": null`, "selector_not_found"} {
			if !strings.Contains(got, want) {
				t.Fatalf("RunInspect(not found) output = %q, want it to contain %q", got, want)
			}
		}
		for _, notWant := range []string{".command.sh", `"inventory"`} {
			if strings.Contains(got, notWant) {
				t.Fatalf("RunInspect(not found) output = %q, did not want inventory marker %q", got, notWant)
			}
		}
	})
}

func TestRunInspectTraceBackedExactUnknownWorkdirRendersDiagnosticWithoutInventory(t *testing.T) {
	runDir := t.TempDir()
	tracePath := filepath.Join(runDir, "trace.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AB/C123DEF,FAILED,ALIGN,ALIGN (sample-01),sample-01,,137,2m,120s,4,8 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write trace fixture: %v", err)
	}

	var output strings.Builder
	err := RunInspect(context.Background(), InspectOptions{
		Global:   GlobalOptions{RunDir: runDir, Format: domain.OutputFormatJSON},
		Selector: "ab/c123def",
	}, &output)
	if err != nil {
		t.Fatalf("RunInspect(exact unknown workdir) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		`"kind": "exact"`,
		`"id": "ab/c123def"`,
		`"workdir": ""`,
		`"inventory": {`,
		`"files": []`,
		"inspect_workdir_unknown",
		"command-file inventory skipped",
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunInspect(exact unknown workdir) output = %q, want it to contain %q", got, want)
		}
	}
	if strings.Contains(got, ".command.sh") {
		t.Fatalf("RunInspect(exact unknown workdir) output = %q, did not want command file paths for unknown workdir", got)
	}
}

func TestRunInspectLogOnlyExactRendersPartialDossierFromEvidence(t *testing.T) {
	runDir := t.TempDir()
	workdir := filepath.Join(runDir, "work", "ab", "c123def")
	if err := os.MkdirAll(workdir, 0o755); err != nil {
		t.Fatalf("create workdir fixture: %v", err)
	}

	logPath := filepath.Join(runDir, ".nextflow.log")
	logContent := strings.Join([]string{
		"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
		"ERROR ~ Error executing process > 'PIPE:ALIGN (sample-01)'",
		"",
		"Command exit status:",
		"  137",
		"",
		"Command error:",
		"  killed by scheduler",
		"",
		"Work dir:",
		"  " + workdir,
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}

	sentinel := filepath.Join(runDir, "log-only-command-file-was-executed")
	commandFiles := map[string]string{
		".command.sh":  "echo log-only inspect\ntouch " + sentinel + "\n",
		".command.err": "before\nstderr from failed task\nafter\n",
		".command.run": "#!/usr/bin/env bash\ntouch " + sentinel + "\n",
	}
	for name, content := range commandFiles {
		if err := os.WriteFile(filepath.Join(workdir, name), []byte(content), 0o755); err != nil {
			t.Fatalf("write %s fixture: %v", name, err)
		}
	}

	var output strings.Builder
	err := RunInspect(context.Background(), InspectOptions{
		Global:   GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
		Selector: "ab/c123def",
	}, &output)
	if err != nil {
		t.Fatalf("RunInspect(log-only exact) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		"selector: ab/c123def",
		"resolution: exact",
		"evidence_kind: log-only-partial",
		"evidence:",
		"id: ab/c123def",
		"observed_status: FAILED",
		"process: PIPE:ALIGN",
		"name: sample-01",
		"workdir: " + workdir,
		"exit: 137",
		"command_files_available: true",
		"killed by scheduler",
		"sources:",
		logPath,
		"workdir:",
		"available: true",
		"command_files:",
		"kind=.command.sh path=" + filepath.Join(workdir, ".command.sh") + " exists=true",
		"echo log-only inspect",
		"kind=.command.err path=" + filepath.Join(workdir, ".command.err") + " exists=true",
		"stderr from failed task",
		"kind=.command.run path=" + filepath.Join(workdir, ".command.run") + " exists=true",
		"log_only_partial",
		"Run future Nextflow workflows with -with-trace",
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunInspect(log-only exact) output = %q, want it to contain %q", got, want)
		}
	}
	if _, err := os.Stat(sentinel); !os.IsNotExist(err) {
		t.Fatalf("command file sentinel stat error = %v, want command files to be inventoried but not executed", err)
	}
	assertNoIndexCacheDir(t, runDir)
}

func TestRunInspectLogOnlyExactJSONResolvesFullWorkdirSelector(t *testing.T) {
	runDir := t.TempDir()
	workdir := filepath.Join(runDir, "work", "bc", "234567")
	if err := os.MkdirAll(workdir, 0o755); err != nil {
		t.Fatalf("create workdir fixture: %v", err)
	}

	logPath := filepath.Join(runDir, ".nextflow.log")
	logContent := strings.Join([]string{
		"ERROR ~ Error executing process > 'PIPE:CALL (tumor-02)'",
		"",
		"Command exit status:",
		"  2",
		"",
		"Command error:",
		"  missing input file",
		"",
		"Work dir:",
		"  " + workdir,
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}
	if err := os.WriteFile(filepath.Join(workdir, ".command.sh"), []byte("missing-tool --input tumor\n"), 0o755); err != nil {
		t.Fatalf("write .command.sh fixture: %v", err)
	}

	var output strings.Builder
	err := RunInspect(context.Background(), InspectOptions{
		Global:   GlobalOptions{RunDir: runDir, Format: domain.OutputFormatJSON},
		Selector: workdir,
	}, &output)
	if err != nil {
		t.Fatalf("RunInspect(log-only full workdir JSON) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		`"format": "json"`,
		`"evidence_kind": "log-only-partial"`,
		`"dossier": null`,
		`"log_only_resolution": {`,
		`"kind": "exact"`,
		`"selector": "` + workdir + `"`,
		`"log_only_dossier": {`,
		`"id": "bc/234567"`,
		`"workdir": "` + workdir + `"`,
		`"observed_status": "FAILED"`,
		`"exit": 2`,
		`"inventory": {`,
		`.command.sh`,
		`"task": null`,
		"selector_exact",
		"log_only_partial",
		logPath,
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunInspect(log-only full workdir JSON) output = %q, want it to contain %q", got, want)
		}
	}
	assertNoIndexCacheDir(t, runDir)
}

func TestRunInspectLogOnlyAmbiguousAndNoWorkdirEvidenceRenderWithoutGuessing(t *testing.T) {
	t.Run("ambiguous", func(t *testing.T) {
		runDir := t.TempDir()
		firstWorkdir := filepath.Join(runDir, "work", "aa", "111111")
		secondWorkdir := filepath.Join(runDir, "work", "bb", "222222")
		logPath := filepath.Join(runDir, ".nextflow.log")
		logContent := strings.Join([]string{
			"ERROR ~ Error executing process > 'ALIGN_STAR (tumor-01)'",
			"Work dir:",
			"  " + firstWorkdir,
			"",
			"ERROR ~ Error executing process > 'ALIGN_BWA (tumor-02)'",
			"Work dir:",
			"  " + secondWorkdir,
			"",
		}, "\n")
		if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
			t.Fatalf("write log fixture: %v", err)
		}

		var output strings.Builder
		err := RunInspect(context.Background(), InspectOptions{
			Global:   GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
			Selector: "ALIGN",
		}, &output)
		if err != nil {
			t.Fatalf("RunInspect(log-only ambiguous) returned error: %v", err)
		}

		got := output.String()
		for _, want := range []string{"resolution: ambiguous", "evidence_kind: log-only-partial", "matches:", "aa/111111", "bb/222222", firstWorkdir, secondWorkdir, "selector_ambiguous"} {
			if !strings.Contains(got, want) {
				t.Fatalf("RunInspect(log-only ambiguous) output = %q, want it to contain %q", got, want)
			}
		}
		for _, notWant := range []string{"evidence:\n", "command_files:"} {
			if strings.Contains(got, notWant) {
				t.Fatalf("RunInspect(log-only ambiguous) output = %q, did not want guessed dossier marker %q", got, notWant)
			}
		}
		assertNoIndexCacheDir(t, runDir)
	})

	t.Run("exact without workdir", func(t *testing.T) {
		runDir := t.TempDir()
		logPath := filepath.Join(runDir, ".nextflow.log")
		logContent := strings.Join([]string{
			"ERROR ~ Error executing process > 'NO_WORKDIR (sample-with-log-only-error)'",
			"",
			"Command error:",
			"  process failed before workdir was observed",
			"",
		}, "\n")
		if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
			t.Fatalf("write log fixture: %v", err)
		}

		var output strings.Builder
		err := RunInspect(context.Background(), InspectOptions{
			Global:   GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
			Selector: "NO_WORKDIR",
		}, &output)
		if err != nil {
			t.Fatalf("RunInspect(log-only no workdir) returned error: %v", err)
		}

		got := output.String()
		for _, want := range []string{
			"resolution: exact",
			"evidence_kind: log-only-partial",
			"process: NO_WORKDIR",
			"name: sample-with-log-only-error",
			"workdir: -",
			"available: false",
			"command_files: none",
			"process failed before workdir was observed",
			"inspect_workdir_unknown",
			"command-file inventory unavailable",
			"Use the selected log error block or rerun with -with-trace",
			logPath,
		} {
			if !strings.Contains(got, want) {
				t.Fatalf("RunInspect(log-only no workdir) output = %q, want it to contain %q", got, want)
			}
		}
		assertNoIndexCacheDir(t, runDir)
	})
}

func TestRunInspectLogOnlyWithoutParseableEvidenceRendersNotFoundDiagnostic(t *testing.T) {
	runDir := t.TempDir()
	logPath := filepath.Join(runDir, ".nextflow.log")
	logContent := strings.Join([]string{
		"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
		"Apr-28 12:02:00.000 [main] ERROR nextflow.Session - Pipeline aborted before any process failure was reported",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(logContent), 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}

	var output strings.Builder
	err := RunInspect(context.Background(), InspectOptions{
		Global:   GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
		Selector: "missing-task",
	}, &output)
	if err != nil {
		t.Fatalf("RunInspect(log-only no parseable evidence) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		"selector: missing-task",
		"resolution: not-found",
		"evidence_kind: log-only-partial",
		"matches: none",
		"selector_not_found",
		"log_only_no_parseable_evidence",
		"No parseable task evidence found in selected Nextflow log",
		logPath,
		"Run future Nextflow workflows with -with-trace",
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunInspect(log-only no parseable evidence) output = %q, want it to contain %q", got, want)
		}
	}
	if strings.Contains(got, "command_files:") {
		t.Fatalf("RunInspect(log-only no parseable evidence) output = %q, did not want command-file inventory", got)
	}
	assertNoIndexCacheDir(t, runDir)
}

func TestRunInspectUnsupportedDoesNotFabricateTasks(t *testing.T) {
	runDir := t.TempDir()

	var output strings.Builder
	err := RunInspect(context.Background(), InspectOptions{
		Global:   GlobalOptions{RunDir: runDir, Format: domain.OutputFormatJSON},
		Selector: "ab/c123def",
	}, &output)
	if err == nil {
		t.Fatalf("RunInspect(unsupported) returned nil error")
	}
	if !strings.Contains(err.Error(), "no usable artifacts") {
		t.Fatalf("RunInspect(unsupported) error = %q, want it to mention no usable artifacts", err.Error())
	}

	got := output.String()
	for _, want := range []string{`"format": "json"`, "inspect_unavailable_unsupported", "unsupported_artifacts", "Run future Nextflow workflows with -with-trace"} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunInspect(unsupported) output = %q, want it to contain %q", got, want)
		}
	}
	for _, notWant := range []string{`"resolution"`, `"dossier"`, `"task"`} {
		if strings.Contains(got, notWant) {
			t.Fatalf("RunInspect(unsupported) output = %q, did not want fabricated inspect JSON marker %q", got, notWant)
		}
	}
	assertNoIndexCacheDir(t, runDir)
}

func TestParseTasksOptionsParsesFiltersAndCommandLocalOptions(t *testing.T) {
	tests := []struct {
		name   string
		global GlobalOptions
		args   []string
		want   TasksOptions
	}{
		{
			name:   "inherits run dir and human format with no filters",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   nil,
			want:   TasksOptions{Global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman}},
		},
		{
			name:   "parses separate process name sample and raw status filters",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   []string{"--process", "ALIGN", "--name", "sample-1", "--sample", "tumor", "--status", "failed"},
			want: TasksOptions{
				Global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
				Query: domain.TaskQuery{
					ProcessSubstring: "ALIGN",
					NameSubstring:    "sample-1",
					SampleSubstring:  "tumor",
					StatusRaw:        "failed",
				},
			},
		},
		{
			name:   "parses equals filters and command local run dir and format",
			global: GlobalOptions{RunDir: "/runs/inherited", Format: domain.OutputFormatHuman},
			args:   []string{"--run-dir=relative/run", "--format=json", "--process=CALL", "--name=sample-2", "--sample=normal", "--status=COMPLETED"},
			want: TasksOptions{
				Global: GlobalOptions{RunDir: "relative/run", Format: domain.OutputFormatJSON},
				Query: domain.TaskQuery{
					ProcessSubstring: "CALL",
					NameSubstring:    "sample-2",
					SampleSubstring:  "normal",
					StatusRaw:        "COMPLETED",
				},
			},
		},
		{
			name:   "json shorthand overrides inherited human format",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   []string{"--json"},
			want:   TasksOptions{Global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatJSON}},
		},
		{
			name:   "command local results dir equals parses before filters",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   []string{"--results-dir=relative/results", "--status", "FAILED"},
			want: TasksOptions{
				Global: GlobalOptions{RunDir: "/runs/nf", ResultsDir: "relative/results", Format: domain.OutputFormatHuman},
				Query:  domain.TaskQuery{StatusRaw: "FAILED"},
			},
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			got, err := ParseTasksOptions(test.global, test.args)
			if err != nil {
				t.Fatalf("ParseTasksOptions(%#v, %v) returned error: %v", test.global, test.args, err)
			}
			if got != test.want {
				t.Fatalf("tasks options = %#v, want %#v", got, test.want)
			}
		})
	}
}

func TestParseTasksOptionsRejectsInvalidTasksArgs(t *testing.T) {
	global := GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman}
	tests := []struct {
		name        string
		args        []string
		wantMessage string
	}{
		{name: "unexpected positional arg", args: []string{"extra"}, wantMessage: "unexpected tasks argument"},
		{name: "unexpected arg after separator", args: []string{"--", "extra"}, wantMessage: "unexpected tasks argument"},
		{name: "unknown flag", args: []string{"--verbose"}, wantMessage: "unknown tasks option"},
		{name: "missing process value", args: []string{"--process"}, wantMessage: "--process requires VALUE"},
		{name: "empty process equals", args: []string{"--process="}, wantMessage: "--process requires VALUE"},
		{name: "missing name value", args: []string{"--name"}, wantMessage: "--name requires VALUE"},
		{name: "empty sample equals", args: []string{"--sample="}, wantMessage: "--sample requires VALUE"},
		{name: "missing status value", args: []string{"--status"}, wantMessage: "--status requires STATUS"},
		{name: "empty status equals", args: []string{"--status="}, wantMessage: "--status requires STATUS"},
		{name: "unsupported format value", args: []string{"--format", "xml"}, wantMessage: "unsupported output format"},
		{name: "missing run dir", args: []string{"--run-dir"}, wantMessage: "--run-dir requires DIR"},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			got, err := ParseTasksOptions(global, test.args)
			if err == nil {
				t.Fatalf("ParseTasksOptions(%#v, %v) returned nil error", global, test.args)
			}
			if !strings.Contains(err.Error(), test.wantMessage) {
				t.Fatalf("error = %q, want it to contain %q", err.Error(), test.wantMessage)
			}
			if got != (TasksOptions{}) {
				t.Fatalf("tasks options on error = %#v, want zero value", got)
			}
		})
	}
}

func TestParseInspectOptionsParsesSelectorAndCommandLocalOptions(t *testing.T) {
	tests := []struct {
		name   string
		global GlobalOptions
		args   []string
		want   InspectOptions
	}{
		{
			name:   "inherits global options and records selector",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   []string{"ab/c123def"},
			want:   InspectOptions{Global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman}, Selector: "ab/c123def"},
		},
		{
			name:   "command local json shorthand overrides inherited human format",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   []string{"--json", "ab/c123def"},
			want:   InspectOptions{Global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatJSON}, Selector: "ab/c123def"},
		},
		{
			name:   "command local run dir and format flags override inherited options",
			global: GlobalOptions{RunDir: "/runs/inherited", Format: domain.OutputFormatHuman},
			args:   []string{"--run-dir", "relative/run", "--format", "json", "ab/c123def"},
			want:   InspectOptions{Global: GlobalOptions{RunDir: "relative/run", Format: domain.OutputFormatJSON}, Selector: "ab/c123def"},
		},
		{
			name:   "equals options and short run dir parse before selector",
			global: GlobalOptions{RunDir: "/runs/inherited", Format: domain.OutputFormatJSON},
			args:   []string{"--format=human", "-d", "short-run", "sample-1"},
			want:   InspectOptions{Global: GlobalOptions{RunDir: "short-run", Format: domain.OutputFormatHuman}, Selector: "sample-1"},
		},
		{
			name:   "command local results dir parses before selector",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   []string{"--results-dir", "relative/results", "sample-1"},
			want:   InspectOptions{Global: GlobalOptions{RunDir: "/runs/nf", ResultsDir: "relative/results", Format: domain.OutputFormatHuman}, Selector: "sample-1"},
		},
		{
			name:   "separator stops option parsing before selector",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   []string{"--json", "--", "-leading-dash-selector"},
			want:   InspectOptions{Global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatJSON}, Selector: "-leading-dash-selector"},
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			got, err := ParseInspectOptions(test.global, test.args)
			if err != nil {
				t.Fatalf("ParseInspectOptions(%#v, %v) returned error: %v", test.global, test.args, err)
			}
			if got != test.want {
				t.Fatalf("inspect options = %#v, want %#v", got, test.want)
			}
		})
	}
}

func TestParseInspectOptionsRejectsInvalidInspectArgs(t *testing.T) {
	global := GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman}
	tests := []struct {
		name        string
		args        []string
		wantMessage string
	}{
		{name: "missing selector", args: nil, wantMessage: "inspect requires exactly one selector"},
		{name: "multiple selectors", args: []string{"ab/c123def", "cd/e456ghi"}, wantMessage: "inspect requires exactly one selector"},
		{name: "multiple selectors after separator", args: []string{"--", "ab/c123def", "cd/e456ghi"}, wantMessage: "inspect requires exactly one selector"},
		{name: "unknown flag", args: []string{"--verbose", "ab/c123def"}, wantMessage: "unknown inspect option"},
		{name: "missing run dir", args: []string{"--run-dir"}, wantMessage: "--run-dir requires DIR"},
		{name: "empty run dir equals", args: []string{"--run-dir="}, wantMessage: "--run-dir requires DIR"},
		{name: "missing short run dir", args: []string{"-d"}, wantMessage: "-d requires DIR"},
		{name: "missing format value", args: []string{"--format"}, wantMessage: "--format requires FORMAT"},
		{name: "empty format equals", args: []string{"--format="}, wantMessage: "--format requires FORMAT"},
		{name: "unsupported format value", args: []string{"--format", "xml", "ab/c123def"}, wantMessage: "unsupported output format"},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			got, err := ParseInspectOptions(global, test.args)
			if err == nil {
				t.Fatalf("ParseInspectOptions(%#v, %v) returned nil error", global, test.args)
			}
			if !strings.Contains(err.Error(), test.wantMessage) {
				t.Fatalf("error = %q, want it to contain %q", err.Error(), test.wantMessage)
			}
			if got != (InspectOptions{}) {
				t.Fatalf("inspect options on error = %#v, want zero value", got)
			}
		})
	}
}

func TestParseIndexOptionsParsesRefreshAndCommandLocalOptions(t *testing.T) {
	tests := []struct {
		name   string
		global GlobalOptions
		args   []string
		want   IndexOptions
	}{
		{
			name:   "inherits run dir and human format with refresh disabled",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   nil,
			want:   IndexOptions{Global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman}, Refresh: false},
		},
		{
			name:   "parses refresh flag",
			global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman},
			args:   []string{"--refresh"},
			want:   IndexOptions{Global: GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman}, Refresh: true},
		},
		{
			name:   "command local run dir and format override inherited options",
			global: GlobalOptions{RunDir: "/runs/inherited", Format: domain.OutputFormatHuman},
			args:   []string{"--run-dir", "relative/run", "--format", "json", "--refresh"},
			want:   IndexOptions{Global: GlobalOptions{RunDir: "relative/run", Format: domain.OutputFormatJSON}, Refresh: true},
		},
		{
			name:   "equals options short run dir and json shorthand parse",
			global: GlobalOptions{RunDir: "/runs/inherited", Format: domain.OutputFormatJSON},
			args:   []string{"--format=human", "-d", "short-run", "--json", "--run-dir=final-run"},
			want:   IndexOptions{Global: GlobalOptions{RunDir: "final-run", Format: domain.OutputFormatJSON}, Refresh: false},
		},
		{
			name:   "command local results dir equals parses with refresh",
			global: GlobalOptions{RunDir: "/runs/inherited", Format: domain.OutputFormatHuman},
			args:   []string{"--results-dir=../results", "--refresh"},
			want:   IndexOptions{Global: GlobalOptions{RunDir: "/runs/inherited", ResultsDir: "../results", Format: domain.OutputFormatHuman}, Refresh: true},
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			got, err := ParseIndexOptions(test.global, test.args)
			if err != nil {
				t.Fatalf("ParseIndexOptions(%#v, %v) returned error: %v", test.global, test.args, err)
			}
			if got != test.want {
				t.Fatalf("index options = %#v, want %#v", got, test.want)
			}
		})
	}
}

func TestParseIndexOptionsRejectsInvalidIndexArgs(t *testing.T) {
	global := GlobalOptions{RunDir: "/runs/nf", Format: domain.OutputFormatHuman}
	tests := []struct {
		name        string
		args        []string
		wantMessage string
	}{
		{name: "unexpected positional arg", args: []string{"extra"}, wantMessage: "unexpected index argument"},
		{name: "unexpected arg after separator", args: []string{"--", "extra"}, wantMessage: "unexpected index argument"},
		{name: "unknown flag", args: []string{"--verbose"}, wantMessage: "unknown index option"},
		{name: "refresh does not take equals value", args: []string{"--refresh=true"}, wantMessage: "unknown index option"},
		{name: "missing run dir", args: []string{"--run-dir"}, wantMessage: "--run-dir requires DIR"},
		{name: "empty run dir equals", args: []string{"--run-dir="}, wantMessage: "--run-dir requires DIR"},
		{name: "missing short run dir", args: []string{"-d"}, wantMessage: "-d requires DIR"},
		{name: "missing format value", args: []string{"--format"}, wantMessage: "--format requires FORMAT"},
		{name: "empty format equals", args: []string{"--format="}, wantMessage: "--format requires FORMAT"},
		{name: "unsupported format value", args: []string{"--format", "xml"}, wantMessage: "unsupported output format"},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			got, err := ParseIndexOptions(global, test.args)
			if err == nil {
				t.Fatalf("ParseIndexOptions(%#v, %v) returned nil error", global, test.args)
			}
			if !strings.Contains(err.Error(), test.wantMessage) {
				t.Fatalf("error = %q, want it to contain %q", err.Error(), test.wantMessage)
			}
			if got != (IndexOptions{}) {
				t.Fatalf("index options on error = %#v, want zero value", got)
			}
		})
	}
}

func TestRunIndexWithoutRefreshReportsMissingIndexWithoutCreatingCache(t *testing.T) {
	runDir := t.TempDir()
	tracePath := filepath.Join(runDir, "trace.csv")
	if err := os.WriteFile(tracePath, []byte(""), 0o644); err != nil {
		t.Fatalf("write trace fixture: %v", err)
	}

	var output strings.Builder
	err := RunIndex(context.Background(), IndexOptions{
		Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
	}, &output)
	if err != nil {
		t.Fatalf("RunIndex(no refresh) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{"mode: trace-backed", "freshness: missing", "index_missing", tracePath} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunIndex(no refresh) output = %q, want it to contain %q", got, want)
		}
	}
	assertNoIndexCacheDir(t, runDir)
}

func TestRunIndexWithoutRefreshUsesCustomResultsDirPipelineInfoTrace(t *testing.T) {
	runDir := t.TempDir()
	resultsDir := filepath.Join(runDir, "custom-results")
	pipelineInfoRoot := filepath.Join(resultsDir, run.PipelineInfoDirName)
	if err := os.MkdirAll(pipelineInfoRoot, 0o755); err != nil {
		t.Fatalf("mkdir pipeline_info fixture: %v", err)
	}

	logPath := filepath.Join(runDir, ".nextflow.log")
	if err := os.WriteFile(logPath, []byte("ERROR ~ log-only decoy should not drive index diagnostics\n"), 0o644); err != nil {
		t.Fatalf("write log-only decoy fixture: %v", err)
	}

	tracePath := filepath.Join(pipelineInfoRoot, "execution_trace_2026-04-29.csv")
	if err := os.WriteFile(tracePath, []byte("not parsed by diagnostics-only index report\n"), 0o644); err != nil {
		t.Fatalf("write pipeline-info trace fixture: %v", err)
	}

	var output strings.Builder
	err := RunIndex(context.Background(), IndexOptions{
		Global: GlobalOptions{RunDir: runDir, ResultsDir: "custom-results", Format: domain.OutputFormatHuman},
	}, &output)
	if err != nil {
		t.Fatalf("RunIndex(no refresh custom results dir) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		"mode: trace-backed",
		"freshness: missing",
		"index_missing",
		tracePath,
		pipelineInfoRoot,
		"pipeline_info execution trace files",
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunIndex(no refresh custom results dir) output = %q, want it to contain %q", got, want)
		}
	}
	for _, notWant := range []string{"mode: log-only", "trace: none"} {
		if strings.Contains(got, notWant) {
			t.Fatalf("RunIndex(no refresh custom results dir) output = %q, did not want %q", got, notWant)
		}
	}
	assertNoIndexCacheDir(t, runDir)
}

func TestRunIndexRefreshUsesCustomResultsDirPipelineInfoTrace(t *testing.T) {
	runDir := t.TempDir()
	resultsDir := filepath.Join(runDir, "custom-results")
	pipelineInfoRoot := filepath.Join(resultsDir, run.PipelineInfoDirName)
	failedWorkdir := filepath.Join(runDir, "work", "dd", "444444")
	for _, dir := range []string{pipelineInfoRoot, failedWorkdir} {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			t.Fatalf("mkdir fixture directory %q: %v", dir, err)
		}
	}

	logPath := filepath.Join(runDir, ".nextflow.log")
	if err := os.WriteFile(logPath, []byte("ERROR ~ log-only decoy should not drive index refresh\n"), 0o644); err != nil {
		t.Fatalf("write log-only decoy fixture: %v", err)
	}

	tracePath := filepath.Join(pipelineInfoRoot, "execution_trace_2026-04-29.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"DD/444444,FAILED,TRACE:CALL,TRACE:CALL (tumor-01),tumor-01," + failedWorkdir + ",1,5m,300s,2,6 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write pipeline-info trace fixture: %v", err)
	}

	var output strings.Builder
	err := RunIndex(context.Background(), IndexOptions{
		Global:  GlobalOptions{RunDir: runDir, ResultsDir: "custom-results", Format: domain.OutputFormatHuman},
		Refresh: true,
	}, &output)
	if err != nil {
		t.Fatalf("RunIndex(refresh custom results dir) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{
		"mode: trace-backed",
		"freshness: fresh",
		"task_count: 1",
		tracePath,
		pipelineInfoRoot,
		"pipeline_info execution trace files",
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunIndex(refresh custom results dir) output = %q, want it to contain %q", got, want)
		}
	}
	for _, notWant := range []string{"mode: log-only", "trace: none"} {
		if strings.Contains(got, notWant) {
			t.Fatalf("RunIndex(refresh custom results dir) output = %q, did not want %q", got, notWant)
		}
	}

	store, err := indexdb.OpenStore(context.Background(), domain.RunDir{Path: runDir})
	if err != nil {
		t.Fatalf("open store after RunIndex custom results dir refresh: %v", err)
	}
	defer store.Close()

	metadata, err := indexdb.ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata after RunIndex custom results dir refresh: %v", err)
	}
	if metadata == nil {
		t.Fatalf("ReadMetadata after RunIndex custom results dir refresh = nil, want metadata")
	}
	if metadata.Mode != domain.IndexModeTraceBacked || metadata.Freshness != domain.IndexFreshnessFresh || metadata.TaskCount != 1 {
		t.Fatalf("metadata mode/freshness/task_count = %q/%q/%d, want %q/%q/1", metadata.Mode, metadata.Freshness, metadata.TaskCount, domain.IndexModeTraceBacked, domain.IndexFreshnessFresh)
	}
	if metadata.Trace == nil || metadata.Trace.Path != tracePath {
		t.Fatalf("metadata trace = %#v, want path %q", metadata.Trace, tracePath)
	}

	tasks, err := indexdb.QueryTasks(context.Background(), store, domain.TaskQuery{})
	if err != nil {
		t.Fatalf("QueryTasks after RunIndex custom results dir refresh: %v", err)
	}
	if len(tasks) != 1 || tasks[0].ID != "dd/444444" || tasks[0].Workdir != failedWorkdir {
		t.Fatalf("task rows after RunIndex custom results dir refresh = %#v, want one dd/444444 task in %q", tasks, failedWorkdir)
	}
}

func TestRunIndexRefreshTraceBackedRebuildsTasksAndReportsRealCount(t *testing.T) {
	runDir, failedID := writeExecuteTraceFixture(t)
	tracePath := filepath.Join(runDir, "trace.csv")

	var output strings.Builder
	err := RunIndex(context.Background(), IndexOptions{
		Global:  GlobalOptions{RunDir: runDir, Format: domain.OutputFormatJSON},
		Refresh: true,
	}, &output)
	if err != nil {
		t.Fatalf("RunIndex(refresh trace-backed) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{`"format": "json"`, `"mode": "trace-backed"`, `"freshness": "fresh"`, `"task_count": 2`, tracePath} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunIndex(refresh trace-backed) output = %q, want it to contain %q", got, want)
		}
	}

	store, err := indexdb.OpenStore(context.Background(), domain.RunDir{Path: runDir})
	if err != nil {
		t.Fatalf("open store after RunIndex refresh: %v", err)
	}
	defer store.Close()

	metadata, err := indexdb.ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata after RunIndex refresh: %v", err)
	}
	if metadata == nil {
		t.Fatalf("ReadMetadata after RunIndex refresh = nil, want metadata")
	}
	if metadata.Mode != domain.IndexModeTraceBacked || metadata.Freshness != domain.IndexFreshnessFresh || metadata.TaskCount != 2 {
		t.Fatalf("metadata mode/freshness/task_count = %q/%q/%d, want %q/%q/2", metadata.Mode, metadata.Freshness, metadata.TaskCount, domain.IndexModeTraceBacked, domain.IndexFreshnessFresh)
	}
	if metadata.Trace == nil || metadata.Trace.Path != tracePath {
		t.Fatalf("metadata trace = %#v, want path %q", metadata.Trace, tracePath)
	}

	tasks, err := indexdb.QueryTasks(context.Background(), store, domain.TaskQuery{})
	if err != nil {
		t.Fatalf("QueryTasks after RunIndex refresh: %v", err)
	}
	if len(tasks) != 2 {
		t.Fatalf("task rows after RunIndex refresh = %d, want 2: %#v", len(tasks), tasks)
	}
	if tasks[0].ID != "aa/111111" || tasks[1].ID != failedID {
		t.Fatalf("task ids after RunIndex refresh = %q, %q; want %q, %q", tasks[0].ID, tasks[1].ID, "aa/111111", failedID)
	}
}

func TestRunIndexRefreshTraceBackedSurfacesMalformedTrace(t *testing.T) {
	runDir := t.TempDir()
	tracePath := filepath.Join(runDir, "trace.csv")
	if err := os.WriteFile(tracePath, []byte(""), 0o644); err != nil {
		t.Fatalf("write malformed trace fixture: %v", err)
	}

	var output strings.Builder
	err := RunIndex(context.Background(), IndexOptions{
		Global:  GlobalOptions{RunDir: runDir, Format: domain.OutputFormatJSON},
		Refresh: true,
	}, &output)
	if err == nil {
		t.Fatalf("RunIndex(refresh malformed trace-backed) returned nil error; output = %q", output.String())
	}
	if !strings.Contains(err.Error(), "missing header row") {
		t.Fatalf("RunIndex(refresh malformed trace-backed) error = %q, want missing header row", err.Error())
	}
}

func TestRunIndexRefreshLogOnlyWritesMetadataOnly(t *testing.T) {
	runDir := t.TempDir()
	logPath := filepath.Join(runDir, ".nextflow.log")
	if err := os.WriteFile(logPath, []byte("Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start\n"), 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}

	var output strings.Builder
	err := RunIndex(context.Background(), IndexOptions{
		Global:  GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
		Refresh: true,
	}, &output)
	if err != nil {
		t.Fatalf("RunIndex(refresh log-only metadata) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{"mode: log-only", "freshness: fresh", "task_count: 0", logPath} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunIndex(refresh log-only metadata) output = %q, want it to contain %q", got, want)
		}
	}

	store, err := indexdb.OpenStore(context.Background(), domain.RunDir{Path: runDir})
	if err != nil {
		t.Fatalf("open store after RunIndex log-only refresh: %v", err)
	}
	defer store.Close()

	metadata, err := indexdb.ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata after RunIndex log-only refresh: %v", err)
	}
	if metadata == nil {
		t.Fatalf("ReadMetadata after RunIndex log-only refresh = nil, want metadata")
	}
	if metadata.Mode != domain.IndexModeLogOnly || metadata.Freshness != domain.IndexFreshnessFresh || metadata.TaskCount != 0 {
		t.Fatalf("metadata mode/freshness/task_count = %q/%q/%d, want %q/%q/0", metadata.Mode, metadata.Freshness, metadata.TaskCount, domain.IndexModeLogOnly, domain.IndexFreshnessFresh)
	}
	if metadata.Log == nil || metadata.Log.Path != logPath {
		t.Fatalf("metadata log = %#v, want path %q", metadata.Log, logPath)
	}
}

func TestRunIndexUnsupportedWithoutRefreshReportsDiagnosticsWithoutCreatingCache(t *testing.T) {
	runDir := t.TempDir()

	var output strings.Builder
	err := RunIndex(context.Background(), IndexOptions{
		Global: GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
	}, &output)
	if err != nil {
		t.Fatalf("RunIndex(unsupported no refresh) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{"unsupported_artifacts", "Use `nextflow run ... -with-trace`"} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunIndex(unsupported no refresh) output = %q, want it to contain %q", got, want)
		}
	}
	assertNoIndexCacheDir(t, runDir)
}

func TestRunIndexRefreshWritesUnsupportedMetadata(t *testing.T) {
	runDir := t.TempDir()

	var output strings.Builder
	err := RunIndex(context.Background(), IndexOptions{
		Global:  GlobalOptions{RunDir: runDir, Format: domain.OutputFormatHuman},
		Refresh: true,
	}, &output)
	if err != nil {
		t.Fatalf("RunIndex(refresh unsupported metadata) returned error: %v", err)
	}

	got := output.String()
	for _, want := range []string{"mode: unsupported", "freshness: unsupported", "unsupported_artifacts"} {
		if !strings.Contains(got, want) {
			t.Fatalf("RunIndex(refresh unsupported metadata) output = %q, want it to contain %q", got, want)
		}
	}
	if _, err := os.Stat(run.IndexPath(domain.RunDir{Path: runDir})); err != nil {
		t.Fatalf("stat index after unsupported refresh: %v", err)
	}

	store, err := indexdb.OpenStore(context.Background(), domain.RunDir{Path: runDir})
	if err != nil {
		t.Fatalf("open store after unsupported refresh: %v", err)
	}
	defer store.Close()

	metadata, err := indexdb.ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata after unsupported refresh: %v", err)
	}
	if metadata == nil {
		t.Fatalf("ReadMetadata after unsupported refresh = nil, want metadata")
	}
	if metadata.Mode != domain.IndexModeUnsupported || metadata.Freshness != domain.IndexFreshnessUnsupported || metadata.StaleReason != "no supported artifacts" {
		t.Fatalf("metadata mode/freshness/stale_reason = %q/%q/%q, want %q/%q/%q", metadata.Mode, metadata.Freshness, metadata.StaleReason, domain.IndexModeUnsupported, domain.IndexFreshnessUnsupported, "no supported artifacts")
	}
}

func writeExecuteTraceFixture(t *testing.T) (string, string) {
	t.Helper()

	runDir := t.TempDir()
	completedWorkdir := filepath.Join(runDir, "work", "aa", "111111")
	failedWorkdir := filepath.Join(runDir, "work", "bb", "222222")
	for _, dir := range []string{completedWorkdir, failedWorkdir} {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			t.Fatalf("create workdir fixture %s: %v", dir, err)
		}
	}
	if err := os.WriteFile(filepath.Join(failedWorkdir, ".command.sh"), []byte("#!/usr/bin/env bash\necho failed\n"), 0o644); err != nil {
		t.Fatalf("write command fixture: %v", err)
	}

	tracePath := filepath.Join(runDir, "trace.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AA/111111,COMPLETED,ALIGN,ALIGN (sample-ok),sample-ok," + completedWorkdir + ",0,1m,60s,2,4 GB",
		"BB/222222,FAILED,QUANT,QUANT (sample-failed),sample-failed," + failedWorkdir + ",1,2m,120s,4,8 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write trace fixture: %v", err)
	}

	return runDir, "bb/222222"
}

func writeExecutePipelineInfoTraceFixture(t *testing.T) (string, string, string, string) {
	t.Helper()

	runDir := t.TempDir()
	resultsDir := filepath.Join(runDir, "custom-results")
	pipelineInfoRoot := filepath.Join(resultsDir, run.PipelineInfoDirName)
	failedWorkdir := filepath.Join(runDir, "work", "cc", "333333")
	for _, dir := range []string{pipelineInfoRoot, failedWorkdir} {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			t.Fatalf("create pipeline-info fixture dir %s: %v", dir, err)
		}
	}
	if err := os.WriteFile(filepath.Join(failedWorkdir, ".command.sh"), []byte("#!/usr/bin/env bash\necho dispatch failed\n"), 0o644); err != nil {
		t.Fatalf("write command fixture: %v", err)
	}
	if err := os.WriteFile(filepath.Join(runDir, ".nextflow.log"), []byte("ERROR ~ Error executing process > 'LOG:ONLY (should-not-drive-dispatch)'\n"), 0o644); err != nil {
		t.Fatalf("write log fallback fixture: %v", err)
	}

	tracePath := filepath.Join(pipelineInfoRoot, "execution_trace_dispatch.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"CC/333333,FAILED,TRACE:DISPATCH,TRACE:DISPATCH (case-01),case-01," + failedWorkdir + ",1,3m,180s,2,5 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write pipeline-info trace fixture: %v", err)
	}

	return runDir, resultsDir, "cc/333333", tracePath
}

func assertNoIndexCacheDir(t *testing.T, runDir string) {
	t.Helper()
	if _, err := os.Stat(filepath.Join(runDir, run.IndexDirName)); !os.IsNotExist(err) {
		t.Fatalf("cache directory stat error = %v, want cache directory not created", err)
	}
}

func TestParseGlobalOptionsRejectsInvalidGlobalFlags(t *testing.T) {
	tests := []struct {
		name        string
		args        []string
		wantMessage string
	}{
		{name: "missing run dir after long flag", args: []string{"--run-dir"}, wantMessage: "--run-dir requires DIR"},
		{name: "empty run dir after equals", args: []string{"--run-dir="}, wantMessage: "--run-dir requires DIR"},
		{name: "missing run dir after short flag", args: []string{"-d"}, wantMessage: "-d requires DIR"},
		{name: "missing results dir after long flag", args: []string{"--results-dir"}, wantMessage: "--results-dir requires DIR"},
		{name: "empty results dir after equals", args: []string{"--results-dir="}, wantMessage: "--results-dir requires DIR"},
		{name: "missing format value", args: []string{"--format"}, wantMessage: "--format requires FORMAT"},
		{name: "unsupported format value", args: []string{"--format", "xml"}, wantMessage: "unsupported output format"},
		{name: "unknown global flag before command", args: []string{"--verbose", "status"}, wantMessage: "unknown global option"},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			gotGlobal, gotRemaining, err := ParseGlobalOptions(test.args)
			if err == nil {
				t.Fatalf("ParseGlobalOptions(%v) returned nil error", test.args)
			}
			if !strings.Contains(err.Error(), test.wantMessage) {
				t.Fatalf("error = %q, want it to contain %q", err.Error(), test.wantMessage)
			}
			if gotGlobal != (GlobalOptions{}) {
				t.Fatalf("global options on error = %#v, want zero value", gotGlobal)
			}
			if gotRemaining != nil {
				t.Fatalf("remaining args on error = %#v, want nil", gotRemaining)
			}
		})
	}
}
