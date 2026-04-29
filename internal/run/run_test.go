package run

import (
	"context"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/mskilab-org/gosh/internal/domain"
)

func TestResolveRunDirUsesCurrentDirectoryForDefault(t *testing.T) {
	runRoot := t.TempDir()
	t.Chdir(runRoot)

	resolved, err := ResolveRunDir(DefaultRunDir)
	if err != nil {
		t.Fatalf("ResolveRunDir(%q) returned error: %v", DefaultRunDir, err)
	}

	expected, err := filepath.Abs(runRoot)
	if err != nil {
		t.Fatalf("filepath.Abs(%q): %v", runRoot, err)
	}
	if resolved.Path != expected {
		t.Fatalf("resolved path = %q, want %q", resolved.Path, expected)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestResolveRunDirResolvesRelativeDirectory(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "runs", "nf-run")
	if err := os.MkdirAll(runRoot, 0o755); err != nil {
		t.Fatalf("mkdir run root: %v", err)
	}
	t.Chdir(workspace)

	resolved, err := ResolveRunDir(filepath.Join("runs", "..", "runs", "nf-run"))
	if err != nil {
		t.Fatalf("ResolveRunDir(relative run dir) returned error: %v", err)
	}

	expected, err := filepath.Abs(runRoot)
	if err != nil {
		t.Fatalf("filepath.Abs(%q): %v", runRoot, err)
	}
	if resolved.Path != expected {
		t.Fatalf("resolved path = %q, want %q", resolved.Path, expected)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestResolveRunDirRejectsMissingPath(t *testing.T) {
	missing := filepath.Join(t.TempDir(), "does-not-exist")

	resolved, err := ResolveRunDir(missing)
	if err == nil {
		t.Fatalf("ResolveRunDir(%q) returned nil error", missing)
	}
	if resolved.Path != "" {
		t.Fatalf("resolved path on error = %q, want empty", resolved.Path)
	}
}

func TestResolveRunDirRejectsRegularFile(t *testing.T) {
	workspace := t.TempDir()
	filePath := filepath.Join(workspace, "trace.txt")
	if err := os.WriteFile(filePath, []byte("task_id\n"), 0o644); err != nil {
		t.Fatalf("write file fixture: %v", err)
	}

	resolved, err := ResolveRunDir(filePath)
	if err == nil {
		t.Fatalf("ResolveRunDir(%q) returned nil error", filePath)
	}
	if !strings.Contains(err.Error(), "not a directory") {
		t.Fatalf("error = %q, want it to mention not a directory", err.Error())
	}
	if resolved.Path != "" {
		t.Fatalf("resolved path on error = %q, want empty", resolved.Path)
	}
	assertNoGoshDirectory(t, workspace)
}

func TestResolveResultsDirDefaultsUnderRunDir(t *testing.T) {
	runRoot := filepath.Join(t.TempDir(), "nf-run")
	if err := os.MkdirAll(runRoot, 0o755); err != nil {
		t.Fatalf("mkdir run root: %v", err)
	}

	resolved, err := ResolveResultsDir(domain.RunDir{Path: runRoot}, "")
	if err != nil {
		t.Fatalf("ResolveResultsDir(default) returned error: %v", err)
	}

	want := filepath.Join(runRoot, DefaultResultsDirName)
	if resolved.Path != want {
		t.Fatalf("resolved path = %q, want %q", resolved.Path, want)
	}
	if _, err := os.Stat(want); !os.IsNotExist(err) {
		t.Fatalf("ResolveResultsDir(default) should not create %q; stat err = %v", want, err)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestResolveResultsDirResolvesRelativeInputAgainstRunDir(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "runs", "nf-run")
	cwd := filepath.Join(workspace, "cwd")
	if err := os.MkdirAll(runRoot, 0o755); err != nil {
		t.Fatalf("mkdir run root: %v", err)
	}
	if err := os.MkdirAll(cwd, 0o755); err != nil {
		t.Fatalf("mkdir cwd fixture: %v", err)
	}
	t.Chdir(cwd)

	input := filepath.Join("custom", "..", "external", "results")
	resolved, err := ResolveResultsDir(domain.RunDir{Path: runRoot}, input)
	if err != nil {
		t.Fatalf("ResolveResultsDir(relative) returned error: %v", err)
	}

	want := filepath.Clean(filepath.Join(runRoot, input))
	if resolved.Path != want {
		t.Fatalf("resolved path = %q, want %q", resolved.Path, want)
	}
	cwdRelative := filepath.Clean(filepath.Join(cwd, input))
	if resolved.Path == cwdRelative {
		t.Fatalf("resolved path = %q, want relative input resolved against run dir instead of cwd %q", resolved.Path, cwdRelative)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestResolveResultsDirPreservesAbsoluteInput(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "runs", "nf-run")
	if err := os.MkdirAll(runRoot, 0o755); err != nil {
		t.Fatalf("mkdir run root: %v", err)
	}

	separator := string(os.PathSeparator)
	input := workspace + separator + "external" + separator + "." + separator + "results" + separator + ".." + separator + "chosen"
	resolved, err := ResolveResultsDir(domain.RunDir{Path: runRoot}, input)
	if err != nil {
		t.Fatalf("ResolveResultsDir(absolute) returned error: %v", err)
	}

	want := filepath.Clean(input)
	if resolved.Path != want {
		t.Fatalf("resolved path = %q, want %q", resolved.Path, want)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestResolveResultsDirRejectsEmptyRunDirForDefault(t *testing.T) {
	resolved, err := ResolveResultsDir(domain.RunDir{}, "")
	if err == nil {
		t.Fatalf("ResolveResultsDir(empty run dir) returned nil error")
	}
	if !strings.Contains(err.Error(), "empty run dir") {
		t.Fatalf("error = %q, want it to mention empty run dir", err.Error())
	}
	if resolved.Path != "" {
		t.Fatalf("resolved path on error = %q, want empty", resolved.Path)
	}
}

func TestBuildArtifactSearchLocationsDescribesRootAndPipelineInfoSearches(t *testing.T) {
	runRoot := filepath.Join(t.TempDir(), "nf-run")
	resultsRoot := filepath.Join(runRoot, DefaultResultsDirName)

	got := BuildArtifactSearchLocations(domain.RunDir{Path: runRoot}, domain.ResultsDir{Path: resultsRoot})

	if len(got) != 3 {
		t.Fatalf("location count = %d, want 3: %+v", len(got), got)
	}
	assertArtifactSearchLocation(t, got[0], domain.SourceKindTrace, runRoot, TracePatterns, "run directory", "trace")
	assertArtifactSearchLocation(t, got[1], domain.SourceKindTrace, filepath.Join(resultsRoot, PipelineInfoDirName), PipelineInfoTracePatterns, "pipeline_info", "execution trace")
	assertArtifactSearchLocation(t, got[2], domain.SourceKindLog, runRoot, LogPatterns, "run directory", "log")
}

func TestBuildArtifactSearchLocationsUsesConfiguredResultsDirForPipelineInfo(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "runs", "nf-run")
	resultsRoot := filepath.Join(workspace, "external", "nf-results")

	got := BuildArtifactSearchLocations(domain.RunDir{Path: runRoot}, domain.ResultsDir{Path: resultsRoot})

	if len(got) != 3 {
		t.Fatalf("location count = %d, want 3: %+v", len(got), got)
	}
	if got[0].BaseDir != runRoot || got[2].BaseDir != runRoot {
		t.Fatalf("root search base dirs = %q and %q, want both %q", got[0].BaseDir, got[2].BaseDir, runRoot)
	}
	wantPipelineInfo := filepath.Join(resultsRoot, PipelineInfoDirName)
	if got[1].BaseDir != wantPipelineInfo {
		t.Fatalf("pipeline-info search base dir = %q, want %q", got[1].BaseDir, wantPipelineInfo)
	}
}

func TestBuildArtifactSearchLocationsReturnsIndependentPatternSlices(t *testing.T) {
	originalTracePatterns := append([]string{}, TracePatterns...)
	originalPipelineInfoTracePatterns := append([]string{}, PipelineInfoTracePatterns...)
	originalLogPatterns := append([]string{}, LogPatterns...)
	t.Cleanup(func() {
		TracePatterns = originalTracePatterns
		PipelineInfoTracePatterns = originalPipelineInfoTracePatterns
		LogPatterns = originalLogPatterns
	})

	runRoot := filepath.Join(t.TempDir(), "nf-run")
	resultsRoot := filepath.Join(runRoot, DefaultResultsDirName)
	got := BuildArtifactSearchLocations(domain.RunDir{Path: runRoot}, domain.ResultsDir{Path: resultsRoot})
	if len(got) != 3 {
		t.Fatalf("location count = %d, want 3: %+v", len(got), got)
	}

	got[0].Patterns[0] = "mutated-root-trace-pattern"
	got[1].Patterns[0] = "mutated-pipeline-info-pattern"
	got[2].Patterns[0] = "mutated-log-pattern"

	assertPatterns(t, TracePatterns, originalTracePatterns)
	assertPatterns(t, PipelineInfoTracePatterns, originalPipelineInfoTracePatterns)
	assertPatterns(t, LogPatterns, originalLogPatterns)

	next := BuildArtifactSearchLocations(domain.RunDir{Path: runRoot}, domain.ResultsDir{Path: resultsRoot})
	assertPatterns(t, next[0].Patterns, originalTracePatterns)
	assertPatterns(t, next[1].Patterns, originalPipelineInfoTracePatterns)
	assertPatterns(t, next[2].Patterns, originalLogPatterns)
}

func TestFindCandidateSourcesInLocationsReturnsSortedFingerprintsAcrossLocations(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "run")
	resultsRoot := filepath.Join(workspace, "results")
	pipelineInfoRoot := filepath.Join(resultsRoot, PipelineInfoDirName)
	if err := os.MkdirAll(runRoot, 0o755); err != nil {
		t.Fatalf("mkdir run root: %v", err)
	}
	if err := os.MkdirAll(pipelineInfoRoot, 0o755); err != nil {
		t.Fatalf("mkdir pipeline_info fixture dir: %v", err)
	}

	base := time.Unix(1_700_000_000, 0)
	rootTrace := writeSourceFixture(t, runRoot, "trace-root.txt", []byte("root\n"), base.Add(2*time.Second))
	pipelineTrace := writeSourceFixture(t, pipelineInfoRoot, "execution_trace_2026-02-23.txt", []byte("pipeline\n"), base.Add(3*time.Second))
	log := writeSourceFixture(t, runRoot, ".nextflow.log", []byte("log\n"), base.Add(time.Second))
	writeSourceFixture(t, pipelineInfoRoot, "trace-not-matched.txt", []byte("ignored\n"), base.Add(4*time.Second))

	locations := BuildArtifactSearchLocations(domain.RunDir{Path: runRoot}, domain.ResultsDir{Path: resultsRoot})
	got, err := FindCandidateSourcesInLocations(locations)
	if err != nil {
		t.Fatalf("FindCandidateSourcesInLocations(configured locations) returned error: %v", err)
	}

	want := []domain.SourceFingerprint{
		{Kind: domain.SourceKindTrace, Path: pipelineTrace, ModTime: base.Add(3 * time.Second), Size: int64(len("pipeline\n"))},
		{Kind: domain.SourceKindLog, Path: log, ModTime: base.Add(time.Second), Size: int64(len("log\n"))},
		{Kind: domain.SourceKindTrace, Path: rootTrace, ModTime: base.Add(2 * time.Second), Size: int64(len("root\n"))},
	}
	assertSourceFingerprints(t, got, want)
	assertNoGoshDirectory(t, runRoot)
}

func TestFindCandidateSourcesInLocationsReturnsEmptyForNoLocations(t *testing.T) {
	got, err := FindCandidateSourcesInLocations(nil)
	if err != nil {
		t.Fatalf("FindCandidateSourcesInLocations(nil) returned error: %v", err)
	}
	if len(got) != 0 {
		t.Fatalf("candidate count = %d, want 0: %+v", len(got), got)
	}
}

func TestFindCandidateSourcesInLocationsDeduplicatesCleanPaths(t *testing.T) {
	runRoot := t.TempDir()
	base := time.Unix(1_700_000_000, 0)
	trace := writeSourceFixture(t, runRoot, "trace-duplicate.txt", []byte("trace\n"), base)

	locations := []domain.ArtifactSearchLocation{
		{Kind: domain.SourceKindTrace, BaseDir: runRoot, Patterns: []string{"trace*.txt"}, Description: "root trace files"},
		{Kind: domain.SourceKindTrace, BaseDir: filepath.Join(runRoot, "."), Patterns: []string{"trace-duplicate.txt"}, Description: "same trace via clean path"},
	}
	got, err := FindCandidateSourcesInLocations(locations)
	if err != nil {
		t.Fatalf("FindCandidateSourcesInLocations(duplicate paths) returned error: %v", err)
	}

	want := []domain.SourceFingerprint{
		{Kind: domain.SourceKindTrace, Path: trace, ModTime: base, Size: int64(len("trace\n"))},
	}
	assertSourceFingerprints(t, got, want)
	assertNoGoshDirectory(t, runRoot)
}

func TestFindCandidateSourcesInLocationsIgnoresNonRegularMatches(t *testing.T) {
	runRoot := t.TempDir()
	base := time.Unix(1_700_000_000, 0)
	trace := writeSourceFixture(t, runRoot, "trace-ok.txt", []byte("ok\n"), base)
	matchingDir := filepath.Join(runRoot, "trace-dir.txt")
	if err := os.Mkdir(matchingDir, 0o755); err != nil {
		t.Fatalf("mkdir matching directory fixture: %v", err)
	}

	locations := []domain.ArtifactSearchLocation{
		{Kind: domain.SourceKindTrace, BaseDir: runRoot, Patterns: []string{"trace*.txt"}, Description: "root trace files"},
	}
	got, err := FindCandidateSourcesInLocations(locations)
	if err != nil {
		t.Fatalf("FindCandidateSourcesInLocations(with matching directory) returned error: %v", err)
	}

	want := []domain.SourceFingerprint{
		{Kind: domain.SourceKindTrace, Path: trace, ModTime: base, Size: int64(len("ok\n"))},
	}
	assertSourceFingerprints(t, got, want)
	assertNoGoshDirectory(t, runRoot)
}

func TestFindCandidateSourcesInLocationsRejectsBadPattern(t *testing.T) {
	runRoot := t.TempDir()
	locations := []domain.ArtifactSearchLocation{
		{Kind: domain.SourceKindTrace, BaseDir: runRoot, Patterns: []string{"["}, Description: "broken trace files"},
	}

	got, err := FindCandidateSourcesInLocations(locations)
	if err == nil {
		t.Fatalf("FindCandidateSourcesInLocations(bad pattern) returned nil error")
	}
	if !strings.Contains(err.Error(), "glob source pattern") {
		t.Fatalf("error = %q, want it to mention glob source pattern", err.Error())
	}
	if len(got) != 0 {
		t.Fatalf("candidates on error = %+v, want none", got)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestIndexPathComputesRunLocalSQLitePath(t *testing.T) {
	runRoot := filepath.Join(t.TempDir(), "nf-run")
	if err := os.MkdirAll(runRoot, 0o755); err != nil {
		t.Fatalf("mkdir run root: %v", err)
	}

	got := IndexPath(domain.RunDir{Path: runRoot})
	want := filepath.Join(runRoot, IndexDirName, IndexFileName)
	if got != want {
		t.Fatalf("IndexPath() = %q, want %q", got, want)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestDiscoverArtifactsWithResultsDirSelectsNewestTraceAcrossRootAndPipelineInfoAndPairsLog(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "run")
	resultsRoot := filepath.Join(workspace, "external-results")
	pipelineInfoRoot := filepath.Join(resultsRoot, PipelineInfoDirName)
	if err := os.MkdirAll(runRoot, 0o755); err != nil {
		t.Fatalf("mkdir run root: %v", err)
	}
	if err := os.MkdirAll(pipelineInfoRoot, 0o755); err != nil {
		t.Fatalf("mkdir pipeline_info root: %v", err)
	}

	base := time.Unix(1_700_000_000, 0)
	writeSourceFixture(t, runRoot, "trace-root.txt", []byte("root older\n"), base.Add(time.Second))
	selectedTrace := writeSourceFixture(t, pipelineInfoRoot, "execution_trace_2026-02-23.txt", []byte("pipeline newer\n"), base.Add(4*time.Second))
	writeSourceFixture(t, pipelineInfoRoot, "timeline.html", []byte("ignored timeline\n"), base.Add(5*time.Second))
	writeSourceFixture(t, runRoot, ".nextflow_1.log", []byte("old log\n"), base.Add(2*time.Second))
	selectedLog := writeSourceFixture(t, runRoot, ".nextflow.log", []byte("new log\n"), base.Add(3*time.Second))

	before := time.Now().UTC()
	got, err := DiscoverArtifactsWithResultsDir(context.Background(), domain.RunDir{Path: runRoot}, domain.ResultsDir{Path: resultsRoot})
	after := time.Now().UTC()
	if err != nil {
		t.Fatalf("DiscoverArtifactsWithResultsDir(pipeline trace) returned error: %v", err)
	}

	if got.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("mode = %q, want %q", got.Mode, domain.IndexModeTraceBacked)
	}
	assertResultsDiscoverySnapshot(t, got, runRoot, resultsRoot, before, after)
	assertArtifactSource(t, got.Trace, domain.SourceKindTrace, selectedTrace, base.Add(4*time.Second), int64(len("pipeline newer\n")))
	assertArtifactSource(t, got.Log, domain.SourceKindLog, selectedLog, base.Add(3*time.Second), int64(len("new log\n")))
	if len(got.Diagnostics) != 0 {
		t.Fatalf("diagnostics = %+v, want none", got.Diagnostics)
	}
	assertNoGoshDirectory(t, runRoot)
	assertNoGoshDirectory(t, resultsRoot)
}

func TestDiscoverArtifactsWithResultsDirBreaksTraceTiesByPathAcrossLocations(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "run")
	resultsRoot := filepath.Join(workspace, "results")
	pipelineInfoRoot := filepath.Join(resultsRoot, PipelineInfoDirName)
	if err := os.MkdirAll(runRoot, 0o755); err != nil {
		t.Fatalf("mkdir run root: %v", err)
	}
	if err := os.MkdirAll(pipelineInfoRoot, 0o755); err != nil {
		t.Fatalf("mkdir pipeline_info root: %v", err)
	}

	tied := time.Unix(1_700_000_000, 0)
	writeSourceFixture(t, runRoot, "trace-root.txt", []byte("root\n"), tied)
	selectedTrace := writeSourceFixture(t, pipelineInfoRoot, "execution_trace_2026-02-23.txt", []byte("pipeline\n"), tied)

	before := time.Now().UTC()
	got, err := DiscoverArtifactsWithResultsDir(context.Background(), domain.RunDir{Path: runRoot}, domain.ResultsDir{Path: resultsRoot})
	after := time.Now().UTC()
	if err != nil {
		t.Fatalf("DiscoverArtifactsWithResultsDir(tied traces) returned error: %v", err)
	}

	if got.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("mode = %q, want %q", got.Mode, domain.IndexModeTraceBacked)
	}
	assertResultsDiscoverySnapshot(t, got, runRoot, resultsRoot, before, after)
	assertArtifactSource(t, got.Trace, domain.SourceKindTrace, selectedTrace, tied, int64(len("pipeline\n")))
	if got.Log != nil {
		t.Fatalf("log = %+v, want nil when no supported log exists", got.Log)
	}
	assertNoGoshDirectory(t, runRoot)
	assertNoGoshDirectory(t, resultsRoot)
}

func TestDiscoverArtifactsWithResultsDirUsesNewestLogOnlyWhenNoTraceExistsAnywhere(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "run")
	resultsRoot := filepath.Join(workspace, "results")
	pipelineInfoRoot := filepath.Join(resultsRoot, PipelineInfoDirName)
	if err := os.MkdirAll(runRoot, 0o755); err != nil {
		t.Fatalf("mkdir run root: %v", err)
	}
	if err := os.MkdirAll(pipelineInfoRoot, 0o755); err != nil {
		t.Fatalf("mkdir pipeline_info root: %v", err)
	}

	base := time.Unix(1_700_000_000, 0)
	writeSourceFixture(t, runRoot, "pipeline.log", []byte("unsupported log name\n"), base.Add(4*time.Second))
	writeSourceFixture(t, pipelineInfoRoot, "trace-not-matched.txt", []byte("ignored trace name\n"), base.Add(5*time.Second))
	writeSourceFixture(t, runRoot, ".nextflow_1.log", []byte("old log\n"), base.Add(time.Second))
	selectedLog := writeSourceFixture(t, runRoot, ".nextflow_2.log", []byte("new log\n"), base.Add(2*time.Second))

	before := time.Now().UTC()
	got, err := DiscoverArtifactsWithResultsDir(context.Background(), domain.RunDir{Path: runRoot}, domain.ResultsDir{Path: resultsRoot})
	after := time.Now().UTC()
	if err != nil {
		t.Fatalf("DiscoverArtifactsWithResultsDir(log-only) returned error: %v", err)
	}

	if got.Mode != domain.IndexModeLogOnly {
		t.Fatalf("mode = %q, want %q", got.Mode, domain.IndexModeLogOnly)
	}
	assertResultsDiscoverySnapshot(t, got, runRoot, resultsRoot, before, after)
	if got.Trace != nil {
		t.Fatalf("trace = %+v, want nil in log-only mode", got.Trace)
	}
	assertArtifactSource(t, got.Log, domain.SourceKindLog, selectedLog, base.Add(2*time.Second), int64(len("new log\n")))
	if len(got.Diagnostics) != 0 {
		t.Fatalf("diagnostics = %+v, want none", got.Diagnostics)
	}
	assertNoGoshDirectory(t, runRoot)
	assertNoGoshDirectory(t, resultsRoot)
}

func TestDiscoverArtifactsUsesDefaultResultsPipelineInfoTrace(t *testing.T) {
	workspace := t.TempDir()
	cwd := filepath.Join(workspace, "cwd")
	runRoot := filepath.Join(workspace, "runs", "nf-run")
	defaultResultsRoot := filepath.Join(runRoot, DefaultResultsDirName)
	pipelineInfoRoot := filepath.Join(defaultResultsRoot, PipelineInfoDirName)
	cwdPipelineInfoRoot := filepath.Join(cwd, DefaultResultsDirName, PipelineInfoDirName)
	if err := os.MkdirAll(pipelineInfoRoot, 0o755); err != nil {
		t.Fatalf("mkdir default pipeline_info root: %v", err)
	}
	if err := os.MkdirAll(cwdPipelineInfoRoot, 0o755); err != nil {
		t.Fatalf("mkdir cwd pipeline_info decoy root: %v", err)
	}
	t.Chdir(cwd)

	base := time.Unix(1_700_000_000, 0)
	selectedTrace := writeSourceFixture(t, pipelineInfoRoot, "execution_trace_2026-02-23.txt", []byte("pipeline trace\n"), base.Add(2*time.Second))
	selectedLog := writeSourceFixture(t, runRoot, ".nextflow.log", []byte("run log\n"), base.Add(time.Second))
	writeSourceFixture(t, cwdPipelineInfoRoot, "execution_trace_decoy.txt", []byte("cwd decoy\n"), base.Add(time.Hour))

	before := time.Now().UTC()
	got, err := DiscoverArtifacts(context.Background(), domain.RunDir{Path: runRoot})
	after := time.Now().UTC()
	if err != nil {
		t.Fatalf("DiscoverArtifacts(default pipeline-info trace) returned error: %v", err)
	}

	if got.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("mode = %q, want %q", got.Mode, domain.IndexModeTraceBacked)
	}
	assertDiscoverySnapshot(t, got, runRoot, before, after)
	assertArtifactSource(t, got.Trace, domain.SourceKindTrace, selectedTrace, base.Add(2*time.Second), int64(len("pipeline trace\n")))
	assertArtifactSource(t, got.Log, domain.SourceKindLog, selectedLog, base.Add(time.Second), int64(len("run log\n")))
	if len(got.Diagnostics) != 0 {
		t.Fatalf("diagnostics = %+v, want none", got.Diagnostics)
	}
	assertNoGoshDirectory(t, runRoot)
	assertNoGoshDirectory(t, defaultResultsRoot)
	assertNoGoshDirectory(t, cwd)
}

func TestDiscoverArtifactsWithNilContextReturnsDefaultDiscoverySnapshot(t *testing.T) {
	runRoot := t.TempDir()

	before := time.Now().UTC()
	got, err := DiscoverArtifacts(nil, domain.RunDir{Path: runRoot})
	after := time.Now().UTC()
	if err == nil {
		t.Fatalf("DiscoverArtifacts(nil context) returned nil error")
	}
	if err.Error() != "discover artifacts: nil context" {
		t.Fatalf("error = %q, want %q", err.Error(), "discover artifacts: nil context")
	}
	assertDiscoverySnapshot(t, got, runRoot, before, after)
	if got.Mode != "" {
		t.Fatalf("mode = %q, want empty on nil context", got.Mode)
	}
	if got.Trace != nil || got.Log != nil {
		t.Fatalf("trace/log = %+v/%+v, want nil on nil context", got.Trace, got.Log)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestDiscoverArtifactsSelectsNewestTraceAndPairsNewestLog(t *testing.T) {
	workspace := t.TempDir()
	cwd := filepath.Join(workspace, "cwd")
	runRoot := filepath.Join(workspace, "runs", "nf-run")
	if err := os.MkdirAll(cwd, 0o755); err != nil {
		t.Fatalf("mkdir cwd fixture: %v", err)
	}
	if err := os.MkdirAll(runRoot, 0o755); err != nil {
		t.Fatalf("mkdir run root: %v", err)
	}
	t.Chdir(cwd)

	base := time.Unix(1_700_000_000, 0)
	writeSourceFixture(t, runRoot, "trace-old.txt", []byte("old\n"), base)
	newestTrace := writeSourceFixture(t, runRoot, "trace-new.csv", []byte("newer\n"), base.Add(4*time.Second))
	writeSourceFixture(t, runRoot, ".nextflow_1.log", []byte("old log\n"), base.Add(time.Second))
	newestLog := writeSourceFixture(t, runRoot, ".nextflow.log", []byte("new log\n"), base.Add(3*time.Second))
	writeSourceFixture(t, cwd, "trace-decoy.txt", []byte("decoy trace\n"), base.Add(time.Hour))
	writeSourceFixture(t, cwd, ".nextflow.log", []byte("decoy log\n"), base.Add(time.Hour))

	before := time.Now().UTC()
	got, err := DiscoverArtifacts(context.Background(), domain.RunDir{Path: runRoot})
	after := time.Now().UTC()
	if err != nil {
		t.Fatalf("DiscoverArtifacts(trace-backed) returned error: %v", err)
	}

	if got.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("mode = %q, want %q", got.Mode, domain.IndexModeTraceBacked)
	}
	assertDiscoverySnapshot(t, got, runRoot, before, after)
	assertArtifactSource(t, got.Trace, domain.SourceKindTrace, newestTrace, base.Add(4*time.Second), int64(len("newer\n")))
	assertArtifactSource(t, got.Log, domain.SourceKindLog, newestLog, base.Add(3*time.Second), int64(len("new log\n")))
	if len(got.Diagnostics) != 0 {
		t.Fatalf("diagnostics = %+v, want none", got.Diagnostics)
	}
	assertNoGoshDirectory(t, runRoot)
	assertNoGoshDirectory(t, cwd)
}

func TestDiscoverArtifactsUsesNewestLogWhenTraceMissing(t *testing.T) {
	runRoot := t.TempDir()
	base := time.Unix(1_700_000_000, 0)
	writeSourceFixture(t, runRoot, ".nextflow_1.log", []byte("old log\n"), base)
	newestLog := writeSourceFixture(t, runRoot, ".nextflow_2.log", []byte("new log\n"), base.Add(2*time.Second))
	writeSourceFixture(t, runRoot, "pipeline.log", []byte("not a supported log\n"), base.Add(time.Hour))

	before := time.Now().UTC()
	got, err := DiscoverArtifacts(context.Background(), domain.RunDir{Path: runRoot})
	after := time.Now().UTC()
	if err != nil {
		t.Fatalf("DiscoverArtifacts(log-only) returned error: %v", err)
	}

	if got.Mode != domain.IndexModeLogOnly {
		t.Fatalf("mode = %q, want %q", got.Mode, domain.IndexModeLogOnly)
	}
	assertDiscoverySnapshot(t, got, runRoot, before, after)
	if got.Trace != nil {
		t.Fatalf("trace = %+v, want nil in log-only mode", got.Trace)
	}
	assertArtifactSource(t, got.Log, domain.SourceKindLog, newestLog, base.Add(2*time.Second), int64(len("new log\n")))
	if len(got.Diagnostics) != 0 {
		t.Fatalf("diagnostics = %+v, want none", got.Diagnostics)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestDiscoverArtifactsBreaksNewestTiesByPath(t *testing.T) {
	runRoot := t.TempDir()
	tied := time.Unix(1_700_000_000, 0)
	writeSourceFixture(t, runRoot, "trace-b.txt", []byte("b\n"), tied)
	traceA := writeSourceFixture(t, runRoot, "trace-a.txt", []byte("a\n"), tied)
	writeSourceFixture(t, runRoot, ".nextflow_2.log", []byte("two\n"), tied)
	logOne := writeSourceFixture(t, runRoot, ".nextflow_1.log", []byte("one\n"), tied)

	before := time.Now().UTC()
	got, err := DiscoverArtifacts(context.Background(), domain.RunDir{Path: runRoot})
	after := time.Now().UTC()
	if err != nil {
		t.Fatalf("DiscoverArtifacts(tied sources) returned error: %v", err)
	}

	if got.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("mode = %q, want %q", got.Mode, domain.IndexModeTraceBacked)
	}
	assertDiscoverySnapshot(t, got, runRoot, before, after)
	assertArtifactSource(t, got.Trace, domain.SourceKindTrace, traceA, tied, int64(len("a\n")))
	assertArtifactSource(t, got.Log, domain.SourceKindLog, logOne, tied, int64(len("one\n")))
	assertNoGoshDirectory(t, runRoot)
}

func TestDiscoverArtifactsReportsUnsupportedWithDiagnostics(t *testing.T) {
	runRoot := t.TempDir()
	writeSourceFixture(t, runRoot, "pipeline.log", []byte("not a supported artifact\n"), time.Unix(1_700_000_000, 0))

	before := time.Now().UTC()
	got, err := DiscoverArtifacts(context.Background(), domain.RunDir{Path: runRoot})
	after := time.Now().UTC()
	if err != nil {
		t.Fatalf("DiscoverArtifacts(unsupported) returned error: %v", err)
	}

	if got.Mode != domain.IndexModeUnsupported {
		t.Fatalf("mode = %q, want %q", got.Mode, domain.IndexModeUnsupported)
	}
	assertDiscoverySnapshot(t, got, runRoot, before, after)
	if got.Trace != nil {
		t.Fatalf("trace = %+v, want nil in unsupported mode", got.Trace)
	}
	if got.Log != nil {
		t.Fatalf("log = %+v, want nil in unsupported mode", got.Log)
	}
	wantDiagnostics := UnsupportedArtifactDiagnostics(domain.RunDir{Path: runRoot})
	assertDiagnostics(t, got.Diagnostics, wantDiagnostics)
	assertNoGoshDirectory(t, runRoot)
}

func TestSourceFingerprintForPathStatsCleanAbsolutePath(t *testing.T) {
	workspace := t.TempDir()
	runRoot := filepath.Join(workspace, "runs", "nf-run")
	if err := os.MkdirAll(runRoot, 0o755); err != nil {
		t.Fatalf("mkdir run root: %v", err)
	}
	filePath := filepath.Join(runRoot, "trace.txt")
	content := []byte("task_id\tstatus\n1\tFAILED\n")
	if err := os.WriteFile(filePath, content, 0o644); err != nil {
		t.Fatalf("write trace fixture: %v", err)
	}
	modTime := time.Unix(1_700_000_000, 0)
	if err := os.Chtimes(filePath, modTime, modTime); err != nil {
		t.Fatalf("set fixture mtime: %v", err)
	}
	t.Chdir(workspace)

	got, err := SourceFingerprintForPath(domain.SourceKindTrace, filepath.Join("runs", "nf-run", "..", "nf-run", "trace.txt"))
	if err != nil {
		t.Fatalf("SourceFingerprintForPath(trace) returned error: %v", err)
	}

	wantPath, err := filepath.Abs(filePath)
	if err != nil {
		t.Fatalf("filepath.Abs(%q): %v", filePath, err)
	}
	if got.Kind != domain.SourceKindTrace {
		t.Fatalf("kind = %q, want %q", got.Kind, domain.SourceKindTrace)
	}
	if got.Path != filepath.Clean(wantPath) {
		t.Fatalf("path = %q, want %q", got.Path, filepath.Clean(wantPath))
	}
	if got.Size != int64(len(content)) {
		t.Fatalf("size = %d, want %d", got.Size, len(content))
	}
	if !got.ModTime.Equal(modTime) {
		t.Fatalf("mod time = %s, want %s", got.ModTime, modTime)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestSourceFingerprintForPathUsesStatOnly(t *testing.T) {
	filePath := filepath.Join(t.TempDir(), ".nextflow.log")
	content := []byte("workflow failed\n")
	if err := os.WriteFile(filePath, content, 0o644); err != nil {
		t.Fatalf("write log fixture: %v", err)
	}
	if err := os.Chmod(filePath, 0o000); err != nil {
		t.Fatalf("remove fixture read permissions: %v", err)
	}
	t.Cleanup(func() { _ = os.Chmod(filePath, 0o644) })

	got, err := SourceFingerprintForPath(domain.SourceKindLog, filePath)
	if err != nil {
		t.Fatalf("SourceFingerprintForPath(log without read permission) returned error: %v", err)
	}
	if got.Kind != domain.SourceKindLog {
		t.Fatalf("kind = %q, want %q", got.Kind, domain.SourceKindLog)
	}
	if got.Size != int64(len(content)) {
		t.Fatalf("size = %d, want %d", got.Size, len(content))
	}
}

func TestSourceFingerprintForPathRejectsMissingSource(t *testing.T) {
	missing := filepath.Join(t.TempDir(), "missing-trace.txt")

	got, err := SourceFingerprintForPath(domain.SourceKindTrace, missing)
	if err == nil {
		t.Fatalf("SourceFingerprintForPath(%q) returned nil error", missing)
	}
	if got != (domain.SourceFingerprint{}) {
		t.Fatalf("fingerprint on error = %+v, want zero value", got)
	}
}

func TestSourceFingerprintForPathRejectsDirectorySource(t *testing.T) {
	dirPath := t.TempDir()

	got, err := SourceFingerprintForPath(domain.SourceKindLog, dirPath)
	if err == nil {
		t.Fatalf("SourceFingerprintForPath(directory %q) returned nil error", dirPath)
	}
	if !strings.Contains(err.Error(), "not a regular file") {
		t.Fatalf("error = %q, want it to mention not a regular file", err.Error())
	}
	if got != (domain.SourceFingerprint{}) {
		t.Fatalf("fingerprint on error = %+v, want zero value", got)
	}
}

func TestFindCandidateSourcesReturnsSortedDirectTraceFingerprints(t *testing.T) {
	runRoot := t.TempDir()
	base := time.Unix(1_700_000_000, 0)

	traceZ := writeSourceFixture(t, runRoot, "trace-z.txt", []byte("z\n"), base.Add(3*time.Second))
	traceA := writeSourceFixture(t, runRoot, "trace-a.csv", []byte("a,b\n"), base.Add(time.Second))
	traceM := writeSourceFixture(t, runRoot, "trace-m.tsv", []byte("a\tb\n"), base.Add(2*time.Second))
	writeSourceFixture(t, runRoot, ".nextflow.log", []byte("not a trace\n"), base.Add(4*time.Second))
	nestedDir := filepath.Join(runRoot, "nested")
	if err := os.MkdirAll(nestedDir, 0o755); err != nil {
		t.Fatalf("mkdir nested fixture dir: %v", err)
	}
	writeSourceFixture(t, nestedDir, "trace-nested.txt", []byte("nested\n"), base.Add(5*time.Second))

	got, err := FindCandidateSources(domain.RunDir{Path: runRoot}, domain.SourceKindTrace, TracePatterns)
	if err != nil {
		t.Fatalf("FindCandidateSources(trace) returned error: %v", err)
	}

	want := []domain.SourceFingerprint{
		{Kind: domain.SourceKindTrace, Path: traceA, ModTime: base.Add(time.Second), Size: int64(len("a,b\n"))},
		{Kind: domain.SourceKindTrace, Path: traceM, ModTime: base.Add(2 * time.Second), Size: int64(len("a\tb\n"))},
		{Kind: domain.SourceKindTrace, Path: traceZ, ModTime: base.Add(3 * time.Second), Size: int64(len("z\n"))},
	}
	assertSourceFingerprints(t, got, want)
	assertNoGoshDirectory(t, runRoot)
}

func TestFindCandidateSourcesReturnsSortedLogFingerprints(t *testing.T) {
	runRoot := t.TempDir()
	base := time.Unix(1_700_000_000, 0)

	currentLog := writeSourceFixture(t, runRoot, ".nextflow.log", []byte("current\n"), base.Add(3*time.Second))
	rotatedTwo := writeSourceFixture(t, runRoot, ".nextflow_2.log", []byte("two\n"), base.Add(time.Second))
	rotatedOne := writeSourceFixture(t, runRoot, ".nextflow_1.log", []byte("one\n"), base.Add(2*time.Second))
	writeSourceFixture(t, runRoot, "trace.txt", []byte("not a log\n"), base.Add(4*time.Second))

	got, err := FindCandidateSources(domain.RunDir{Path: runRoot}, domain.SourceKindLog, LogPatterns)
	if err != nil {
		t.Fatalf("FindCandidateSources(log) returned error: %v", err)
	}

	want := []domain.SourceFingerprint{
		{Kind: domain.SourceKindLog, Path: currentLog, ModTime: base.Add(3 * time.Second), Size: int64(len("current\n"))},
		{Kind: domain.SourceKindLog, Path: rotatedOne, ModTime: base.Add(2 * time.Second), Size: int64(len("one\n"))},
		{Kind: domain.SourceKindLog, Path: rotatedTwo, ModTime: base.Add(time.Second), Size: int64(len("two\n"))},
	}
	assertSourceFingerprints(t, got, want)
	assertNoGoshDirectory(t, runRoot)
}

func TestFindCandidateSourcesReturnsEmptyWhenNoPatternsMatch(t *testing.T) {
	runRoot := t.TempDir()
	writeSourceFixture(t, runRoot, "not-a-trace.log", []byte("ignored\n"), time.Unix(1_700_000_000, 0))

	got, err := FindCandidateSources(domain.RunDir{Path: runRoot}, domain.SourceKindTrace, TracePatterns)
	if err != nil {
		t.Fatalf("FindCandidateSources(no matches) returned error: %v", err)
	}
	if len(got) != 0 {
		t.Fatalf("candidate count = %d, want 0: %+v", len(got), got)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestFindCandidateSourcesIgnoresNonRegularMatches(t *testing.T) {
	runRoot := t.TempDir()
	base := time.Unix(1_700_000_000, 0)
	regularTrace := writeSourceFixture(t, runRoot, "trace-ok.txt", []byte("ok\n"), base)
	matchingDir := filepath.Join(runRoot, "trace-dir.txt")
	if err := os.Mkdir(matchingDir, 0o755); err != nil {
		t.Fatalf("mkdir matching directory fixture: %v", err)
	}

	got, err := FindCandidateSources(domain.RunDir{Path: runRoot}, domain.SourceKindTrace, TracePatterns)
	if err != nil {
		t.Fatalf("FindCandidateSources(with matching directory) returned error: %v", err)
	}

	want := []domain.SourceFingerprint{
		{Kind: domain.SourceKindTrace, Path: regularTrace, ModTime: base, Size: int64(len("ok\n"))},
	}
	assertSourceFingerprints(t, got, want)
	assertNoGoshDirectory(t, runRoot)
}

func TestFindCandidateSourcesRejectsBadPattern(t *testing.T) {
	runRoot := t.TempDir()

	got, err := FindCandidateSources(domain.RunDir{Path: runRoot}, domain.SourceKindTrace, []string{"["})
	if err == nil {
		t.Fatalf("FindCandidateSources(bad pattern) returned nil error")
	}
	if !strings.Contains(err.Error(), "glob source pattern") {
		t.Fatalf("error = %q, want it to mention glob source pattern", err.Error())
	}
	if len(got) != 0 {
		t.Fatalf("candidates on error = %+v, want none", got)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestChooseNewestSourceRejectsEmptyInput(t *testing.T) {
	got, err := ChooseNewestSource(nil)
	if err == nil {
		t.Fatalf("ChooseNewestSource(nil) returned nil error")
	}
	if !strings.Contains(err.Error(), "no candidate sources") {
		t.Fatalf("error = %q, want it to mention no candidate sources", err.Error())
	}
	if got != nil {
		t.Fatalf("source on error = %+v, want nil", got)
	}
}

func TestChooseNewestSourceSelectsLatestModTime(t *testing.T) {
	base := time.Unix(1_700_000_000, 0)
	sources := []domain.SourceFingerprint{
		{Kind: domain.SourceKindTrace, Path: "/run/trace-old.txt", ModTime: base, Size: 100},
		{Kind: domain.SourceKindTrace, Path: "/run/trace-new.txt", ModTime: base.Add(2 * time.Minute), Size: 200},
		{Kind: domain.SourceKindTrace, Path: "/run/trace-middle.txt", ModTime: base.Add(time.Minute), Size: 150},
	}

	got, err := ChooseNewestSource(sources)
	if err != nil {
		t.Fatalf("ChooseNewestSource(sources) returned error: %v", err)
	}
	if got == nil {
		t.Fatalf("ChooseNewestSource(sources) returned nil source")
	}
	if got.Path != "/run/trace-new.txt" {
		t.Fatalf("selected path = %q, want %q", got.Path, "/run/trace-new.txt")
	}
	if !got.ModTime.Equal(base.Add(2 * time.Minute)) {
		t.Fatalf("selected mtime = %s, want %s", got.ModTime, base.Add(2*time.Minute))
	}
}

func TestChooseNewestSourceBreaksMTimeTieByPath(t *testing.T) {
	tied := time.Unix(1_700_000_000, 0)
	sources := []domain.SourceFingerprint{
		{Kind: domain.SourceKindLog, Path: "/run/.nextflow_2.log", ModTime: tied, Size: 200},
		{Kind: domain.SourceKindLog, Path: "/run/.nextflow.log", ModTime: tied, Size: 100},
		{Kind: domain.SourceKindLog, Path: "/run/.nextflow_1.log", ModTime: tied, Size: 150},
	}

	got, err := ChooseNewestSource(sources)
	if err != nil {
		t.Fatalf("ChooseNewestSource(tied sources) returned error: %v", err)
	}
	if got == nil {
		t.Fatalf("ChooseNewestSource(tied sources) returned nil source")
	}
	// Equal mtimes are intentionally resolved by ascending path order so the
	// result does not depend on filesystem glob order or caller slice order.
	if got.Path != "/run/.nextflow.log" {
		t.Fatalf("selected path = %q, want stable smallest path %q", got.Path, "/run/.nextflow.log")
	}
}

func TestChooseNewestSourceTieSelectionIsIndependentOfInputOrder(t *testing.T) {
	tied := time.Unix(1_700_000_000, 0)
	forward := []domain.SourceFingerprint{
		{Kind: domain.SourceKindTrace, Path: "/run/trace-b.txt", ModTime: tied, Size: 2},
		{Kind: domain.SourceKindTrace, Path: "/run/trace-a.txt", ModTime: tied, Size: 1},
	}
	reversed := []domain.SourceFingerprint{
		{Kind: domain.SourceKindTrace, Path: "/run/trace-a.txt", ModTime: tied, Size: 1},
		{Kind: domain.SourceKindTrace, Path: "/run/trace-b.txt", ModTime: tied, Size: 2},
	}

	gotForward, err := ChooseNewestSource(forward)
	if err != nil {
		t.Fatalf("ChooseNewestSource(forward) returned error: %v", err)
	}
	gotReversed, err := ChooseNewestSource(reversed)
	if err != nil {
		t.Fatalf("ChooseNewestSource(reversed) returned error: %v", err)
	}
	if gotForward == nil || gotReversed == nil {
		t.Fatalf("ChooseNewestSource returned nil source for tied inputs: forward=%+v reversed=%+v", gotForward, gotReversed)
	}
	if gotForward.Path != "/run/trace-a.txt" || gotReversed.Path != "/run/trace-a.txt" {
		t.Fatalf("tie selection paths = %q and %q, want both %q", gotForward.Path, gotReversed.Path, "/run/trace-a.txt")
	}
}

func TestUnsupportedArtifactDiagnosticsReportsUnsupportedRunDirAndSearchedPatterns(t *testing.T) {
	runRoot := t.TempDir()

	got := UnsupportedArtifactDiagnostics(domain.RunDir{Path: runRoot})

	if len(got) != 2 {
		t.Fatalf("diagnostic count = %d, want 2: %+v", len(got), got)
	}
	missing := got[0]
	if missing.Severity != domain.DiagnosticError {
		t.Fatalf("missing artifact severity = %q, want %q", missing.Severity, domain.DiagnosticError)
	}
	if missing.Code != "unsupported_artifacts" {
		t.Fatalf("missing artifact code = %q, want %q", missing.Code, "unsupported_artifacts")
	}
	wantMessage := "No supported Nextflow trace or log artifacts found in " + runRoot
	if missing.Message != wantMessage {
		t.Fatalf("missing artifact message = %q, want %q", missing.Message, wantMessage)
	}
	wantDetail := "Searched trace locations:\n" +
		"- " + runRoot + " (run directory trace files): " + strings.Join(TracePatterns, ", ") + "\n" +
		"- " + filepath.Join(runRoot, DefaultResultsDirName, PipelineInfoDirName) + " (pipeline_info execution trace files): " + strings.Join(PipelineInfoTracePatterns, ", ") + "\n" +
		"Searched log locations:\n" +
		"- " + runRoot + " (run directory log files): " + strings.Join(LogPatterns, ", ")
	if missing.Detail != wantDetail {
		t.Fatalf("missing artifact detail = %q, want %q", missing.Detail, wantDetail)
	}
	assertNoGoshDirectory(t, runRoot)
}

func TestUnsupportedArtifactDiagnosticsDoNotClaimLogOnlyCompleteness(t *testing.T) {
	got := UnsupportedArtifactDiagnostics(domain.RunDir{Path: "/tmp/nf-run"})

	var combined strings.Builder
	for _, diagnostic := range got {
		combined.WriteString(diagnostic.Message)
		combined.WriteByte('\n')
		combined.WriteString(diagnostic.Detail)
		combined.WriteByte('\n')
	}
	for _, forbidden := range []string{
		"log-only completeness",
		"complete log-only",
		"log-only complete",
		"logs can build a complete",
		"complete task/resource/status data from logs",
	} {
		if strings.Contains(combined.String(), forbidden) {
			t.Fatalf("unsupported diagnostics contain unsupported log-only completeness claim %q in %q", forbidden, combined.String())
		}
	}
}

func TestUnsupportedArtifactDiagnosticsIncludesWithTraceRecommendation(t *testing.T) {
	got := UnsupportedArtifactDiagnostics(domain.RunDir{Path: "/tmp/nf-run"})

	if len(got) != 2 {
		t.Fatalf("diagnostic count = %d, want 2: %+v", len(got), got)
	}
	recommendation := got[1]
	if recommendation.Severity != domain.DiagnosticInfo {
		t.Fatalf("recommendation severity = %q, want %q", recommendation.Severity, domain.DiagnosticInfo)
	}
	if recommendation.Code != "nextflow_with_trace_recommended" {
		t.Fatalf("recommendation code = %q, want %q", recommendation.Code, "nextflow_with_trace_recommended")
	}
	wantMessage := "Run future Nextflow workflows with -with-trace"
	if recommendation.Message != wantMessage {
		t.Fatalf("recommendation message = %q, want %q", recommendation.Message, wantMessage)
	}
	wantDetail := "Use `nextflow run ... -with-trace` for future runs so gosh can build a complete trace-backed task index."
	if recommendation.Detail != wantDetail {
		t.Fatalf("recommendation detail = %q, want %q", recommendation.Detail, wantDetail)
	}
}

func TestUnsupportedArtifactDiagnosticsAreDeterministic(t *testing.T) {
	runDir := domain.RunDir{Path: "/tmp/nf-run"}
	first := UnsupportedArtifactDiagnostics(runDir)
	second := UnsupportedArtifactDiagnostics(runDir)

	if len(first) != len(second) {
		t.Fatalf("diagnostic counts differ: first=%d second=%d", len(first), len(second))
	}
	for i := range first {
		if first[i] != second[i] {
			t.Fatalf("diagnostic %d differs between calls: first=%+v second=%+v", i, first[i], second[i])
		}
	}
}

func assertDiscoverySnapshot(t *testing.T, got domain.ArtifactSet, wantRunRoot string, before time.Time, after time.Time) {
	t.Helper()
	assertResultsDiscoverySnapshot(t, got, wantRunRoot, filepath.Join(wantRunRoot, DefaultResultsDirName), before, after)
}

func assertSearchedPatterns(t *testing.T, got []string) {
	t.Helper()
	want := append([]string{}, TracePatterns...)
	want = append(want, LogPatterns...)
	if len(got) != len(want) {
		t.Fatalf("searched patterns count = %d, want %d: got=%+v want=%+v", len(got), len(want), got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("searched pattern %d = %q, want %q", i, got[i], want[i])
		}
	}
}

func assertResultsDiscoverySnapshot(t *testing.T, got domain.ArtifactSet, wantRunRoot string, wantResultsRoot string, before time.Time, after time.Time) {
	t.Helper()
	if got.RunDir.Path != wantRunRoot {
		t.Fatalf("run dir = %q, want %q", got.RunDir.Path, wantRunRoot)
	}
	if got.ResultsDir.Path != wantResultsRoot {
		t.Fatalf("results dir = %q, want %q", got.ResultsDir.Path, wantResultsRoot)
	}
	if got.SelectedAt.IsZero() {
		t.Fatalf("selected_at is zero")
	}
	if got.SelectedAt.Before(before) || got.SelectedAt.After(after) {
		t.Fatalf("selected_at = %s, want between %s and %s", got.SelectedAt, before, after)
	}
	assertResultsSearchedPatterns(t, got.SearchedPatterns)
	if len(got.SearchLocations) != 3 {
		t.Fatalf("search locations count = %d, want 3: %+v", len(got.SearchLocations), got.SearchLocations)
	}
	assertArtifactSearchLocation(t, got.SearchLocations[0], domain.SourceKindTrace, wantRunRoot, TracePatterns, "run directory", "trace")
	assertArtifactSearchLocation(t, got.SearchLocations[1], domain.SourceKindTrace, filepath.Join(wantResultsRoot, PipelineInfoDirName), PipelineInfoTracePatterns, "pipeline_info", "execution trace")
	assertArtifactSearchLocation(t, got.SearchLocations[2], domain.SourceKindLog, wantRunRoot, LogPatterns, "run directory", "log")
}

func assertResultsSearchedPatterns(t *testing.T, got []string) {
	t.Helper()
	want := append([]string{}, TracePatterns...)
	want = append(want, PipelineInfoTracePatterns...)
	want = append(want, LogPatterns...)
	if len(got) != len(want) {
		t.Fatalf("searched patterns count = %d, want %d: got=%+v want=%+v", len(got), len(want), got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("searched pattern %d = %q, want %q", i, got[i], want[i])
		}
	}
}

func assertArtifactSearchLocation(t *testing.T, got domain.ArtifactSearchLocation, wantKind domain.SourceKind, wantBaseDir string, wantPatterns []string, wantDescriptionSubstrings ...string) {
	t.Helper()
	if got.Kind != wantKind {
		t.Fatalf("location kind = %q, want %q", got.Kind, wantKind)
	}
	if got.BaseDir != wantBaseDir {
		t.Fatalf("location base dir = %q, want %q", got.BaseDir, wantBaseDir)
	}
	assertPatterns(t, got.Patterns, wantPatterns)
	if got.Description == "" {
		t.Fatalf("location description is empty")
	}
	for _, substring := range wantDescriptionSubstrings {
		if !strings.Contains(got.Description, substring) {
			t.Fatalf("location description = %q, want it to contain %q", got.Description, substring)
		}
	}
}

func assertPatterns(t *testing.T, got []string, want []string) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("pattern count = %d, want %d: got=%+v want=%+v", len(got), len(want), got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("pattern %d = %q, want %q", i, got[i], want[i])
		}
	}
}

func assertArtifactSource(t *testing.T, got *domain.SourceFingerprint, wantKind domain.SourceKind, wantPath string, wantModTime time.Time, wantSize int64) {
	t.Helper()
	if got == nil {
		t.Fatalf("source = nil, want %q source at %q", wantKind, wantPath)
	}
	if got.Kind != wantKind {
		t.Fatalf("source kind = %q, want %q", got.Kind, wantKind)
	}
	if got.Path != wantPath {
		t.Fatalf("source path = %q, want %q", got.Path, wantPath)
	}
	if !got.ModTime.Equal(wantModTime) {
		t.Fatalf("source mod time = %s, want %s", got.ModTime, wantModTime)
	}
	if got.Size != wantSize {
		t.Fatalf("source size = %d, want %d", got.Size, wantSize)
	}
}

func assertDiagnostics(t *testing.T, got []domain.Diagnostic, want []domain.Diagnostic) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("diagnostic count = %d, want %d: got=%+v want=%+v", len(got), len(want), got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("diagnostic %d = %+v, want %+v", i, got[i], want[i])
		}
	}
}

func writeSourceFixture(t *testing.T, dir string, name string, content []byte, modTime time.Time) string {
	t.Helper()
	path := filepath.Join(dir, name)
	if err := os.WriteFile(path, content, 0o644); err != nil {
		t.Fatalf("write source fixture %q: %v", path, err)
	}
	if err := os.Chtimes(path, modTime, modTime); err != nil {
		t.Fatalf("set source fixture mtime %q: %v", path, err)
	}
	return filepath.Clean(path)
}

func assertSourceFingerprints(t *testing.T, got []domain.SourceFingerprint, want []domain.SourceFingerprint) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("candidate count = %d, want %d: got=%+v want=%+v", len(got), len(want), got, want)
	}
	for i := range want {
		if got[i].Kind != want[i].Kind {
			t.Fatalf("candidate %d kind = %q, want %q", i, got[i].Kind, want[i].Kind)
		}
		if got[i].Path != want[i].Path {
			t.Fatalf("candidate %d path = %q, want %q", i, got[i].Path, want[i].Path)
		}
		if got[i].Size != want[i].Size {
			t.Fatalf("candidate %d size = %d, want %d", i, got[i].Size, want[i].Size)
		}
		if !got[i].ModTime.Equal(want[i].ModTime) {
			t.Fatalf("candidate %d mod time = %s, want %s", i, got[i].ModTime, want[i].ModTime)
		}
	}
}

func assertNoGoshDirectory(t *testing.T, runRoot string) {
	t.Helper()
	goshPath := filepath.Join(runRoot, IndexDirName)
	if _, err := os.Stat(goshPath); err == nil {
		t.Fatalf("ResolveRunDir mutated run root by creating %q", goshPath)
	} else if !os.IsNotExist(err) {
		t.Fatalf("stat %q: %v", goshPath, err)
	}
}
