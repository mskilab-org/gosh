package tasks

import (
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
	"time"

	"github.com/mskilab-org/gosh/internal/domain"
)

func TestNormalizeTaskQueryTrimsTextFiltersWithoutChangingCase(t *testing.T) {
	query := domain.TaskQuery{
		ProcessSubstring: "  Align_STAR\t",
		NameSubstring:    "\nSample-Tumor  ",
		SampleSubstring:  "\tSAMPLE_T1 ",
	}

	got, err := NormalizeTaskQuery(query)
	if err != nil {
		t.Fatalf("NormalizeTaskQuery() returned error: %v", err)
	}

	want := domain.TaskQuery{
		ProcessSubstring: "Align_STAR",
		NameSubstring:    "Sample-Tumor",
		SampleSubstring:  "SAMPLE_T1",
	}
	if got != want {
		t.Fatalf("NormalizeTaskQuery() = %#v, want %#v", got, want)
	}
}

func TestNormalizeTaskQueryNormalizesStatusRaw(t *testing.T) {
	query := domain.TaskQuery{StatusRaw: "  failed\t"}

	got, err := NormalizeTaskQuery(query)
	if err != nil {
		t.Fatalf("NormalizeTaskQuery() returned error: %v", err)
	}

	want := domain.TaskQuery{Status: domain.TaskStatusFailed, StatusRaw: "failed"}
	if got != want {
		t.Fatalf("NormalizeTaskQuery() = %#v, want %#v", got, want)
	}
}

func TestNormalizeTaskQueryDoesNotTreatBlankStatusRawAsUnknownFilter(t *testing.T) {
	query := domain.TaskQuery{StatusRaw: " \t\n ", ProcessSubstring: " align "}

	got, err := NormalizeTaskQuery(query)
	if err != nil {
		t.Fatalf("NormalizeTaskQuery() returned error: %v", err)
	}

	want := domain.TaskQuery{ProcessSubstring: "align"}
	if got != want {
		t.Fatalf("NormalizeTaskQuery() = %#v, want %#v", got, want)
	}
}

func TestNormalizeTaskQueryRejectsUnknownStatusRaw(t *testing.T) {
	_, err := NormalizeTaskQuery(domain.TaskQuery{StatusRaw: "finished"})

	assertErrorContains(t, err, "unknown task status")
	assertErrorContains(t, err, "finished")
}

func TestNormalizeTaskQueryNormalizesExistingStatusWhenStatusRawAbsent(t *testing.T) {
	query := domain.TaskQuery{Status: domain.TaskStatus("completed")}

	got, err := NormalizeTaskQuery(query)
	if err != nil {
		t.Fatalf("NormalizeTaskQuery() returned error: %v", err)
	}

	want := domain.TaskQuery{Status: domain.TaskStatusCompleted}
	if got != want {
		t.Fatalf("NormalizeTaskQuery() = %#v, want %#v", got, want)
	}
}

func TestApplyTaskQueryReturnsAllTasksSortedBySourceOrderWithoutMutatingInput(t *testing.T) {
	taskList := []domain.Task{
		{RowOrder: 30, ID: "cc/333333", Status: domain.TaskStatusCompleted, Process: "LATE"},
		{RowOrder: 10, ID: "aa/111111", Status: domain.TaskStatusFailed, Process: "EARLY"},
		{RowOrder: 20, ID: "bb/222222", Status: domain.TaskStatusCached, Process: "MIDDLE"},
	}

	got, err := ApplyTaskQuery(taskList, domain.TaskQuery{})
	if err != nil {
		t.Fatalf("ApplyTaskQuery() returned error: %v", err)
	}

	want := []domain.Task{taskList[1], taskList[2], taskList[0]}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("ApplyTaskQuery() = %#v, want %#v", got, want)
	}

	gotInputIDs := []string{taskList[0].ID, taskList[1].ID, taskList[2].ID}
	wantInputIDs := []string{"cc/333333", "aa/111111", "bb/222222"}
	if !reflect.DeepEqual(gotInputIDs, wantInputIDs) {
		t.Fatalf("input task IDs after ApplyTaskQuery() = %#v, want %#v", gotInputIDs, wantInputIDs)
	}
}

func TestApplyTaskQueryAppliesCaseInsensitiveProcessNameAndStatusFilters(t *testing.T) {
	taskList := []domain.Task{
		{RowOrder: 30, ID: "cc/333333", Status: domain.TaskStatusFailed, Process: "CALL_VARIANTS", Name: "sample-tumor"},
		{RowOrder: 10, ID: "aa/111111", Status: domain.TaskStatusFailed, Process: "ALIGN_STAR", Name: "Sample-Tumor"},
		{RowOrder: 20, ID: "bb/222222", Status: domain.TaskStatusCompleted, Process: "align_star", Name: "sample-tumor"},
		{RowOrder: 40, ID: "dd/444444", Status: domain.TaskStatusFailed, Process: "ALIGN_STAR", Name: "control"},
	}
	query := domain.TaskQuery{
		ProcessSubstring: " align ",
		NameSubstring:    "TUMOR",
		StatusRaw:        " failed ",
	}

	got, err := ApplyTaskQuery(taskList, query)
	if err != nil {
		t.Fatalf("ApplyTaskQuery() returned error: %v", err)
	}

	want := []domain.Task{taskList[1]}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("ApplyTaskQuery() = %#v, want %#v", got, want)
	}
}

func TestApplyTaskQueryAppliesSampleAliasToNameOrTagOnly(t *testing.T) {
	taskList := []domain.Task{
		{RowOrder: 30, ID: "cc/333333", Process: "ALIGN_STAR", Name: "Sample-Tumor", Tag: "lane-a"},
		{RowOrder: 10, ID: "aa/111111", Process: "ALIGN_STAR", Name: "sample-control", Tag: "Tumor-Replicate"},
		{RowOrder: 20, ID: "bb/222222", Process: "TUMOR_PROCESS", Name: "control", Tag: "normal"},
	}

	got, err := ApplyTaskQuery(taskList, domain.TaskQuery{SampleSubstring: " tumor "})
	if err != nil {
		t.Fatalf("ApplyTaskQuery() returned error: %v", err)
	}

	want := []domain.Task{taskList[1], taskList[0]}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("ApplyTaskQuery() = %#v, want %#v", got, want)
	}
}

func TestApplyTaskQueryReturnsEmptySliceWhenNoRowsMatch(t *testing.T) {
	taskList := []domain.Task{{RowOrder: 1, ID: "aa/111111", Process: "ALIGN_STAR"}}

	got, err := ApplyTaskQuery(taskList, domain.TaskQuery{ProcessSubstring: "CALL"})
	if err != nil {
		t.Fatalf("ApplyTaskQuery() returned error: %v", err)
	}

	if got == nil || len(got) != 0 {
		t.Fatalf("ApplyTaskQuery() = %#v, want empty non-nil slice", got)
	}
}

func TestApplyTaskQueryRejectsUnknownStatusFilter(t *testing.T) {
	_, err := ApplyTaskQuery(nil, domain.TaskQuery{StatusRaw: "finished"})

	assertErrorContains(t, err, "unknown task status")
	assertErrorContains(t, err, "finished")
}

func TestMatchCanonicalIDMatchesExactCanonicalID(t *testing.T) {
	task := domain.Task{ID: "ab/c123def"}

	if !MatchCanonicalID("ab/c123def", task) {
		t.Fatalf("MatchCanonicalID() = false, want true for exact canonical ID")
	}
}

func TestMatchCanonicalIDNormalizesCaseAndSpace(t *testing.T) {
	task := domain.Task{ID: "ab/c123def"}

	if !MatchCanonicalID("  AB/C123DEF\t", task) {
		t.Fatalf("MatchCanonicalID() = false, want true after canonical ID casing normalization")
	}
}

func TestMatchCanonicalIDRequiresWholeCanonicalID(t *testing.T) {
	task := domain.Task{ID: "ab/c123def"}
	selectors := []string{
		"ab/c123",
		"c123def",
		"xxab/c123def",
		"ab/c123def00",
		"ab/c123def/.command.err",
	}

	for _, selector := range selectors {
		t.Run(selector, func(t *testing.T) {
			if MatchCanonicalID(selector, task) {
				t.Fatalf("MatchCanonicalID(%q) = true, want false for non-exact selector", selector)
			}
		})
	}
}

func TestMatchCanonicalIDRejectsNonCanonicalOrEmptySelectors(t *testing.T) {
	task := domain.Task{ID: "ab/c123def"}
	selectors := []string{
		"",
		" \t\n ",
		"abc123def",
		"/runs/example/work/ab/c123def",
	}

	for _, selector := range selectors {
		t.Run(selector, func(t *testing.T) {
			if MatchCanonicalID(selector, task) {
				t.Fatalf("MatchCanonicalID(%q) = true, want false for non-canonical selector", selector)
			}
		})
	}

	if MatchCanonicalID("ab/c123def", domain.Task{}) {
		t.Fatalf("MatchCanonicalID() = true, want false when task has no canonical ID")
	}
}

func TestMatchWorkdirPathMatchesCleanedFullPath(t *testing.T) {
	workdir := filepath.Join(t.TempDir(), "work", "ab", "c123def")
	separator := string(os.PathSeparator)
	selector := "  " + filepath.Dir(workdir) + separator + "." + separator + filepath.Base(workdir) + separator + "\t"
	task := domain.Task{Workdir: workdir}

	if !MatchWorkdirPath(selector, task) {
		t.Fatalf("MatchWorkdirPath() = false, want true for cleaned full workdir path")
	}
}

func TestMatchWorkdirPathMatchesRelativeSelectorAfterAbsoluteResolution(t *testing.T) {
	root := t.TempDir()
	t.Chdir(root)

	selector := filepath.Join("work", "de", "f456")
	workdir, err := filepath.Abs(selector)
	if err != nil {
		t.Fatalf("filepath.Abs(%q): %v", selector, err)
	}
	task := domain.Task{Workdir: workdir}

	if !MatchWorkdirPath(selector, task) {
		t.Fatalf("MatchWorkdirPath() = false, want true when relative selector resolves to task workdir")
	}
}

func TestMatchWorkdirPathRequiresWholePath(t *testing.T) {
	workdir := filepath.Join(t.TempDir(), "work", "ab", "c123def")
	task := domain.Task{Workdir: workdir}
	selectors := []string{
		filepath.Join("ab", "c123def"),
		filepath.Base(workdir),
		filepath.Dir(workdir),
		workdir + "-extra",
		filepath.Join(workdir, ".command.err"),
	}

	for _, selector := range selectors {
		t.Run(selector, func(t *testing.T) {
			if MatchWorkdirPath(selector, task) {
				t.Fatalf("MatchWorkdirPath(%q) = true, want false for non-exact workdir selector", selector)
			}
		})
	}
}

func TestMatchWorkdirPathRejectsEmptySelectorOrUnknownWorkdir(t *testing.T) {
	workdir := filepath.Join(t.TempDir(), "work", "ab", "c123def")
	task := domain.Task{Workdir: workdir}

	for _, selector := range []string{"", " \t\n "} {
		t.Run("selector", func(t *testing.T) {
			if MatchWorkdirPath(selector, task) {
				t.Fatalf("MatchWorkdirPath(%q) = true, want false for empty selector", selector)
			}
		})
	}

	for _, unknownWorkdir := range []string{"", " \t\n "} {
		t.Run("workdir", func(t *testing.T) {
			if MatchWorkdirPath(workdir, domain.Task{Workdir: unknownWorkdir}) {
				t.Fatalf("MatchWorkdirPath() = true, want false for unknown task workdir %q", unknownWorkdir)
			}
		})
	}
}

func TestMatchHumanSelectorMatchesProcessNameOrTagCaseInsensitively(t *testing.T) {
	task := domain.Task{
		Process: "NFCORE_RNA:ALIGN_STAR",
		Name:    "sample-Tumor-A",
		Tag:     "SAMPLE_T1 replicate",
	}

	selectors := []string{
		"align_star", // process substring
		"TUMOR",      // name substring with different case
		"sample_t1",  // tag substring with different case
	}

	for _, selector := range selectors {
		t.Run(selector, func(t *testing.T) {
			if !MatchHumanSelector(selector, task) {
				t.Fatalf("MatchHumanSelector(%q) = false, want true", selector)
			}
		})
	}
}

func TestMatchHumanSelectorTrimsSelectorWhitespace(t *testing.T) {
	task := domain.Task{Process: "ALIGN_STAR"}

	if !MatchHumanSelector("  align_star\t", task) {
		t.Fatalf("MatchHumanSelector() = false, want true after trimming selector whitespace")
	}
}

func TestMatchHumanSelectorRejectsEmptySelectors(t *testing.T) {
	task := domain.Task{Process: "ALIGN_STAR", Name: "sample", Tag: "tumor"}

	for _, selector := range []string{"", " \t\n "} {
		t.Run(selector, func(t *testing.T) {
			if MatchHumanSelector(selector, task) {
				t.Fatalf("MatchHumanSelector(%q) = true, want false for empty selector", selector)
			}
		})
	}
}

func TestMatchHumanSelectorDoesNotMatchIDOrWorkdir(t *testing.T) {
	task := domain.Task{
		ID:      "ab/c123def",
		Process: "ALIGN_STAR",
		Name:    "control",
		Tag:     "normal",
		Workdir: "/runs/sample/work/ab/c123def",
	}

	selectors := []string{
		"ab/c123def",
		"sample/work",
	}

	for _, selector := range selectors {
		t.Run(selector, func(t *testing.T) {
			if MatchHumanSelector(selector, task) {
				t.Fatalf("MatchHumanSelector(%q) = true, want false for ID/workdir-only selector", selector)
			}
		})
	}
}

func TestMatchHumanSelectorRejectsAbsentProcessNameTagSubstring(t *testing.T) {
	task := domain.Task{Process: "ALIGN_STAR", Name: "control", Tag: "normal"}

	if MatchHumanSelector("variant", task) {
		t.Fatalf("MatchHumanSelector() = true, want false when process/name/tag do not contain selector")
	}
}

func TestResolveSelectorPrefersCanonicalIDBeforeHumanMatches(t *testing.T) {
	taskList := []domain.Task{
		{RowOrder: 1, ID: "aa/111111", Process: "ALIGN_STAR", Name: "tumor", Workdir: filepath.FromSlash("/runs/example/work/aa/111111")},
		{RowOrder: 2, ID: "bb/222222", Process: "qa aa/111111 review", Name: "control", Workdir: filepath.FromSlash("/runs/example/work/bb/222222")},
	}

	got, err := ResolveSelector("  AA/111111\t", taskList)
	if err != nil {
		t.Fatalf("ResolveSelector() returned error: %v", err)
	}

	assertExactResolution(t, got, "AA/111111", taskList[0], "matched canonical id")
}

func TestResolveSelectorPrefersWorkdirPathBeforeHumanMatches(t *testing.T) {
	root := t.TempDir()
	workdir := filepath.Join(root, "work", "aa", "111111")
	selector := filepath.Join(root, "work", "aa", ".", "111111")
	taskList := []domain.Task{
		{RowOrder: 1, ID: "aa/111111", Process: "ALIGN_STAR", Name: "tumor", Workdir: workdir},
		{RowOrder: 2, ID: "bb/222222", Process: "mentions " + selector, Name: "control", Workdir: filepath.Join(root, "work", "bb", "222222")},
	}

	got, err := ResolveSelector(selector, taskList)
	if err != nil {
		t.Fatalf("ResolveSelector() returned error: %v", err)
	}

	assertExactResolution(t, got, selector, taskList[0], "matched full workdir path")
}

func TestResolveSelectorReturnsExactForSingleHumanMatch(t *testing.T) {
	taskList := []domain.Task{
		{RowOrder: 1, ID: "aa/111111", Process: "ALIGN_STAR", Name: "control", Tag: "normal", Workdir: filepath.FromSlash("/runs/example/work/aa/111111")},
		{RowOrder: 2, ID: "bb/222222", Process: "CALL_VARIANTS", Name: "tumor", Tag: "somatic", Workdir: filepath.FromSlash("/runs/example/work/bb/222222")},
	}

	got, err := ResolveSelector("  TUMOR ", taskList)
	if err != nil {
		t.Fatalf("ResolveSelector() returned error: %v", err)
	}

	assertExactResolution(t, got, "TUMOR", taskList[1], "matched process/name/tag")
}

func TestResolveSelectorReturnsAmbiguousForMultipleHumanMatches(t *testing.T) {
	taskList := []domain.Task{
		{RowOrder: 1, ID: "aa/111111", Process: "ALIGN_STAR", Name: "tumor", Workdir: filepath.FromSlash("/runs/example/work/aa/111111")},
		{RowOrder: 2, ID: "bb/222222", Process: "ALIGN_BWA", Name: "control", Workdir: filepath.FromSlash("/runs/example/work/bb/222222")},
		{RowOrder: 3, ID: "cc/333333", Process: "CALL_VARIANTS", Name: "tumor", Workdir: filepath.FromSlash("/runs/example/work/cc/333333")},
	}

	got, err := ResolveSelector("align", taskList)
	if err != nil {
		t.Fatalf("ResolveSelector() returned error: %v", err)
	}

	if got.Kind != domain.SelectorResolutionAmbiguous {
		t.Fatalf("ResolveSelector() kind = %q, want %q", got.Kind, domain.SelectorResolutionAmbiguous)
	}
	if got.Selector != "align" {
		t.Fatalf("ResolveSelector() selector = %q, want %q", got.Selector, "align")
	}
	if got.Task != nil {
		t.Fatalf("ResolveSelector() task = %#v, want nil for ambiguous selector", got.Task)
	}
	wantMatches := []domain.Task{taskList[0], taskList[1]}
	if !reflect.DeepEqual(got.Matches, wantMatches) {
		t.Fatalf("ResolveSelector() matches = %#v, want %#v", got.Matches, wantMatches)
	}
	assertResolutionDiagnostic(t, got.Diagnostics, domain.DiagnosticWarning, "selector_ambiguous", "selector matched more than one task")
}

func TestResolveSelectorReturnsNotFoundForMissingOrBlankSelector(t *testing.T) {
	taskList := []domain.Task{{RowOrder: 1, ID: "aa/111111", Process: "ALIGN_STAR", Name: "tumor", Workdir: filepath.FromSlash("/runs/example/work/aa/111111")}}

	for _, selector := range []string{"missing", " \t\n "} {
		t.Run(selector, func(t *testing.T) {
			got, err := ResolveSelector(selector, taskList)
			if err != nil {
				t.Fatalf("ResolveSelector() returned error: %v", err)
			}

			if got.Kind != domain.SelectorResolutionNotFound {
				t.Fatalf("ResolveSelector() kind = %q, want %q", got.Kind, domain.SelectorResolutionNotFound)
			}
			if got.Selector != strings.TrimSpace(selector) {
				t.Fatalf("ResolveSelector() selector = %q, want %q", got.Selector, strings.TrimSpace(selector))
			}
			if got.Task != nil {
				t.Fatalf("ResolveSelector() task = %#v, want nil for not-found selector", got.Task)
			}
			if got.Matches == nil || len(got.Matches) != 0 {
				t.Fatalf("ResolveSelector() matches = %#v, want empty non-nil slice", got.Matches)
			}
			assertResolutionDiagnostic(t, got.Diagnostics, domain.DiagnosticError, "selector_not_found", "selector did not match any indexed task")
		})
	}
}

func TestBuildStatusSummaryCombinesMetadataCountsAndSources(t *testing.T) {
	builtAt := time.Date(2024, time.April, 28, 12, 34, 56, 0, time.UTC)
	traceMod := time.Date(2024, time.April, 28, 12, 0, 0, 0, time.UTC)
	logMod := time.Date(2024, time.April, 28, 12, 1, 0, 0, time.UTC)
	runDir := domain.RunDir{Path: "/runs/example"}
	metadata := domain.IndexMetadata{
		RunDir:    runDir.Path,
		IndexPath: filepath.FromSlash("/runs/example/.gosh/index.sqlite"),
		Mode:      domain.IndexModeTraceBacked,
		Trace:     &domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: filepath.FromSlash("/runs/example/trace.txt"), ModTime: traceMod, Size: 1200},
		Log:       &domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: filepath.FromSlash("/runs/example/.nextflow.log"), ModTime: logMod, Size: 3400},
		BuiltAt:   builtAt,
		Freshness: domain.IndexFreshnessFresh,
	}
	counts := []domain.StatusCount{
		{Status: domain.TaskStatusCompleted, Count: 7},
		{Status: domain.TaskStatusFailed, Count: 2},
	}

	got, err := BuildStatusSummary(runDir, metadata, counts, nil)
	if err != nil {
		t.Fatalf("BuildStatusSummary() returned error: %v", err)
	}

	if got.RunDir != runDir {
		t.Fatalf("RunDir = %#v, want %#v", got.RunDir, runDir)
	}
	if got.Mode != metadata.Mode {
		t.Fatalf("Mode = %q, want %q", got.Mode, metadata.Mode)
	}
	if got.IndexPath != metadata.IndexPath {
		t.Fatalf("IndexPath = %q, want %q", got.IndexPath, metadata.IndexPath)
	}
	if got.Freshness != metadata.Freshness {
		t.Fatalf("Freshness = %q, want %q", got.Freshness, metadata.Freshness)
	}
	if got.BuiltAt == nil || !got.BuiltAt.Equal(builtAt) {
		t.Fatalf("BuiltAt = %#v, want %s", got.BuiltAt, builtAt.Format(time.RFC3339Nano))
	}
	wantSources := domain.ArtifactSet{
		RunDir: runDir,
		Mode:   metadata.Mode,
		Trace:  metadata.Trace,
		Log:    metadata.Log,
	}
	if !reflect.DeepEqual(got.Sources, wantSources) {
		t.Fatalf("Sources = %#v, want %#v", got.Sources, wantSources)
	}
	if !reflect.DeepEqual(got.Counts, counts) {
		t.Fatalf("Counts = %#v, want %#v", got.Counts, counts)
	}
	if got.FailedCount != 0 {
		t.Fatalf("FailedCount = %d, want 0 with no provided failed tasks", got.FailedCount)
	}
	if len(got.FailedPreview) != 0 {
		t.Fatalf("FailedPreview = %#v, want empty", got.FailedPreview)
	}
	if len(got.LogOnlyFailures) != 0 {
		t.Fatalf("LogOnlyFailures = %#v, want empty for trace-backed summary", got.LogOnlyFailures)
	}
}

func TestBuildStatusSummaryUsesFailedTasksForFailedCountAndLimitsPreviewToThree(t *testing.T) {
	exitOne := 1
	exitTwo := 2
	exitThree := 3
	runDir := domain.RunDir{Path: "/runs/example"}
	metadata := domain.IndexMetadata{
		RunDir:    runDir.Path,
		IndexPath: filepath.FromSlash("/runs/example/.gosh/index.sqlite"),
		Mode:      domain.IndexModeTraceBacked,
		Freshness: domain.IndexFreshnessFresh,
	}
	counts := []domain.StatusCount{
		{Status: domain.TaskStatusFailed, Count: 99},
		{Status: domain.TaskStatusCompleted, Count: 1},
	}
	failedTasks := []domain.Task{
		{RowOrder: 40, ID: "dd/444444", Status: domain.TaskStatusFailed, Process: "FOURTH", Name: "sample-4", Exit: nil, ErrorSummary: "fourth failure"},
		{RowOrder: 10, ID: "aa/111111", Status: domain.TaskStatusFailed, Process: "FIRST", Name: "sample-1", Tag: "tumor", Workdir: filepath.FromSlash("/runs/example/work/aa/111111"), Exit: &exitOne, ErrorSummary: "first failure"},
		{RowOrder: 30, ID: "cc/333333", Status: domain.TaskStatusFailed, Process: "THIRD", Name: "sample-3", Exit: &exitThree, ErrorSummary: "third failure"},
		{RowOrder: 20, ID: "bb/222222", Status: domain.TaskStatusAborted, Process: "SECOND", Name: "sample-2", Exit: &exitTwo, ErrorSummary: "second abort"},
	}

	got, err := BuildStatusSummary(runDir, metadata, counts, failedTasks)
	if err != nil {
		t.Fatalf("BuildStatusSummary() returned error: %v", err)
	}

	if got.FailedCount != len(failedTasks) {
		t.Fatalf("FailedCount = %d, want %d from provided failed/error-like tasks, not status counts", got.FailedCount, len(failedTasks))
	}
	wantPreview := []domain.FailedTaskPreview{
		{ID: "aa/111111", Status: domain.TaskStatusFailed, Process: "FIRST", Name: "sample-1", Tag: "tumor", Workdir: filepath.FromSlash("/runs/example/work/aa/111111"), Exit: &exitOne, ErrorSummary: "first failure"},
		{ID: "bb/222222", Status: domain.TaskStatusAborted, Process: "SECOND", Name: "sample-2", Exit: &exitTwo, ErrorSummary: "second abort"},
		{ID: "cc/333333", Status: domain.TaskStatusFailed, Process: "THIRD", Name: "sample-3", Exit: &exitThree, ErrorSummary: "third failure"},
	}
	if !reflect.DeepEqual(got.FailedPreview, wantPreview) {
		t.Fatalf("FailedPreview = %#v, want first three failed/error-like tasks by source order: %#v", got.FailedPreview, wantPreview)
	}
}

func TestBuildStatusSummaryDoesNotAliasMutableInputs(t *testing.T) {
	builtAt := time.Date(2024, time.April, 28, 12, 34, 56, 0, time.UTC)
	trace := &domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: filepath.FromSlash("/runs/example/trace.txt"), Size: 1200}
	metadata := domain.IndexMetadata{
		RunDir:    "/runs/example",
		IndexPath: filepath.FromSlash("/runs/example/.gosh/index.sqlite"),
		Mode:      domain.IndexModeTraceBacked,
		Trace:     trace,
		BuiltAt:   builtAt,
		Freshness: domain.IndexFreshnessFresh,
	}
	counts := []domain.StatusCount{{Status: domain.TaskStatusCompleted, Count: 7}}

	got, err := BuildStatusSummary(domain.RunDir{Path: "/runs/example"}, metadata, counts, nil)
	if err != nil {
		t.Fatalf("BuildStatusSummary() returned error: %v", err)
	}

	counts[0].Count = 99
	trace.Path = filepath.FromSlash("/runs/example/changed-trace.txt")
	metadata.BuiltAt = builtAt.Add(time.Hour)

	if got.Counts[0].Count != 7 {
		t.Fatalf("Counts[0].Count = %d after mutating input counts, want 7", got.Counts[0].Count)
	}
	if got.Sources.Trace == nil || got.Sources.Trace.Path != filepath.FromSlash("/runs/example/trace.txt") {
		t.Fatalf("Sources.Trace = %#v after mutating input trace, want original path", got.Sources.Trace)
	}
	if got.BuiltAt == nil || !got.BuiltAt.Equal(builtAt) {
		t.Fatalf("BuiltAt = %#v after mutating metadata, want original %s", got.BuiltAt, builtAt.Format(time.RFC3339Nano))
	}
}

func TestSelectFailedPreviewSortsBySourceOrderAndLimits(t *testing.T) {
	exitTwo := 2
	exitOne := 1
	failedTasks := []domain.Task{
		{
			RowOrder:     30,
			ID:           "cc/333333",
			Status:       domain.TaskStatusFailed,
			Process:      "LATE_PROCESS",
			Name:         "late-sample",
			Tag:          "late-tag",
			Workdir:      "/run/work/cc/333333",
			Exit:         &exitTwo,
			Duration:     "1h",
			Realtime:     "1h",
			CPUs:         "8",
			Memory:       "32 GB",
			ErrorSummary: "late failure",
		},
		{
			RowOrder:     10,
			ID:           "aa/111111",
			Status:       domain.TaskStatusFailed,
			Process:      "EARLY_PROCESS",
			Name:         "early-sample",
			Tag:          "early-tag",
			Workdir:      "/run/work/aa/111111",
			Exit:         &exitOne,
			ErrorSummary: "early failure",
		},
		{
			RowOrder:     20,
			ID:           "bb/222222",
			Status:       domain.TaskStatusAborted,
			Process:      "MIDDLE_PROCESS",
			Name:         "middle-sample",
			Tag:          "middle-tag",
			Workdir:      "/run/work/bb/222222",
			Exit:         nil,
			ErrorSummary: "middle abort",
		},
	}

	got := SelectFailedPreview(failedTasks, 2)
	want := []domain.FailedTaskPreview{
		{
			ID:           "aa/111111",
			Status:       domain.TaskStatusFailed,
			Process:      "EARLY_PROCESS",
			Name:         "early-sample",
			Tag:          "early-tag",
			Workdir:      "/run/work/aa/111111",
			Exit:         &exitOne,
			ErrorSummary: "early failure",
		},
		{
			ID:           "bb/222222",
			Status:       domain.TaskStatusAborted,
			Process:      "MIDDLE_PROCESS",
			Name:         "middle-sample",
			Tag:          "middle-tag",
			Workdir:      "/run/work/bb/222222",
			Exit:         nil,
			ErrorSummary: "middle abort",
		},
	}

	if !reflect.DeepEqual(got, want) {
		t.Fatalf("SelectFailedPreview() = %#v, want %#v", got, want)
	}
}

func TestSelectFailedPreviewReturnsAllWhenLimitExceedsTaskCount(t *testing.T) {
	failedTasks := []domain.Task{
		{RowOrder: 2, ID: "bb/222222", Status: domain.TaskStatusFailed, Process: "B"},
		{RowOrder: 1, ID: "aa/111111", Status: domain.TaskStatusFailed, Process: "A"},
	}

	got := SelectFailedPreview(failedTasks, 10)
	want := []domain.FailedTaskPreview{
		{ID: "aa/111111", Status: domain.TaskStatusFailed, Process: "A"},
		{ID: "bb/222222", Status: domain.TaskStatusFailed, Process: "B"},
	}

	if !reflect.DeepEqual(got, want) {
		t.Fatalf("SelectFailedPreview() = %#v, want %#v", got, want)
	}
}

func TestSelectFailedPreviewReturnsEmptyForNonPositiveLimit(t *testing.T) {
	failedTasks := []domain.Task{{RowOrder: 1, ID: "aa/111111", Status: domain.TaskStatusFailed}}

	for _, limit := range []int{0, -1} {
		got := SelectFailedPreview(failedTasks, limit)
		if len(got) != 0 {
			t.Fatalf("SelectFailedPreview(_, %d) length = %d, want 0", limit, len(got))
		}
	}
}

func TestSelectFailedPreviewDoesNotMutateInputOrder(t *testing.T) {
	failedTasks := []domain.Task{
		{RowOrder: 2, ID: "bb/222222", Status: domain.TaskStatusFailed},
		{RowOrder: 1, ID: "aa/111111", Status: domain.TaskStatusFailed},
	}

	_ = SelectFailedPreview(failedTasks, 2)

	gotIDs := []string{failedTasks[0].ID, failedTasks[1].ID}
	wantIDs := []string{"bb/222222", "aa/111111"}
	if !reflect.DeepEqual(gotIDs, wantIDs) {
		t.Fatalf("input task IDs after SelectFailedPreview() = %#v, want %#v", gotIDs, wantIDs)
	}
}

func TestBuildTaskDossierCombinesExactResolutionInventoryAndDiagnostics(t *testing.T) {
	exitCode := 137
	task := domain.Task{
		ID:           "ab/c123def",
		Status:       domain.TaskStatusFailed,
		Process:      "ALIGN_STAR",
		Name:         "tumor-sample",
		Tag:          "tumor replicate 1",
		Workdir:      filepath.FromSlash("/runs/example/work/ab/c123def"),
		Exit:         &exitCode,
		Duration:     "1h",
		Realtime:     "58m",
		CPUs:         "8",
		Memory:       "32 GB",
		ErrorSummary: "process failed: command exited with 137",
	}
	inventory := domain.CommandFileInventory{
		Workdir: task.Workdir,
		Files: []domain.CommandFile{
			{
				Kind:   domain.CommandFileShell,
				Path:   filepath.Join(task.Workdir, string(domain.CommandFileShell)),
				Exists: true,
				Size:   42,
				Snippet: &domain.Snippet{
					Path:      filepath.Join(task.Workdir, string(domain.CommandFileShell)),
					Strategy:  domain.SnippetStrategyHead,
					StartLine: 1,
					EndLine:   1,
					Content:   "echo hello",
					MaxBytes:  4096,
				},
			},
			{
				Kind:   domain.CommandFileLog,
				Path:   filepath.Join(task.Workdir, string(domain.CommandFileLog)),
				Exists: true,
				Size:   512,
			},
		},
	}
	diagnostics := []domain.Diagnostic{
		{Severity: domain.DiagnosticInfo, Code: "selector_exact", Message: "selector resolved exactly"},
	}
	resolution := domain.SelectorResolution{
		Kind:        domain.SelectorResolutionExact,
		Selector:    task.ID,
		Task:        &task,
		Diagnostics: diagnostics,
	}

	got, err := BuildTaskDossier(resolution, inventory)
	if err != nil {
		t.Fatalf("BuildTaskDossier() error = %v, want nil", err)
	}

	want := domain.TaskDossier{Task: task, Inventory: inventory, Diagnostics: diagnostics}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("BuildTaskDossier() = %#v, want %#v", got, want)
	}
}

func TestBuildTaskDossierRejectsExactResolutionWithoutTask(t *testing.T) {
	diagnostics := []domain.Diagnostic{
		{Severity: domain.DiagnosticError, Code: "selector_exact_missing_task", Message: "exact selector did not include task data"},
	}

	got, err := BuildTaskDossier(domain.SelectorResolution{
		Kind:        domain.SelectorResolutionExact,
		Selector:    "ab/c123def",
		Diagnostics: diagnostics,
	}, domain.CommandFileInventory{Workdir: filepath.FromSlash("/runs/example/work/ab/c123def")})
	if err == nil {
		t.Fatalf("BuildTaskDossier() error = nil, want error for exact resolution without task")
	}
	assertErrorContains(t, err, "exact")
	assertErrorContains(t, err, "task")
	assertEmptyDossierExceptDiagnostics(t, got, diagnostics)
}

func TestBuildTaskDossierRejectsAmbiguousResolutionWithoutGuessing(t *testing.T) {
	diagnostics := []domain.Diagnostic{
		{Severity: domain.DiagnosticWarning, Code: "selector_ambiguous", Message: "selector matched more than one task"},
	}
	resolution := domain.SelectorResolution{
		Kind:     domain.SelectorResolutionAmbiguous,
		Selector: "ALIGN",
		Matches: []domain.Task{
			{ID: "aa/111111", Process: "ALIGN_STAR", Workdir: filepath.FromSlash("/runs/example/work/aa/111111")},
			{ID: "bb/222222", Process: "ALIGN_STAR", Workdir: filepath.FromSlash("/runs/example/work/bb/222222")},
		},
		Diagnostics: diagnostics,
	}
	inventory := domain.CommandFileInventory{
		Workdir: filepath.FromSlash("/runs/example/work/aa/111111"),
		Files:   []domain.CommandFile{{Kind: domain.CommandFileLog, Exists: true}},
	}

	got, err := BuildTaskDossier(resolution, inventory)
	if err == nil {
		t.Fatalf("BuildTaskDossier() error = nil, want error for ambiguous selector")
	}
	assertErrorContains(t, err, "ambiguous")
	assertErrorContains(t, err, "ALIGN")
	assertErrorContains(t, err, "2")
	assertEmptyDossierExceptDiagnostics(t, got, diagnostics)
}

func TestBuildTaskDossierRejectsNotFoundResolution(t *testing.T) {
	diagnostics := []domain.Diagnostic{
		{Severity: domain.DiagnosticError, Code: "selector_not_found", Message: "selector did not match any indexed task"},
	}
	resolution := domain.SelectorResolution{
		Kind:        domain.SelectorResolutionNotFound,
		Selector:    "missing-task",
		Diagnostics: diagnostics,
	}

	got, err := BuildTaskDossier(resolution, domain.CommandFileInventory{Workdir: filepath.FromSlash("/runs/example/work/cc/333333")})
	if err == nil {
		t.Fatalf("BuildTaskDossier() error = nil, want error for not-found selector")
	}
	assertErrorContains(t, err, "not found")
	assertErrorContains(t, err, "missing-task")
	assertEmptyDossierExceptDiagnostics(t, got, diagnostics)
}

func TestBuildTaskDossierRejectsUnknownResolutionKind(t *testing.T) {
	diagnostics := []domain.Diagnostic{
		{Severity: domain.DiagnosticError, Code: "selector_unknown_kind", Message: "selector resolution kind was not recognized"},
	}
	resolution := domain.SelectorResolution{
		Kind:        domain.SelectorResolutionKind("maybe"),
		Selector:    "ab/c123def",
		Diagnostics: diagnostics,
	}

	got, err := BuildTaskDossier(resolution, domain.CommandFileInventory{})
	if err == nil {
		t.Fatalf("BuildTaskDossier() error = nil, want error for unknown selector resolution kind")
	}
	assertErrorContains(t, err, "unsupported")
	assertErrorContains(t, err, "maybe")
	assertEmptyDossierExceptDiagnostics(t, got, diagnostics)
}

func assertExactResolution(t *testing.T, got domain.SelectorResolution, wantSelector string, wantTask domain.Task, wantDetail string) {
	t.Helper()
	if got.Kind != domain.SelectorResolutionExact {
		t.Fatalf("ResolveSelector() kind = %q, want %q", got.Kind, domain.SelectorResolutionExact)
	}
	if got.Selector != wantSelector {
		t.Fatalf("ResolveSelector() selector = %q, want %q", got.Selector, wantSelector)
	}
	if got.Task == nil {
		t.Fatalf("ResolveSelector() task = nil, want %#v", wantTask)
	}
	if !reflect.DeepEqual(*got.Task, wantTask) {
		t.Fatalf("ResolveSelector() task = %#v, want %#v", *got.Task, wantTask)
	}
	if got.Matches == nil || len(got.Matches) != 0 {
		t.Fatalf("ResolveSelector() matches = %#v, want empty non-nil slice", got.Matches)
	}
	assertResolutionDiagnostic(t, got.Diagnostics, domain.DiagnosticInfo, "selector_exact", "selector resolved exactly")
	if got.Diagnostics[0].Detail != wantDetail {
		t.Fatalf("ResolveSelector() diagnostic detail = %q, want %q", got.Diagnostics[0].Detail, wantDetail)
	}
}

func assertResolutionDiagnostic(t *testing.T, diagnostics []domain.Diagnostic, wantSeverity domain.DiagnosticSeverity, wantCode string, wantMessage string) {
	t.Helper()
	if len(diagnostics) != 1 {
		t.Fatalf("diagnostics length = %d, want 1: %#v", len(diagnostics), diagnostics)
	}
	want := domain.Diagnostic{Severity: wantSeverity, Code: wantCode, Message: wantMessage}
	if diagnostics[0].Severity != want.Severity || diagnostics[0].Code != want.Code || diagnostics[0].Message != want.Message {
		t.Fatalf("diagnostic = %#v, want severity/code/message %#v", diagnostics[0], want)
	}
}

func assertErrorContains(t *testing.T, err error, want string) {
	t.Helper()
	if err == nil {
		t.Fatalf("error = nil, want substring %q", want)
	}
	if !strings.Contains(err.Error(), want) {
		t.Fatalf("error = %q, want substring %q", err.Error(), want)
	}
}

func assertEmptyDossierExceptDiagnostics(t *testing.T, got domain.TaskDossier, wantDiagnostics []domain.Diagnostic) {
	t.Helper()
	if !reflect.DeepEqual(got.Task, domain.Task{}) {
		t.Fatalf("dossier task = %#v, want zero task", got.Task)
	}
	if !reflect.DeepEqual(got.Inventory, domain.CommandFileInventory{}) {
		t.Fatalf("dossier inventory = %#v, want zero inventory", got.Inventory)
	}
	if !reflect.DeepEqual(got.Diagnostics, wantDiagnostics) {
		t.Fatalf("dossier diagnostics = %#v, want %#v", got.Diagnostics, wantDiagnostics)
	}
}
