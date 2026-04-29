package trace

import (
	"context"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/mskilab-org/gosh/internal/domain"
)

func TestDefaultColumnAliasSetIncludesRootTraceAliases(t *testing.T) {
	aliases := DefaultColumnAliasSet()
	tests := []struct {
		name string
		got  []string
		want string
	}{
		{name: "hash", got: aliases.Hash, want: "hash"},
		{name: "workdir", got: aliases.Workdir, want: "workdir"},
		{name: "process", got: aliases.Process, want: "process"},
		{name: "name", got: aliases.Name, want: "name"},
		{name: "tag", got: aliases.Tag, want: "tag"},
		{name: "status", got: aliases.Status, want: "status"},
		{name: "exit", got: aliases.Exit, want: "exit"},
		{name: "duration", got: aliases.Duration, want: "duration"},
		{name: "realtime", got: aliases.Realtime, want: "realtime"},
		{name: "cpus", got: aliases.CPU, want: "cpus"},
		{name: "memory", got: aliases.Memory, want: "memory"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			requireAlias(t, tt.got, tt.want)
		})
	}
}

func TestDefaultColumnAliasSetIncludesNFCoreExecutionTraceAliases(t *testing.T) {
	aliases := DefaultColumnAliasSet()
	tests := []struct {
		name string
		got  []string
		want string
	}{
		{name: "task_id", got: aliases.TaskID, want: "task_id"},
		{name: "hash", got: aliases.Hash, want: "hash"},
		{name: "native_id", got: aliases.NativeID, want: "native_id"},
		{name: "name", got: aliases.Name, want: "name"},
		{name: "status", got: aliases.Status, want: "status"},
		{name: "exit", got: aliases.Exit, want: "exit"},
		{name: "duration", got: aliases.Duration, want: "duration"},
		{name: "realtime", got: aliases.Realtime, want: "realtime"},
		{name: "%cpu", got: aliases.CPU, want: "%cpu"},
		{name: "peak_rss", got: aliases.PeakRSS, want: "peak_rss"},
		{name: "peak_vmem", got: aliases.PeakVMem, want: "peak_vmem"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			requireAlias(t, tt.got, tt.want)
		})
	}
}

func TestDefaultColumnAliasSetHasNoBlankOrDuplicateAliases(t *testing.T) {
	aliases := DefaultColumnAliasSet()
	groups := []struct {
		name    string
		aliases []string
	}{
		{name: "TaskID", aliases: aliases.TaskID},
		{name: "Hash", aliases: aliases.Hash},
		{name: "NativeID", aliases: aliases.NativeID},
		{name: "Workdir", aliases: aliases.Workdir},
		{name: "Process", aliases: aliases.Process},
		{name: "Name", aliases: aliases.Name},
		{name: "Tag", aliases: aliases.Tag},
		{name: "Status", aliases: aliases.Status},
		{name: "Exit", aliases: aliases.Exit},
		{name: "Duration", aliases: aliases.Duration},
		{name: "Realtime", aliases: aliases.Realtime},
		{name: "CPU", aliases: aliases.CPU},
		{name: "Memory", aliases: aliases.Memory},
		{name: "PeakRSS", aliases: aliases.PeakRSS},
		{name: "PeakVMem", aliases: aliases.PeakVMem},
	}

	for _, group := range groups {
		t.Run(group.name, func(t *testing.T) {
			if len(group.aliases) == 0 {
				t.Fatalf("%s aliases are empty", group.name)
			}
			seen := make(map[string]bool, len(group.aliases))
			for _, alias := range group.aliases {
				if strings.TrimSpace(alias) == "" {
					t.Fatalf("%s aliases contain blank alias: %#v", group.name, group.aliases)
				}
				if seen[alias] {
					t.Fatalf("%s aliases contain duplicate %q: %#v", group.name, alias, group.aliases)
				}
				seen[alias] = true
			}
		})
	}
}

func TestDefaultColumnAliasSetReturnsIndependentSlices(t *testing.T) {
	first := DefaultColumnAliasSet()
	if len(first.TaskID) == 0 || len(first.Hash) == 0 {
		t.Fatalf("default aliases must include task and hash aliases: %#v", first)
	}
	first.TaskID[0] = "mutated_task_id"
	first.Hash[0] = "mutated_hash"

	second := DefaultColumnAliasSet()
	requireAlias(t, second.TaskID, "task_id")
	requireAlias(t, second.Hash, "hash")
	if second.TaskID[0] == "mutated_task_id" || second.Hash[0] == "mutated_hash" {
		t.Fatalf("DefaultColumnAliasSet returned aliases sharing mutable backing storage: first=%#v second=%#v", first, second)
	}
}

func requireAlias(t *testing.T, aliases []string, want string) {
	t.Helper()

	for _, alias := range aliases {
		if alias == want {
			return
		}
	}
	t.Fatalf("aliases %#v do not include %q", aliases, want)
}

func TestNormalizeRecordColumnsMapsRootTraceAliasesAndCopiesRawColumns(t *testing.T) {
	workdir := filepath.FromSlash("/runs/example/work/AB/C123DEF")
	record := RawRecord{
		RowOrder: 11,
		Columns: map[string]string{
			"hash":        "AB/C123DEF",
			"workdir":     workdir,
			"process":     "ALIGN_STAR",
			"name":        "ALIGN_STAR (sample-1)",
			"tag":         "sample-1",
			"status":      "FAILED",
			"exit":        "137",
			"duration":    "1h 2m",
			"realtime":    "62m",
			"cpus":        "8",
			"memory":      "16 GB",
			"custom_note": "kept verbatim",
		},
	}

	got, err := NormalizeRecordColumns(record, DefaultColumnAliasSet())
	if err != nil {
		t.Fatalf("NormalizeRecordColumns() returned error: %v", err)
	}

	want := NormalizedRecord{
		RowOrder:   11,
		Hash:       "AB/C123DEF",
		Workdir:    workdir,
		Process:    "ALIGN_STAR",
		Name:       "ALIGN_STAR (sample-1)",
		Tag:        "sample-1",
		Status:     "FAILED",
		Exit:       "137",
		Duration:   "1h 2m",
		Realtime:   "62m",
		CPUDisplay: "8",
		Memory:     "16 GB",
		Columns:    record.Columns,
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("NormalizeRecordColumns() = %#v, want %#v", got, want)
	}

	got.Columns["custom_note"] = "changed"
	if record.Columns["custom_note"] != "kept verbatim" {
		t.Fatalf("NormalizeRecordColumns() reused raw column map; record.Columns = %#v", record.Columns)
	}
}

func TestNormalizeRecordColumnsMapsExecutionTraceAliases(t *testing.T) {
	record := RawRecord{
		RowOrder: 12,
		Columns: map[string]string{
			"task_id":   "42",
			"hash":      "de/f456",
			"native_id": "12345",
			"name":      "NFGOS:AMBER_STEP:BAM_AMBER:AMBER (WG-26-03_vs_WG-26-04)",
			"status":    "CACHED",
			"exit":      "0",
			"duration":  "3m",
			"realtime":  "2m",
			"%cpu":      "87.5%",
			"peak_rss":  "4.5 GB",
			"peak_vmem": "10 GB",
			"rchar":     "1000",
		},
	}

	got, err := NormalizeRecordColumns(record, DefaultColumnAliasSet())
	if err != nil {
		t.Fatalf("NormalizeRecordColumns() returned error: %v", err)
	}

	want := NormalizedRecord{
		RowOrder:   12,
		TaskID:     "42",
		Hash:       "de/f456",
		NativeID:   "12345",
		Name:       "NFGOS:AMBER_STEP:BAM_AMBER:AMBER (WG-26-03_vs_WG-26-04)",
		Status:     "CACHED",
		Exit:       "0",
		Duration:   "3m",
		Realtime:   "2m",
		CPUDisplay: "87.5%",
		PeakRSS:    "4.5 GB",
		PeakVMem:   "10 GB",
		Columns:    record.Columns,
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("NormalizeRecordColumns() = %#v, want %#v", got, want)
	}
}

func TestNormalizeRecordColumnsSupportsCustomAliasesAndDefaultsMissingFields(t *testing.T) {
	record := RawRecord{
		RowOrder: 5,
		Columns: map[string]string{
			"task hash":     "12/345",
			"state":         "COMPLETED",
			"operator note": "manual review",
		},
	}
	aliases := ColumnAliasSet{
		Hash:   []string{"task hash"},
		Status: []string{"state"},
	}

	got, err := NormalizeRecordColumns(record, aliases)
	if err != nil {
		t.Fatalf("NormalizeRecordColumns() returned error: %v", err)
	}

	want := NormalizedRecord{
		RowOrder: 5,
		Hash:     "12/345",
		Status:   "COMPLETED",
		Columns:  record.Columns,
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("NormalizeRecordColumns() = %#v, want %#v", got, want)
	}
}

func TestNormalizeRecordColumnsHandlesEmptyColumns(t *testing.T) {
	got, err := NormalizeRecordColumns(RawRecord{RowOrder: 6, Columns: map[string]string{}}, DefaultColumnAliasSet())
	if err != nil {
		t.Fatalf("NormalizeRecordColumns() returned error: %v", err)
	}

	want := NormalizedRecord{RowOrder: 6, Columns: map[string]string{}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("NormalizeRecordColumns() = %#v, want %#v", got, want)
	}
}

func TestNormalizeRecordColumnsRejectsAmbiguousAliases(t *testing.T) {
	got, err := NormalizeRecordColumns(
		RawRecord{RowOrder: 7, Columns: map[string]string{"id": "ab/c123def"}},
		ColumnAliasSet{
			TaskID: []string{"id"},
			Hash:   []string{"id"},
		},
	)
	if err == nil {
		t.Fatalf("NormalizeRecordColumns() returned nil error and record %#v", got)
	}
	if !reflect.DeepEqual(got, NormalizedRecord{}) {
		t.Fatalf("record on error = %#v, want zero value", got)
	}
	for _, want := range []string{"normalize record columns", "alias", "id"} {
		if !strings.Contains(err.Error(), want) {
			t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
		}
	}
}

func TestDeriveNamePartsSplitsScopedProcessAndTag(t *testing.T) {
	fullName := "NFGOS:AMBER_STEP:BAM_AMBER:AMBER (WG-26-03_vs_WG-26-04)"
	process := "NFGOS:AMBER_STEP:BAM_AMBER:AMBER"
	tag := "WG-26-03_vs_WG-26-04"

	got := DeriveNameParts(fullName)
	want := NameParts{
		FullName:      fullName,
		Process:       process,
		Tag:           tag,
		SelectorTerms: []string{fullName, process, tag},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("DeriveNameParts() = %#v, want %#v", got, want)
	}
}

func TestDeriveNamePartsTrimsFullNameAndKeepsSpacedTags(t *testing.T) {
	got := DeriveNameParts("  ALIGN_STAR (tumor replicate 1)\n")
	want := NameParts{
		FullName:      "ALIGN_STAR (tumor replicate 1)",
		Process:       "ALIGN_STAR",
		Tag:           "tumor replicate 1",
		SelectorTerms: []string{"ALIGN_STAR (tumor replicate 1)", "ALIGN_STAR", "tumor replicate 1"},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("DeriveNameParts() = %#v, want %#v", got, want)
	}
}

func TestDeriveNamePartsDefaultsUntaggedNamesToProcessOnly(t *testing.T) {
	got := DeriveNameParts("PIPE:QC")
	want := NameParts{
		FullName:      "PIPE:QC",
		Process:       "PIPE:QC",
		SelectorTerms: []string{"PIPE:QC"},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("DeriveNameParts() = %#v, want %#v", got, want)
	}
}

func TestDeriveNamePartsLeavesAmbiguousParenthesesUnsplit(t *testing.T) {
	tests := []string{
		"PIPE:ALIGN (sample-1",
		"PIPE:ALIGN (sample(1))",
		"PIPE:ALIGN ()",
		"PIPE:ALIGN(sample-1)",
	}

	for _, fullName := range tests {
		t.Run(fullName, func(t *testing.T) {
			got := DeriveNameParts(fullName)
			want := NameParts{
				FullName:      fullName,
				Process:       fullName,
				SelectorTerms: []string{fullName},
			}
			if !reflect.DeepEqual(got, want) {
				t.Fatalf("DeriveNameParts(%q) = %#v, want %#v", fullName, got, want)
			}
		})
	}
}

func TestDeriveNamePartsReturnsEmptyPartsForBlankName(t *testing.T) {
	got := DeriveNameParts(" \t\n ")
	if got.FullName != "" || got.Process != "" || got.Tag != "" || len(got.SelectorTerms) != 0 {
		t.Fatalf("DeriveNameParts(blank) = %#v, want empty name parts", got)
	}
}

func TestTaskFromNormalizedRecordBuildsTaskFromExecutionTraceColumns(t *testing.T) {
	runDir := t.TempDir()
	wantWorkdir := filepath.Join(runDir, "work", "de", "f456")
	if err := os.MkdirAll(wantWorkdir, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q): %v", wantWorkdir, err)
	}

	fullName := "NFGOS:AMBER_STEP:BAM_AMBER:AMBER (WG-26-03_vs_WG-26-04)"
	got, err := TaskFromNormalizedRecord(domain.RunDir{Path: runDir}, NormalizedRecord{
		RowOrder:   12,
		TaskID:     "42",
		Hash:       "DE/F456",
		NativeID:   "12345",
		Name:       fullName,
		Status:     " cached ",
		Exit:       "0",
		Duration:   "3m",
		Realtime:   "2m",
		CPUDisplay: "87.5%",
		PeakRSS:    "4.5 GB",
		PeakVMem:   "10 GB",
	})
	if err != nil {
		t.Fatalf("TaskFromNormalizedRecord() returned error: %v", err)
	}

	wantExit := 0
	want := domain.Task{
		RowOrder: 12,
		ID:       "de/f456",
		Status:   domain.TaskStatusCached,
		Process:  "NFGOS:AMBER_STEP:BAM_AMBER:AMBER",
		Name:     fullName,
		Tag:      "WG-26-03_vs_WG-26-04",
		Workdir:  wantWorkdir,
		Exit:     &wantExit,
		Duration: "3m",
		Realtime: "2m",
		CPUs:     "87.5%",
		Memory:   "peak_rss=4.5 GB; peak_vmem=10 GB",
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("TaskFromNormalizedRecord() = %#v, want %#v", got, want)
	}
}

func TestTaskFromNormalizedRecordPreservesExplicitDisplayFields(t *testing.T) {
	workdir := filepath.Join(t.TempDir(), "custom-work-root", "AB", "C123DEF")
	got, err := TaskFromNormalizedRecord(domain.RunDir{Path: t.TempDir()}, NormalizedRecord{
		RowOrder:   7,
		Hash:       "AB/C123DEF",
		Workdir:    workdir,
		Process:    "ALIGN_STAR",
		Name:       "ALIGN_STAR (sample-1)",
		Tag:        "sample-1",
		Status:     " failed ",
		Exit:       " 137 ",
		Duration:   "1h 2m",
		Realtime:   "62m",
		CPUDisplay: "8",
		Memory:     "16 GB",
		PeakRSS:    "4.5 GB",
		PeakVMem:   "10 GB",
	})
	if err != nil {
		t.Fatalf("TaskFromNormalizedRecord() returned error: %v", err)
	}

	wantExit := 137
	want := domain.Task{
		RowOrder: 7,
		ID:       "ab/c123def",
		Status:   domain.TaskStatusFailed,
		Process:  "ALIGN_STAR",
		Name:     "ALIGN_STAR (sample-1)",
		Tag:      "sample-1",
		Workdir:  filepath.Clean(workdir),
		Exit:     &wantExit,
		Duration: "1h 2m",
		Realtime: "62m",
		CPUs:     "8",
		Memory:   "16 GB",
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("TaskFromNormalizedRecord() = %#v, want %#v", got, want)
	}
}

func TestTaskFromNormalizedRecordUsesWorkdirWhenHashColumnIsMissing(t *testing.T) {
	workdir := filepath.Join(t.TempDir(), "work", "DE", "F456")
	got, err := TaskFromNormalizedRecord(domain.RunDir{Path: t.TempDir()}, NormalizedRecord{
		RowOrder: 2,
		Workdir:  workdir,
		Status:   "COMPLETED",
	})
	if err != nil {
		t.Fatalf("TaskFromNormalizedRecord() returned error: %v", err)
	}
	if got.ID != "de/f456" {
		t.Fatalf("task ID = %q, want canonical ID from workdir", got.ID)
	}
	if got.Workdir != filepath.Clean(workdir) {
		t.Fatalf("task workdir = %q, want %q", got.Workdir, filepath.Clean(workdir))
	}
	if got.Status != domain.TaskStatusCompleted {
		t.Fatalf("task status = %q, want %q", got.Status, domain.TaskStatusCompleted)
	}
	if got.Exit != nil {
		t.Fatalf("task exit = %#v, want nil for missing exit", got.Exit)
	}
}

func TestTaskFromNormalizedRecordRejectsRecordsWithoutParseableHashOrWorkdir(t *testing.T) {
	tests := []struct {
		name   string
		record NormalizedRecord
	}{
		{name: "missing id columns", record: NormalizedRecord{Status: "COMPLETED"}},
		{name: "invalid id columns", record: NormalizedRecord{TaskID: "42", Hash: "not-a-hash", Workdir: "also-not-a-path"}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tt.record.RowOrder = 9
			got, err := TaskFromNormalizedRecord(domain.RunDir{Path: t.TempDir()}, tt.record)
			if err == nil {
				t.Fatalf("TaskFromNormalizedRecord() returned nil error and task %#v", got)
			}
			if !reflect.DeepEqual(got, domain.Task{}) {
				t.Fatalf("task on error = %#v, want zero value", got)
			}
			for _, want := range []string{"task from normalized record", "row 9", "hash", "workdir"} {
				if !strings.Contains(err.Error(), want) {
					t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
				}
			}
		})
	}
}

func TestTaskFromNormalizedRecordPropagatesInvalidExitValues(t *testing.T) {
	got, err := TaskFromNormalizedRecord(domain.RunDir{Path: t.TempDir()}, NormalizedRecord{
		RowOrder: 4,
		Hash:     "ab/c123def",
		Exit:     "killed",
	})
	if err == nil {
		t.Fatalf("TaskFromNormalizedRecord() returned nil error and task %#v", got)
	}
	if !reflect.DeepEqual(got, domain.Task{}) {
		t.Fatalf("task on error = %#v, want zero value", got)
	}
	for _, want := range []string{"task from normalized record", "row 4", "exit", "killed"} {
		if !strings.Contains(err.Error(), want) {
			t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
		}
	}
}

func TestNormalizeTaskStatusMapsKnownStatuses(t *testing.T) {
	tests := []struct {
		raw  string
		want domain.TaskStatus
	}{
		{raw: "FAILED", want: domain.TaskStatusFailed},
		{raw: "COMPLETED", want: domain.TaskStatusCompleted},
		{raw: "CACHED", want: domain.TaskStatusCached},
		{raw: "ABORTED", want: domain.TaskStatusAborted},
		{raw: "SUBMITTED", want: domain.TaskStatusSubmitted},
		{raw: "RUNNING", want: domain.TaskStatusRunning},
	}

	for _, tt := range tests {
		t.Run(tt.raw, func(t *testing.T) {
			got := NormalizeTaskStatus(tt.raw)
			if got != tt.want {
				t.Fatalf("NormalizeTaskStatus(%q) = %q, want %q", tt.raw, got, tt.want)
			}
		})
	}
}

func TestNormalizeTaskStatusIgnoresCaseAndSurroundingWhitespace(t *testing.T) {
	tests := []struct {
		raw  string
		want domain.TaskStatus
	}{
		{raw: " failed ", want: domain.TaskStatusFailed},
		{raw: "\tCompleted\n", want: domain.TaskStatusCompleted},
		{raw: "CaChEd", want: domain.TaskStatusCached},
		{raw: " aborted ", want: domain.TaskStatusAborted},
		{raw: " submitted ", want: domain.TaskStatusSubmitted},
		{raw: " running ", want: domain.TaskStatusRunning},
	}

	for _, tt := range tests {
		t.Run(tt.raw, func(t *testing.T) {
			got := NormalizeTaskStatus(tt.raw)
			if got != tt.want {
				t.Fatalf("NormalizeTaskStatus(%q) = %q, want %q", tt.raw, got, tt.want)
			}
		})
	}
}

func TestNormalizeTaskStatusReturnsUnknownForBlankAndUnknownValues(t *testing.T) {
	tests := []string{"", "   ", "PENDING", "DONE", "FAILED (retry)"}

	for _, raw := range tests {
		t.Run(raw, func(t *testing.T) {
			got := NormalizeTaskStatus(raw)
			if got != domain.TaskStatusUnknown {
				t.Fatalf("NormalizeTaskStatus(%q) = %q, want %q", raw, got, domain.TaskStatusUnknown)
			}
		})
	}
}

func TestParseNullableExitReturnsNilForBlankAndMissingLikeValues(t *testing.T) {
	tests := []string{"", "   ", "\t\n", "-", "NA", "n/a", "NULL"}

	for _, raw := range tests {
		t.Run(raw, func(t *testing.T) {
			got, err := ParseNullableExit(raw)
			if err != nil {
				t.Fatalf("ParseNullableExit(%q) returned error: %v", raw, err)
			}
			if got != nil {
				t.Fatalf("ParseNullableExit(%q) = %d, want nil", raw, *got)
			}
		})
	}
}

func TestParseNullableExitParsesIntegers(t *testing.T) {
	tests := []struct {
		raw  string
		want int
	}{
		{raw: "0", want: 0},
		{raw: "1", want: 1},
		{raw: " 137\n", want: 137},
		{raw: "-1", want: -1},
	}

	for _, tt := range tests {
		t.Run(tt.raw, func(t *testing.T) {
			got, err := ParseNullableExit(tt.raw)
			if err != nil {
				t.Fatalf("ParseNullableExit(%q) returned error: %v", tt.raw, err)
			}
			if got == nil {
				t.Fatalf("ParseNullableExit(%q) = nil, want %d", tt.raw, tt.want)
			}
			if *got != tt.want {
				t.Fatalf("ParseNullableExit(%q) = %d, want %d", tt.raw, *got, tt.want)
			}
		})
	}
}

func TestParseNullableExitRejectsInvalidNonblankValues(t *testing.T) {
	tests := []string{"failed", "1.5", "137 killed"}

	for _, raw := range tests {
		t.Run(raw, func(t *testing.T) {
			got, err := ParseNullableExit(raw)
			if err == nil {
				t.Fatalf("ParseNullableExit(%q) returned nil error and value %#v", raw, got)
			}
			if got != nil {
				t.Fatalf("value on error = %d, want nil", *got)
			}
			for _, want := range []string{"parse nullable exit", raw} {
				if !strings.Contains(err.Error(), want) {
					t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
				}
			}
		})
	}
}

func TestDetectDelimiterChoosesDelimiterBySupportedExtension(t *testing.T) {
	tests := []struct {
		name string
		path string
		want Delimiter
	}{
		{name: "csv uses comma", path: "/runs/example/trace.csv", want: DelimiterComma},
		{name: "tsv uses tab", path: "/runs/example/trace.tsv", want: DelimiterTab},
		{name: "txt uses tab", path: "/runs/example/trace.txt", want: DelimiterTab},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := DetectDelimiter(tt.path)
			if err != nil {
				t.Fatalf("DetectDelimiter(%q) returned error: %v", tt.path, err)
			}
			if got != tt.want {
				t.Fatalf("DetectDelimiter(%q) = %q, want %q", tt.path, rune(got), rune(tt.want))
			}
		})
	}
}

func TestDetectDelimiterRejectsUnsupportedExtensionClearly(t *testing.T) {
	path := "/runs/example/trace.json"

	got, err := DetectDelimiter(path)
	if err == nil {
		t.Fatalf("DetectDelimiter(%q) returned nil error", path)
	}
	if got != 0 {
		t.Fatalf("delimiter on error = %q, want zero value", rune(got))
	}
	for _, want := range []string{"unsupported", ".json", ".csv", ".tsv", ".txt"} {
		if !strings.Contains(err.Error(), want) {
			t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
		}
	}
}

func TestDetectDelimiterRejectsEmptyPathClearly(t *testing.T) {
	got, err := DetectDelimiter("")
	if err == nil {
		t.Fatalf("DetectDelimiter(empty path) returned nil error")
	}
	if got != 0 {
		t.Fatalf("delimiter on error = %q, want zero value", rune(got))
	}
	if !strings.Contains(err.Error(), "empty path") {
		t.Fatalf("error = %q, want it to mention empty path", err.Error())
	}
}

func TestParseTraceParsesCSVFixtureInSourceOrder(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	firstWorkdir := filepath.Join(runDir.Path, "custom-work", "AB", "C123DEF")
	secondWorkdir := filepath.Join(runDir.Path, "custom-work", "DE", "F456")
	tracePath := writeTraceFixture(t, runDir.Path, "trace.csv", strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AB/C123DEF,COMPLETED,ALIGN,ALIGN (sample-1),sample-1," + firstWorkdir + ",0,1m,60s,2,4 GB",
		"DE/F456,FAILED,QUANT,QUANT (sample-2),sample-2," + secondWorkdir + ",137,2m,120s,4,8 GB",
		"",
	}, "\n"))

	got, err := ParseTrace(context.Background(), runDir, domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: tracePath})
	if err != nil {
		t.Fatalf("ParseTrace(csv fixture) returned error: %v", err)
	}

	firstExit := 0
	secondExit := 137
	want := []domain.Task{
		{
			RowOrder: 1,
			ID:       "ab/c123def",
			Status:   domain.TaskStatusCompleted,
			Process:  "ALIGN",
			Name:     "ALIGN (sample-1)",
			Tag:      "sample-1",
			Workdir:  filepath.Clean(firstWorkdir),
			Exit:     &firstExit,
			Duration: "1m",
			Realtime: "60s",
			CPUs:     "2",
			Memory:   "4 GB",
		},
		{
			RowOrder: 2,
			ID:       "de/f456",
			Status:   domain.TaskStatusFailed,
			Process:  "QUANT",
			Name:     "QUANT (sample-2)",
			Tag:      "sample-2",
			Workdir:  filepath.Clean(secondWorkdir),
			Exit:     &secondExit,
			Duration: "2m",
			Realtime: "120s",
			CPUs:     "4",
			Memory:   "8 GB",
		},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("ParseTrace(csv fixture) = %#v, want %#v", got, want)
	}
}

func TestParseTraceParsesTSVFixture(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	workdir := filepath.Join(runDir.Path, "work", "12", "ABCDEF")
	tracePath := writeTraceFixture(t, runDir.Path, "trace.tsv", strings.Join([]string{
		"hash\tstatus\tprocess\tname\tworkdir\texit",
		"12/ABCDEF\tCACHED\tALIGN\tALIGN (cached)\t" + workdir + "\t-",
		"",
	}, "\n"))

	got, err := ParseTrace(context.Background(), runDir, domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: tracePath})
	if err != nil {
		t.Fatalf("ParseTrace(tsv fixture) returned error: %v", err)
	}

	want := []domain.Task{
		{
			RowOrder: 1,
			ID:       "12/abcdef",
			Status:   domain.TaskStatusCached,
			Process:  "ALIGN",
			Name:     "ALIGN (cached)",
			Workdir:  filepath.Clean(workdir),
		},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("ParseTrace(tsv fixture) = %#v, want %#v", got, want)
	}
}

func TestParseTraceParsesTXTFixtureWithMissingOptionalFields(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	tracePath := writeTraceFixture(t, runDir.Path, "trace.txt", strings.Join([]string{
		"hash\tstatus",
		"ABC123\tRUNNING",
		"",
	}, "\n"))

	got, err := ParseTrace(context.Background(), runDir, domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: tracePath})
	if err != nil {
		t.Fatalf("ParseTrace(txt fixture with missing optional fields) returned error: %v", err)
	}

	want := []domain.Task{
		{
			RowOrder: 1,
			ID:       "ab/c123",
			Status:   domain.TaskStatusRunning,
		},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("ParseTrace(txt fixture with missing optional fields) = %#v, want %#v", got, want)
	}
}

func TestParseTraceRejectsInvalidSource(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	validTracePath := writeTraceFixture(t, runDir.Path, "trace.txt", "hash\tstatus\nab/c123\tCOMPLETED\n")
	unsupportedTracePath := writeTraceFixture(t, runDir.Path, "trace.json", "hash\tstatus\nab/c123\tCOMPLETED\n")

	tests := []struct {
		name   string
		source domain.SourceFingerprint
		want   []string
	}{
		{
			name:   "wrong source kind",
			source: domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: validTracePath},
			want:   []string{"parse trace", "trace source", string(domain.SourceKindLog), string(domain.SourceKindTrace)},
		},
		{
			name:   "unsupported source extension",
			source: domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: unsupportedTracePath},
			want:   []string{"parse trace", "unsupported extension", ".json"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ParseTrace(context.Background(), runDir, tt.source)
			if err == nil {
				t.Fatalf("ParseTrace(%s) returned nil error and tasks %#v", tt.name, got)
			}
			if got != nil {
				t.Fatalf("tasks on error = %#v, want nil", got)
			}
			for _, want := range tt.want {
				if !strings.Contains(err.Error(), want) {
					t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
				}
			}
		})
	}
}

func TestDeriveCanonicalTaskIDAcceptsFullWorkdirPaths(t *testing.T) {
	tests := []struct {
		name string
		raw  string
		want string
	}{
		{name: "absolute workdir", raw: "/runs/example/work/ab/c123def", want: "ab/c123def"},
		{name: "relative workdir", raw: "work/de/f456", want: "de/f456"},
		{name: "nested file below workdir", raw: "/runs/example/work/12/abcdef/.command.err", want: "12/abcdef"},
		{name: "trailing slash", raw: "/runs/example/work/aa/bbbb/", want: "aa/bbbb"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := DeriveCanonicalTaskID(tt.raw)
			if err != nil {
				t.Fatalf("DeriveCanonicalTaskID(%q) returned error: %v", tt.raw, err)
			}
			if got != tt.want {
				t.Fatalf("DeriveCanonicalTaskID(%q) = %q, want %q", tt.raw, got, tt.want)
			}
		})
	}
}

func TestDeriveCanonicalTaskIDAcceptsDirectAndHashLikeValues(t *testing.T) {
	tests := []struct {
		name string
		raw  string
		want string
	}{
		{name: "direct canonical id", raw: "ab/c123def", want: "ab/c123def"},
		{name: "normalizes uppercase direct id", raw: " AB/C123DEF ", want: "ab/c123def"},
		{name: "unsplit hex hash", raw: "abc123def", want: "ab/c123def"},
		{name: "bracketed trace hash", raw: "[de/f456]", want: "de/f456"},
		{name: "quoted bracketed trace hash", raw: "\"[ABC123]\"", want: "ab/c123"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := DeriveCanonicalTaskID(tt.raw)
			if err != nil {
				t.Fatalf("DeriveCanonicalTaskID(%q) returned error: %v", tt.raw, err)
			}
			if got != tt.want {
				t.Fatalf("DeriveCanonicalTaskID(%q) = %q, want %q", tt.raw, got, tt.want)
			}
		})
	}
}

func TestDeriveCanonicalTaskIDRejectsEmptyAndUnparseableValuesClearly(t *testing.T) {
	tests := []struct {
		name string
		raw  string
		want []string
	}{
		{name: "empty", raw: "", want: []string{"derive canonical task id", "empty"}},
		{name: "blank", raw: "  \t\n ", want: []string{"derive canonical task id", "empty"}},
		{name: "missing work segment", raw: "/runs/example/ab/c123def", want: []string{"derive canonical task id", "workdir"}},
		{name: "non hex prefix", raw: "/runs/example/work/zz/c123def", want: []string{"derive canonical task id", "workdir"}},
		{name: "non hex suffix", raw: "ab/not-a-hash", want: []string{"derive canonical task id", "workdir"}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := DeriveCanonicalTaskID(tt.raw)
			if err == nil {
				t.Fatalf("DeriveCanonicalTaskID(%q) returned nil error and id %q", tt.raw, got)
			}
			if got != "" {
				t.Fatalf("id on error = %q, want empty", got)
			}
			for _, want := range tt.want {
				if !strings.Contains(err.Error(), want) {
					t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
				}
			}
		})
	}
}

func TestResolveTaskWorkdirUsesProvidedFullPathWithoutRequiringItToExist(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	workdir := filepath.Join(t.TempDir(), "custom-work-root", "ab", "c123def")

	got, err := ResolveTaskWorkdir(runDir, "ab/c123def", workdir)
	if err != nil {
		t.Fatalf("ResolveTaskWorkdir() returned error: %v", err)
	}
	if got != workdir {
		t.Fatalf("ResolveTaskWorkdir() = %q, want provided full path %q", got, workdir)
	}
	if _, statErr := os.Stat(workdir); !os.IsNotExist(statErr) {
		t.Fatalf("ResolveTaskWorkdir() should not create provided workdir; stat error = %v", statErr)
	}
}

func TestResolveTaskWorkdirDerivesExistingRunWorkdirFromHashOnlyValue(t *testing.T) {
	runDir := t.TempDir()
	want := filepath.Join(runDir, "work", "ab", "c123def")
	if err := os.MkdirAll(want, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q): %v", want, err)
	}

	got, err := ResolveTaskWorkdir(domain.RunDir{Path: runDir}, "ab/c123def", "ab/c123def")
	if err != nil {
		t.Fatalf("ResolveTaskWorkdir() returned error: %v", err)
	}
	if got != want {
		t.Fatalf("ResolveTaskWorkdir() = %q, want existing derived path %q", got, want)
	}
}

func TestResolveTaskWorkdirDerivesExistingRunWorkdirFromUnsplitRawHash(t *testing.T) {
	runDir := t.TempDir()
	want := filepath.Join(runDir, "work", "de", "f456")
	if err := os.MkdirAll(want, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q): %v", want, err)
	}

	got, err := ResolveTaskWorkdir(domain.RunDir{Path: runDir}, "", "DEF456")
	if err != nil {
		t.Fatalf("ResolveTaskWorkdir() returned error: %v", err)
	}
	if got != want {
		t.Fatalf("ResolveTaskWorkdir() = %q, want existing derived path %q", got, want)
	}
}

func TestResolveTaskWorkdirReturnsUnknownForHashOnlyValueWhenDerivedPathIsMissing(t *testing.T) {
	runDir := t.TempDir()
	missing := filepath.Join(runDir, "work", "ab", "c123def")

	got, err := ResolveTaskWorkdir(domain.RunDir{Path: runDir}, "ab/c123def", "ab/c123def")
	if err != nil {
		t.Fatalf("ResolveTaskWorkdir() returned error: %v", err)
	}
	if got != "" {
		t.Fatalf("ResolveTaskWorkdir() = %q, want empty unknown workdir", got)
	}
	if _, statErr := os.Stat(missing); !os.IsNotExist(statErr) {
		t.Fatalf("ResolveTaskWorkdir() should not create derived workdir; stat error = %v", statErr)
	}
}

func TestResolveTaskWorkdirReturnsUnknownForBlankValueWhenDerivedPathIsMissing(t *testing.T) {
	runDir := t.TempDir()

	got, err := ResolveTaskWorkdir(domain.RunDir{Path: runDir}, "ab/c123def", "  \t\n  ")
	if err != nil {
		t.Fatalf("ResolveTaskWorkdir() returned error: %v", err)
	}
	if got != "" {
		t.Fatalf("ResolveTaskWorkdir() = %q, want empty unknown workdir", got)
	}
}

func TestNormalizeTraceRecordMapsCommonTraceColumns(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	workdir := filepath.Join(t.TempDir(), "custom-work-root", "AB", "C123DEF")
	record := RawRecord{
		RowOrder: 7,
		Columns: map[string]string{
			"hash":     "AB/C123DEF",
			"workdir":  workdir,
			"status":   " failed ",
			"process":  "ALIGN_STAR",
			"name":     "ALIGN_STAR (sample-1)",
			"tag":      "sample-1",
			"exit":     " 137 ",
			"duration": "1h 2m",
			"realtime": "62m",
			"cpus":     "8",
			"memory":   "16 GB",
		},
	}

	got, err := NormalizeTraceRecord(runDir, record)
	if err != nil {
		t.Fatalf("NormalizeTraceRecord() returned error: %v", err)
	}

	wantExit := 137
	want := domain.Task{
		RowOrder: 7,
		ID:       "ab/c123def",
		Status:   domain.TaskStatusFailed,
		Process:  "ALIGN_STAR",
		Name:     "ALIGN_STAR (sample-1)",
		Tag:      "sample-1",
		Workdir:  filepath.Clean(workdir),
		Exit:     &wantExit,
		Duration: "1h 2m",
		Realtime: "62m",
		CPUs:     "8",
		Memory:   "16 GB",
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("NormalizeTraceRecord() = %#v, want %#v", got, want)
	}
}

func TestNormalizeTraceRecordMapsNFCoreExecutionTraceAliases(t *testing.T) {
	runDir := t.TempDir()
	wantWorkdir := filepath.Join(runDir, "work", "de", "f456")
	if err := os.MkdirAll(wantWorkdir, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q): %v", wantWorkdir, err)
	}

	fullName := "NFGOS:AMBER_STEP:BAM_AMBER:AMBER (WG-26-03_vs_WG-26-04)"
	record := RawRecord{
		RowOrder: 12,
		Columns: map[string]string{
			"task_id":   "42",
			"hash":      "DE/F456",
			"native_id": "12345",
			"name":      fullName,
			"status":    " cached ",
			"exit":      "0",
			"duration":  "3m",
			"realtime":  "2m",
			"%cpu":      "87.5%",
			"peak_rss":  "4.5 GB",
			"peak_vmem": "10 GB",
		},
	}

	got, err := NormalizeTraceRecord(domain.RunDir{Path: runDir}, record)
	if err != nil {
		t.Fatalf("NormalizeTraceRecord() returned error: %v", err)
	}

	wantExit := 0
	want := domain.Task{
		RowOrder: 12,
		ID:       "de/f456",
		Status:   domain.TaskStatusCached,
		Process:  "NFGOS:AMBER_STEP:BAM_AMBER:AMBER",
		Name:     fullName,
		Tag:      "WG-26-03_vs_WG-26-04",
		Workdir:  wantWorkdir,
		Exit:     &wantExit,
		Duration: "3m",
		Realtime: "2m",
		CPUs:     "87.5%",
		Memory:   "peak_rss=4.5 GB; peak_vmem=10 GB",
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("NormalizeTraceRecord() = %#v, want %#v", got, want)
	}
}

func TestNormalizeTraceRecordUsesWorkdirWhenHashColumnIsMissing(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	workdir := filepath.Join(t.TempDir(), "work", "DE", "F456")
	record := RawRecord{
		RowOrder: 2,
		Columns: map[string]string{
			"workdir": workdir,
			"status":  "COMPLETED",
		},
	}

	got, err := NormalizeTraceRecord(runDir, record)
	if err != nil {
		t.Fatalf("NormalizeTraceRecord() returned error: %v", err)
	}
	if got.ID != "de/f456" {
		t.Fatalf("task ID = %q, want canonical ID from workdir", got.ID)
	}
	if got.Workdir != filepath.Clean(workdir) {
		t.Fatalf("task workdir = %q, want %q", got.Workdir, filepath.Clean(workdir))
	}
	if got.Status != domain.TaskStatusCompleted {
		t.Fatalf("task status = %q, want %q", got.Status, domain.TaskStatusCompleted)
	}
}

func TestNormalizeTraceRecordDefaultsMissingOptionalColumns(t *testing.T) {
	runDir := t.TempDir()
	wantWorkdir := filepath.Join(runDir, "work", "de", "f456")
	if err := os.MkdirAll(wantWorkdir, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q): %v", wantWorkdir, err)
	}

	got, err := NormalizeTraceRecord(domain.RunDir{Path: runDir}, RawRecord{
		RowOrder: 3,
		Columns: map[string]string{
			"hash": "DEF456",
		},
	})
	if err != nil {
		t.Fatalf("NormalizeTraceRecord() returned error: %v", err)
	}

	want := domain.Task{
		RowOrder: 3,
		ID:       "de/f456",
		Status:   domain.TaskStatusUnknown,
		Workdir:  wantWorkdir,
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("NormalizeTraceRecord() = %#v, want %#v", got, want)
	}
}

func TestNormalizeTraceRecordRejectsRecordsWithoutParseableHashOrWorkdir(t *testing.T) {
	tests := []struct {
		name    string
		columns map[string]string
	}{
		{name: "missing id columns", columns: map[string]string{"status": "COMPLETED"}},
		{name: "invalid id columns", columns: map[string]string{"hash": "not-a-hash", "workdir": "also-not-a-path"}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := NormalizeTraceRecord(domain.RunDir{Path: t.TempDir()}, RawRecord{RowOrder: 9, Columns: tt.columns})
			if err == nil {
				t.Fatalf("NormalizeTraceRecord() returned nil error and task %#v", got)
			}
			if !reflect.DeepEqual(got, domain.Task{}) {
				t.Fatalf("task on error = %#v, want zero value", got)
			}
			for _, want := range []string{"normalize trace record", "row 9", "hash", "workdir"} {
				if !strings.Contains(err.Error(), want) {
					t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
				}
			}
		})
	}
}

func TestNormalizeTraceRecordPropagatesInvalidExitValues(t *testing.T) {
	got, err := NormalizeTraceRecord(domain.RunDir{Path: t.TempDir()}, RawRecord{
		RowOrder: 4,
		Columns: map[string]string{
			"hash": "ab/c123def",
			"exit": "killed",
		},
	})
	if err == nil {
		t.Fatalf("NormalizeTraceRecord() returned nil error and task %#v", got)
	}
	if !reflect.DeepEqual(got, domain.Task{}) {
		t.Fatalf("task on error = %#v, want zero value", got)
	}
	for _, want := range []string{"normalize trace record", "row 4", "exit", "killed"} {
		if !strings.Contains(err.Error(), want) {
			t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
		}
	}
}

func TestParseTraceRecordsParsesCommaRowsInSourceOrder(t *testing.T) {
	input := strings.Join([]string{
		"task_id,status,workdir",
		"1,COMPLETED,/runs/work/ab/c123",
		"2,FAILED,/runs/work/de/f456",
		"",
	}, "\n")

	got, err := ParseTraceRecords(strings.NewReader(input), DelimiterComma)
	if err != nil {
		t.Fatalf("ParseTraceRecords returned error: %v", err)
	}
	if len(got) != 2 {
		t.Fatalf("record count = %d, want 2", len(got))
	}

	wantFirst := RawRecord{
		RowOrder: 1,
		Columns: map[string]string{
			"task_id": "1",
			"status":  "COMPLETED",
			"workdir": "/runs/work/ab/c123",
		},
	}
	wantSecond := RawRecord{
		RowOrder: 2,
		Columns: map[string]string{
			"task_id": "2",
			"status":  "FAILED",
			"workdir": "/runs/work/de/f456",
		},
	}
	if !reflect.DeepEqual(got[0], wantFirst) {
		t.Fatalf("first record = %#v, want %#v", got[0], wantFirst)
	}
	if !reflect.DeepEqual(got[1], wantSecond) {
		t.Fatalf("second record = %#v, want %#v", got[1], wantSecond)
	}
}

func TestParseTraceRecordsParsesTabDelimitedRows(t *testing.T) {
	input := strings.Join([]string{
		"process\tname\tstatus",
		"ALIGN\tSAMPLE_T1\tCACHED",
		"",
	}, "\n")

	got, err := ParseTraceRecords(strings.NewReader(input), DelimiterTab)
	if err != nil {
		t.Fatalf("ParseTraceRecords returned error: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("record count = %d, want 1", len(got))
	}
	want := RawRecord{
		RowOrder: 1,
		Columns: map[string]string{
			"process": "ALIGN",
			"name":    "SAMPLE_T1",
			"status":  "CACHED",
		},
	}
	if !reflect.DeepEqual(got[0], want) {
		t.Fatalf("record = %#v, want %#v", got[0], want)
	}
}

func TestParseTraceRecordsDoesNotRequireOptionalColumns(t *testing.T) {
	input := strings.Join([]string{
		"status",
		"COMPLETED",
		"FAILED",
		"",
	}, "\n")

	got, err := ParseTraceRecords(strings.NewReader(input), DelimiterComma)
	if err != nil {
		t.Fatalf("ParseTraceRecords returned error for minimal trace columns: %v", err)
	}
	if len(got) != 2 {
		t.Fatalf("record count = %d, want 2", len(got))
	}
	if got[0].RowOrder != 1 || got[0].Columns["status"] != "COMPLETED" {
		t.Fatalf("first record = %#v, want row order 1 and status COMPLETED", got[0])
	}
	if got[1].RowOrder != 2 || got[1].Columns["status"] != "FAILED" {
		t.Fatalf("second record = %#v, want row order 2 and status FAILED", got[1])
	}
}

func TestParseTraceRecordsRejectsMissingHeaderClearly(t *testing.T) {
	got, err := ParseTraceRecords(strings.NewReader(""), DelimiterComma)
	if err == nil {
		t.Fatalf("ParseTraceRecords(empty input) returned nil error")
	}
	if got != nil {
		t.Fatalf("records on error = %#v, want nil", got)
	}
	for _, want := range []string{"parse trace records", "header"} {
		if !strings.Contains(err.Error(), want) {
			t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
		}
	}
}

func TestParseTraceRecordsRejectsMalformedHeadersClearly(t *testing.T) {
	tests := []struct {
		name  string
		input string
		want  []string
	}{
		{
			name:  "empty column name",
			input: "task_id,,status\n1,x,COMPLETED\n",
			want:  []string{"header", "empty", "column 2"},
		},
		{
			name:  "duplicate column name",
			input: "status,status\nCOMPLETED,FAILED\n",
			want:  []string{"header", "duplicate", "status"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ParseTraceRecords(strings.NewReader(tt.input), DelimiterComma)
			if err == nil {
				t.Fatalf("ParseTraceRecords(%s) returned nil error and records %#v", tt.name, got)
			}
			if got != nil {
				t.Fatalf("records on error = %#v, want nil", got)
			}
			for _, want := range tt.want {
				if !strings.Contains(err.Error(), want) {
					t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
				}
			}
		})
	}
}

func TestParseTraceRecordsRejectsMalformedRowsClearly(t *testing.T) {
	input := strings.Join([]string{
		"status,workdir",
		"COMPLETED,/runs/work/ab/c123",
		"FAILED",
		"",
	}, "\n")

	got, err := ParseTraceRecords(strings.NewReader(input), DelimiterComma)
	if err == nil {
		t.Fatalf("ParseTraceRecords(malformed row) returned nil error and records %#v", got)
	}
	if got != nil {
		t.Fatalf("records on error = %#v, want nil", got)
	}
	for _, want := range []string{"row 2", "wrong number"} {
		if !strings.Contains(err.Error(), want) {
			t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
		}
	}
}

func TestParseTraceRecordsRejectsUnsupportedDelimiter(t *testing.T) {
	got, err := ParseTraceRecords(strings.NewReader("status\nCOMPLETED\n"), Delimiter('|'))
	if err == nil {
		t.Fatalf("ParseTraceRecords(unsupported delimiter) returned nil error and records %#v", got)
	}
	if got != nil {
		t.Fatalf("records on error = %#v, want nil", got)
	}
	for _, want := range []string{"unsupported delimiter", ",", "tab"} {
		if !strings.Contains(err.Error(), want) {
			t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
		}
	}
}

func writeTraceFixture(t *testing.T, dir string, name string, content string) string {
	t.Helper()

	path := filepath.Join(dir, name)
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile(%q): %v", path, err)
	}
	return path
}
