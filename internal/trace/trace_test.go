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
