package inspect

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/mskilab-org/gosh/internal/domain"
)

func TestCommandFilePathJoinsWorkdirWithEveryKnownCommandKind(t *testing.T) {
	workdir := filepath.FromSlash("/runs/example/work/ab/c123def")

	tests := []struct {
		name string
		kind domain.CommandFileKind
		want string
	}{
		{name: "shell script", kind: domain.CommandFileShell, want: filepath.FromSlash("/runs/example/work/ab/c123def/.command.sh")},
		{name: "combined log", kind: domain.CommandFileLog, want: filepath.FromSlash("/runs/example/work/ab/c123def/.command.log")},
		{name: "stderr", kind: domain.CommandFileErr, want: filepath.FromSlash("/runs/example/work/ab/c123def/.command.err")},
		{name: "stdout", kind: domain.CommandFileOut, want: filepath.FromSlash("/runs/example/work/ab/c123def/.command.out")},
		{name: "run wrapper", kind: domain.CommandFileRun, want: filepath.FromSlash("/runs/example/work/ab/c123def/.command.run")},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := CommandFilePath(workdir, tt.kind)
			if got != tt.want {
				t.Fatalf("CommandFilePath(%q, %q) = %q, want %q", workdir, tt.kind, got, tt.want)
			}
		})
	}
}

func TestCommandFilePathCleansTrailingWorkdirSeparator(t *testing.T) {
	workdir := filepath.FromSlash("/runs/example/work/de/f456/")
	want := filepath.FromSlash("/runs/example/work/de/f456/.command.err")

	got := CommandFilePath(workdir, domain.CommandFileErr)
	if got != want {
		t.Fatalf("CommandFilePath(%q, %q) = %q, want %q", workdir, domain.CommandFileErr, got, want)
	}
}

func TestCommandFilePathHandlesRelativeAndEmptyWorkdirs(t *testing.T) {
	tests := []struct {
		name    string
		workdir string
		kind    domain.CommandFileKind
		want    string
	}{
		{name: "relative workdir", workdir: filepath.FromSlash("work/12/abcdef"), kind: domain.CommandFileOut, want: filepath.FromSlash("work/12/abcdef/.command.out")},
		{name: "empty workdir falls back to relative command file", workdir: "", kind: domain.CommandFileRun, want: ".command.run"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := CommandFilePath(tt.workdir, tt.kind)
			if got != tt.want {
				t.Fatalf("CommandFilePath(%q, %q) = %q, want %q", tt.workdir, tt.kind, got, tt.want)
			}
		})
	}
}

func TestSnippetStrategyForKindMapsKnownCommandKinds(t *testing.T) {
	tests := []struct {
		name string
		kind domain.CommandFileKind
		want domain.SnippetStrategy
	}{
		{name: "shell script uses head", kind: domain.CommandFileShell, want: domain.SnippetStrategyHead},
		{name: "combined log uses error-focused", kind: domain.CommandFileLog, want: domain.SnippetStrategyError},
		{name: "stderr uses error-focused", kind: domain.CommandFileErr, want: domain.SnippetStrategyError},
		{name: "stdout uses tail", kind: domain.CommandFileOut, want: domain.SnippetStrategyTail},
		{name: "run wrapper uses head", kind: domain.CommandFileRun, want: domain.SnippetStrategyHead},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := SnippetStrategyForKind(tt.kind)
			if got != tt.want {
				t.Fatalf("SnippetStrategyForKind(%q) = %q, want %q", tt.kind, got, tt.want)
			}
		})
	}
}

func TestSnippetStrategyForKindDefaultsUnknownKindsToTail(t *testing.T) {
	unknownKind := domain.CommandFileKind(".command.custom")

	got := SnippetStrategyForKind(unknownKind)
	if got != domain.SnippetStrategyTail {
		t.Fatalf("SnippetStrategyForKind(%q) = %q, want %q", unknownKind, got, domain.SnippetStrategyTail)
	}
}

func TestInventoryCommandFilesListsAllExpectedKindsAndSnippetsExistingFiles(t *testing.T) {
	workdir := t.TempDir()
	markerPath := filepath.Join(workdir, "executed-marker")
	shellContent := strings.Join([]string{
		"#!/usr/bin/env bash",
		"touch " + markerPath,
		"echo after",
	}, "\n")
	logContent := strings.Join([]string{
		"setup",
		"fatal ERROR writing sample",
		"cleanup failed later",
		"done",
	}, "\n")
	outContent := strings.Join([]string{"one", "two", "three"}, "\n")

	for _, fixture := range []struct {
		name    string
		content string
		mode    os.FileMode
	}{
		{name: string(domain.CommandFileShell), content: shellContent, mode: 0o755},
		{name: string(domain.CommandFileLog), content: logContent, mode: 0o644},
		{name: string(domain.CommandFileOut), content: outContent, mode: 0o644},
	} {
		path := filepath.Join(workdir, fixture.name)
		if err := os.WriteFile(path, []byte(fixture.content), fixture.mode); err != nil {
			t.Fatalf("write inventory fixture %s: %v", fixture.name, err)
		}
	}

	got, err := InventoryCommandFiles(context.Background(), workdir, SnippetOptions{MaxBytes: 4096, MaxLines: 2})
	if err != nil {
		t.Fatalf("InventoryCommandFiles returned error: %v", err)
	}

	if got.Workdir != workdir {
		t.Fatalf("inventory workdir = %q, want %q", got.Workdir, workdir)
	}
	if len(got.Files) != 5 {
		t.Fatalf("inventory files len = %d, want 5: %#v", len(got.Files), got.Files)
	}

	expectedKinds := []domain.CommandFileKind{
		domain.CommandFileShell,
		domain.CommandFileLog,
		domain.CommandFileErr,
		domain.CommandFileOut,
		domain.CommandFileRun,
	}
	for index, kind := range expectedKinds {
		file := got.Files[index]
		wantPath := CommandFilePath(workdir, kind)
		if file.Kind != kind {
			t.Fatalf("file[%d] kind = %q, want %q", index, file.Kind, kind)
		}
		if file.Path != wantPath {
			t.Fatalf("file[%d] path = %q, want %q", index, file.Path, wantPath)
		}
	}

	if !got.Files[0].Exists || got.Files[0].Size != int64(len([]byte(shellContent))) {
		t.Fatalf("shell file metadata = exists:%t size:%d, want exists:true size:%d", got.Files[0].Exists, got.Files[0].Size, len([]byte(shellContent)))
	}
	assertSnippet(t, got.Files[0].Snippet, CommandFilePath(workdir, domain.CommandFileShell), domain.SnippetStrategyHead, 1, 2, "#!/usr/bin/env bash\ntouch "+markerPath, true, 4096)

	if !got.Files[1].Exists || got.Files[1].Size != int64(len([]byte(logContent))) {
		t.Fatalf("log file metadata = exists:%t size:%d, want exists:true size:%d", got.Files[1].Exists, got.Files[1].Size, len([]byte(logContent)))
	}
	assertSnippet(t, got.Files[1].Snippet, CommandFilePath(workdir, domain.CommandFileLog), domain.SnippetStrategyError, 2, 3, "fatal ERROR writing sample\ncleanup failed later", true, 4096)

	if got.Files[2].Exists || got.Files[2].Size != 0 || got.Files[2].Snippet != nil {
		t.Fatalf("missing err file = %#v, want exists false, zero size, nil snippet", got.Files[2])
	}

	if !got.Files[3].Exists || got.Files[3].Size != int64(len([]byte(outContent))) {
		t.Fatalf("out file metadata = exists:%t size:%d, want exists:true size:%d", got.Files[3].Exists, got.Files[3].Size, len([]byte(outContent)))
	}
	assertSnippet(t, got.Files[3].Snippet, CommandFilePath(workdir, domain.CommandFileOut), domain.SnippetStrategyTail, 2, 3, "two\nthree", true, 4096)

	if got.Files[4].Exists || got.Files[4].Size != 0 || got.Files[4].Snippet != nil {
		t.Fatalf("missing run file = %#v, want exists false, zero size, nil snippet", got.Files[4])
	}
	if _, err := os.Stat(markerPath); !errors.Is(err, os.ErrNotExist) {
		t.Fatalf("marker stat error = %v, want marker not to exist because command files must not be executed", err)
	}
}

func TestInventoryCommandFilesReportsMissingExpectedPathsWhenWorkdirHasNoCommandFiles(t *testing.T) {
	workdir := t.TempDir()

	got, err := InventoryCommandFiles(context.Background(), workdir, SnippetOptions{MaxBytes: 64, MaxLines: 2})
	if err != nil {
		t.Fatalf("InventoryCommandFiles(empty workdir) returned error: %v", err)
	}

	if got.Workdir != workdir {
		t.Fatalf("inventory workdir = %q, want %q", got.Workdir, workdir)
	}
	expectedKinds := []domain.CommandFileKind{
		domain.CommandFileShell,
		domain.CommandFileLog,
		domain.CommandFileErr,
		domain.CommandFileOut,
		domain.CommandFileRun,
	}
	if len(got.Files) != len(expectedKinds) {
		t.Fatalf("inventory files len = %d, want %d: %#v", len(got.Files), len(expectedKinds), got.Files)
	}
	for index, kind := range expectedKinds {
		file := got.Files[index]
		if file.Kind != kind {
			t.Fatalf("file[%d] kind = %q, want %q", index, file.Kind, kind)
		}
		if file.Path != CommandFilePath(workdir, kind) {
			t.Fatalf("file[%d] path = %q, want %q", index, file.Path, CommandFilePath(workdir, kind))
		}
		if file.Exists || file.Size != 0 || file.Snippet != nil {
			t.Fatalf("file[%d] = %#v, want exists false, zero size, nil snippet", index, file)
		}
	}
}

func TestInventoryCommandFilesReturnsContextErrorBeforeReading(t *testing.T) {
	workdir := t.TempDir()
	if err := os.WriteFile(filepath.Join(workdir, string(domain.CommandFileLog)), []byte("error should not be read"), 0o644); err != nil {
		t.Fatalf("write command fixture: %v", err)
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()

	got, err := InventoryCommandFiles(ctx, workdir, SnippetOptions{MaxBytes: 64, MaxLines: 2})
	if err == nil {
		t.Fatalf("InventoryCommandFiles(cancelled context) returned nil error and inventory %#v", got)
	}
	if !errors.Is(err, context.Canceled) {
		t.Fatalf("error = %v, want context.Canceled", err)
	}
}

func TestReadBoundedSnippetHeadRespectsLineLimitAndDoesNotExecute(t *testing.T) {
	workspace := t.TempDir()
	commandPath := filepath.Join(workspace, ".command.sh")
	markerPath := filepath.Join(workspace, "executed-marker")
	content := strings.Join([]string{
		"#!/usr/bin/env bash",
		"touch " + markerPath,
		"echo after",
	}, "\n")
	if err := os.WriteFile(commandPath, []byte(content), 0o755); err != nil {
		t.Fatalf("write command fixture: %v", err)
	}

	got, err := ReadBoundedSnippet(context.Background(), commandPath, domain.SnippetStrategyHead, SnippetOptions{MaxBytes: 4096, MaxLines: 2})
	if err != nil {
		t.Fatalf("ReadBoundedSnippet(head) returned error: %v", err)
	}
	if got == nil {
		t.Fatalf("ReadBoundedSnippet(head) returned nil snippet")
	}

	wantContent := strings.Join([]string{
		"#!/usr/bin/env bash",
		"touch " + markerPath,
	}, "\n")
	assertSnippet(t, got, commandPath, domain.SnippetStrategyHead, 1, 2, wantContent, true, 4096)
	if _, err := os.Stat(markerPath); !errors.Is(err, os.ErrNotExist) {
		t.Fatalf("marker stat error = %v, want marker not to exist because command files must not be executed", err)
	}
}

func TestReadBoundedSnippetTailReturnsLastLines(t *testing.T) {
	commandPath := writeCommandFixture(t, ".command.out", strings.Join([]string{"one", "two", "three", "four"}, "\n"))

	got, err := ReadBoundedSnippet(context.Background(), commandPath, domain.SnippetStrategyTail, SnippetOptions{MaxBytes: 64, MaxLines: 2})
	if err != nil {
		t.Fatalf("ReadBoundedSnippet(tail) returned error: %v", err)
	}

	assertSnippet(t, got, commandPath, domain.SnippetStrategyTail, 3, 4, "three\nfour", true, 64)
}

func TestReadBoundedSnippetErrorFocusesWindowAroundFirstErrorTerm(t *testing.T) {
	commandPath := writeCommandFixture(t, ".command.log", strings.Join([]string{
		"setup",
		"download ok",
		"warning only",
		"fatal ERROR writing sample",
		"cleanup failed later",
		"done",
	}, "\n"))

	got, err := ReadBoundedSnippet(context.Background(), commandPath, domain.SnippetStrategyError, SnippetOptions{MaxBytes: 200, MaxLines: 3})
	if err != nil {
		t.Fatalf("ReadBoundedSnippet(error-focused) returned error: %v", err)
	}

	wantContent := strings.Join([]string{
		"warning only",
		"fatal ERROR writing sample",
		"cleanup failed later",
	}, "\n")
	assertSnippet(t, got, commandPath, domain.SnippetStrategyError, 3, 5, wantContent, true, 200)
}

func TestReadBoundedSnippetErrorFallsBackToTailWhenNoTermMatches(t *testing.T) {
	commandPath := writeCommandFixture(t, ".command.err", strings.Join([]string{"first", "second", "third", "fourth"}, "\n"))

	got, err := ReadBoundedSnippet(context.Background(), commandPath, domain.SnippetStrategyError, SnippetOptions{MaxBytes: 64, MaxLines: 2})
	if err != nil {
		t.Fatalf("ReadBoundedSnippet(error-focused without error term) returned error: %v", err)
	}

	assertSnippet(t, got, commandPath, domain.SnippetStrategyError, 3, 4, "third\nfourth", true, 64)
}

func TestReadBoundedSnippetRespectsByteLimitWithoutSplittingLinesWhenPossible(t *testing.T) {
	commandPath := writeCommandFixture(t, ".command.log", strings.Join([]string{"alpha", "bravo", "charlie"}, "\n"))

	got, err := ReadBoundedSnippet(context.Background(), commandPath, domain.SnippetStrategyHead, SnippetOptions{MaxBytes: 11, MaxLines: 3})
	if err != nil {
		t.Fatalf("ReadBoundedSnippet(head with byte limit) returned error: %v", err)
	}

	assertSnippet(t, got, commandPath, domain.SnippetStrategyHead, 1, 2, "alpha\nbravo", true, 11)
}

func TestReadBoundedSnippetReturnsNilSnippetAndClearErrorForMissingFile(t *testing.T) {
	missingPath := filepath.Join(t.TempDir(), ".command.log")

	got, err := ReadBoundedSnippet(context.Background(), missingPath, domain.SnippetStrategyTail, SnippetOptions{MaxBytes: 64, MaxLines: 2})
	if err == nil {
		t.Fatalf("ReadBoundedSnippet(missing file) returned nil error and snippet %#v", got)
	}
	if got != nil {
		t.Fatalf("snippet for missing file = %#v, want nil", got)
	}
	if !errors.Is(err, os.ErrNotExist) {
		t.Fatalf("error = %v, want it to wrap os.ErrNotExist", err)
	}
	if !strings.Contains(strings.ToLower(err.Error()), "missing") {
		t.Fatalf("error = %q, want it to clearly mention missing file", err.Error())
	}
}

func TestReadBoundedSnippetReturnsContextErrorBeforeReading(t *testing.T) {
	commandPath := writeCommandFixture(t, ".command.log", "line one\nline two")
	ctx, cancel := context.WithCancel(context.Background())
	cancel()

	got, err := ReadBoundedSnippet(ctx, commandPath, domain.SnippetStrategyHead, SnippetOptions{MaxBytes: 64, MaxLines: 2})
	if err == nil {
		t.Fatalf("ReadBoundedSnippet(cancelled context) returned nil error and snippet %#v", got)
	}
	if got != nil {
		t.Fatalf("snippet for cancelled context = %#v, want nil", got)
	}
	if !errors.Is(err, context.Canceled) {
		t.Fatalf("error = %v, want context.Canceled", err)
	}
}

func TestReadBoundedSnippetRejectsUnboundedLimits(t *testing.T) {
	commandPath := writeCommandFixture(t, ".command.log", "line one\nline two")

	tests := []struct {
		name      string
		options   SnippetOptions
		wantError string
	}{
		{name: "zero max bytes", options: SnippetOptions{MaxBytes: 0, MaxLines: 2}, wantError: "max bytes"},
		{name: "negative max bytes", options: SnippetOptions{MaxBytes: -1, MaxLines: 2}, wantError: "max bytes"},
		{name: "zero max lines", options: SnippetOptions{MaxBytes: 64, MaxLines: 0}, wantError: "max lines"},
		{name: "negative max lines", options: SnippetOptions{MaxBytes: 64, MaxLines: -1}, wantError: "max lines"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ReadBoundedSnippet(context.Background(), commandPath, domain.SnippetStrategyHead, tt.options)
			if err == nil {
				t.Fatalf("ReadBoundedSnippet(%+v) returned nil error and snippet %#v", tt.options, got)
			}
			if got != nil {
				t.Fatalf("snippet for invalid options = %#v, want nil", got)
			}
			if !strings.Contains(strings.ToLower(err.Error()), tt.wantError) {
				t.Fatalf("error = %q, want it to mention %q", err.Error(), tt.wantError)
			}
		})
	}
}

func writeCommandFixture(t *testing.T, name string, content string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("write command fixture: %v", err)
	}
	return path
}

func assertSnippet(t *testing.T, got *domain.Snippet, wantPath string, wantStrategy domain.SnippetStrategy, wantStartLine int, wantEndLine int, wantContent string, wantTruncated bool, wantMaxBytes int64) {
	t.Helper()
	if got == nil {
		t.Fatalf("snippet = nil, want populated snippet")
	}
	if got.Path != wantPath {
		t.Fatalf("snippet path = %q, want %q", got.Path, wantPath)
	}
	if got.Strategy != wantStrategy {
		t.Fatalf("snippet strategy = %q, want %q", got.Strategy, wantStrategy)
	}
	if got.StartLine != wantStartLine || got.EndLine != wantEndLine {
		t.Fatalf("snippet lines = %d-%d, want %d-%d", got.StartLine, got.EndLine, wantStartLine, wantEndLine)
	}
	if got.Content != wantContent {
		t.Fatalf("snippet content = %q, want %q", got.Content, wantContent)
	}
	if got.Truncated != wantTruncated {
		t.Fatalf("snippet truncated = %t, want %t", got.Truncated, wantTruncated)
	}
	if got.MaxBytes != wantMaxBytes {
		t.Fatalf("snippet max bytes = %d, want %d", got.MaxBytes, wantMaxBytes)
	}
	if int64(len([]byte(got.Content))) > wantMaxBytes {
		t.Fatalf("snippet content length = %d bytes, want <= %d", len([]byte(got.Content)), wantMaxBytes)
	}
}
