package index

import (
	"context"
	"database/sql"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"reflect"
	"sort"
	"strings"
	"testing"
	"time"

	"github.com/mskilab-org/gosh/internal/domain"
	"github.com/mskilab-org/gosh/internal/run"
)

func TestOpenStoreUsesRunIndexPathCreatesCacheDirectoryAndLeavesArtifactsUntouched(t *testing.T) {
	runRoot := t.TempDir()
	runDir := domain.RunDir{Path: runRoot}
	artifactPath := filepath.Join(runRoot, ".nextflow.log")
	artifactContent := []byte("log stays untouched\n")
	if err := os.WriteFile(artifactPath, artifactContent, 0o644); err != nil {
		t.Fatalf("write existing run artifact: %v", err)
	}

	store, err := OpenStore(context.Background(), runDir)
	if err != nil {
		t.Fatalf("OpenStore() returned error: %v", err)
	}
	defer store.Close()

	wantPath := run.IndexPath(runDir)
	if store.Path != wantPath {
		t.Fatalf("Store.Path = %q, want %q", store.Path, wantPath)
	}
	if store.DB == nil {
		t.Fatalf("Store.DB = nil, want an open sqlite database")
	}
	if err := store.DB.PingContext(context.Background()); err != nil {
		t.Fatalf("ping opened store: %v", err)
	}

	cacheInfo, err := os.Stat(filepath.Join(runRoot, run.IndexDirName))
	if err != nil {
		t.Fatalf("stat cache directory: %v", err)
	}
	if !cacheInfo.IsDir() {
		t.Fatalf("cache path is not a directory")
	}
	indexInfo, err := os.Stat(wantPath)
	if err != nil {
		t.Fatalf("stat sqlite index: %v", err)
	}
	if !indexInfo.Mode().IsRegular() {
		t.Fatalf("sqlite index is not a regular file")
	}

	gotArtifactContent, err := os.ReadFile(artifactPath)
	if err != nil {
		t.Fatalf("read existing run artifact: %v", err)
	}
	if string(gotArtifactContent) != string(artifactContent) {
		t.Fatalf("existing run artifact content = %q, want %q", gotArtifactContent, artifactContent)
	}
}

func TestOpenStoreRejectsCanceledContextWithoutCreatingCacheDirectory(t *testing.T) {
	runRoot := t.TempDir()
	ctx, cancel := context.WithCancel(context.Background())
	cancel()

	store, err := OpenStore(ctx, domain.RunDir{Path: runRoot})
	if err == nil {
		t.Fatalf("OpenStore() returned nil error for canceled context")
	}
	if !errors.Is(err, context.Canceled) {
		t.Fatalf("OpenStore() error = %v, want context.Canceled", err)
	}
	if store != nil {
		defer store.Close()
		t.Fatalf("OpenStore() store = %#v, want nil on error", store)
	}
	if _, statErr := os.Stat(filepath.Join(runRoot, run.IndexDirName)); !errors.Is(statErr, os.ErrNotExist) {
		t.Fatalf("cache directory stat error = %v, want os.ErrNotExist", statErr)
	}
}

func TestOpenStoreRejectsNilContextWithoutCreatingCacheDirectory(t *testing.T) {
	runRoot := t.TempDir()

	store, err := OpenStore(nil, domain.RunDir{Path: runRoot})
	if err == nil {
		t.Fatalf("OpenStore() returned nil error for nil context")
	}
	if !strings.Contains(err.Error(), "nil context") {
		t.Fatalf("OpenStore() error = %v, want nil context error", err)
	}
	if store != nil {
		defer store.Close()
		t.Fatalf("OpenStore() store = %#v, want nil on error", store)
	}
	if _, statErr := os.Stat(filepath.Join(runRoot, run.IndexDirName)); !errors.Is(statErr, os.ErrNotExist) {
		t.Fatalf("cache directory stat error = %v, want os.ErrNotExist", statErr)
	}
}

func TestOpenStoreRejectsEmptyRunDirWithoutCreatingRelativeCache(t *testing.T) {
	originalWorkingDir, err := os.Getwd()
	if err != nil {
		t.Fatalf("get working directory: %v", err)
	}
	tempWorkingDir := t.TempDir()
	if err := os.Chdir(tempWorkingDir); err != nil {
		t.Fatalf("chdir to temp working directory: %v", err)
	}
	t.Cleanup(func() {
		if err := os.Chdir(originalWorkingDir); err != nil {
			t.Errorf("restore working directory: %v", err)
		}
	})

	store, err := OpenStore(context.Background(), domain.RunDir{})
	if err == nil {
		t.Fatalf("OpenStore() returned nil error for empty run dir")
	}
	if !strings.Contains(err.Error(), "empty") || !strings.Contains(err.Error(), "run dir") {
		t.Fatalf("OpenStore() error = %v, want empty run dir error", err)
	}
	if store != nil {
		defer store.Close()
		t.Fatalf("OpenStore() store = %#v, want nil on error", store)
	}
	if _, statErr := os.Stat(filepath.Join(tempWorkingDir, run.IndexDirName)); !errors.Is(statErr, os.ErrNotExist) {
		t.Fatalf("relative cache directory stat error = %v, want os.ErrNotExist", statErr)
	}
}

func TestOpenStoreCloseClosesDatabase(t *testing.T) {
	store, err := OpenStore(context.Background(), domain.RunDir{Path: t.TempDir()})
	if err != nil {
		t.Fatalf("OpenStore() returned error: %v", err)
	}

	if err := store.Close(); err != nil {
		t.Fatalf("Close() returned error: %v", err)
	}
	if err := store.DB.PingContext(context.Background()); err == nil {
		t.Fatalf("PingContext() after Close() returned nil error")
	}
	if err := store.Close(); err != nil {
		t.Fatalf("second Close() returned error: %v", err)
	}
}

func TestInitializeSchemaCreatesVersionedMetadataAndTaskTables(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	var version int
	if err := store.DB.QueryRow("PRAGMA user_version").Scan(&version); err != nil {
		t.Fatalf("query user_version: %v", err)
	}
	if version != SchemaVersion {
		t.Fatalf("schema user_version = %d, want %d", version, SchemaVersion)
	}

	requireTableColumns(t, store.DB, "index_metadata", []string{
		"id",
		"schema_version",
		"run_dir",
		"index_path",
		"mode",
		"built_at",
		"freshness",
		"stale_reason",
		"task_count",
	})
	requireTableColumns(t, store.DB, "source_fingerprints", []string{
		"kind",
		"path",
		"mod_time",
		"size",
	})
	requireTableColumns(t, store.DB, "tasks", []string{
		"row_order",
		"id",
		"status",
		"process",
		"name",
		"tag",
		"workdir",
		"exit_code",
		"duration",
		"realtime",
		"cpus",
		"memory",
		"error_summary",
	})
}

func TestInitializeSchemaSupportsNullableNormalizedTaskRows(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	_, err := store.DB.Exec(`
		INSERT INTO tasks (
			row_order,
			id,
			status,
			process,
			name,
			tag,
			workdir,
			exit_code,
			duration,
			realtime,
			cpus,
			memory,
			error_summary
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
	`, int64(7), "ab/c123def", "FAILED", "CALL_VARIANTS", "sample-1", "tumor", "/run/work/ab/c123def", nil, "1m", "60s", "2", "4 GB", "process failed")
	if err != nil {
		t.Fatalf("insert normalized task row: %v", err)
	}

	var rowOrder int64
	var id string
	var exitCode sql.NullInt64
	if err := store.DB.QueryRow("SELECT row_order, id, exit_code FROM tasks WHERE id = ?", "ab/c123def").Scan(&rowOrder, &id, &exitCode); err != nil {
		t.Fatalf("query task row: %v", err)
	}
	if rowOrder != 7 || id != "ab/c123def" {
		t.Fatalf("task identity = (%d, %q), want (7, %q)", rowOrder, id, "ab/c123def")
	}
	if exitCode.Valid {
		t.Fatalf("exit_code.Valid = true, want false for an unknown/blank exit code")
	}
}

func TestInitializeSchemaIsIdempotentAndPreservesRows(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("first InitializeSchema() returned error: %v", err)
	}
	_, err := store.DB.Exec(`
		INSERT INTO index_metadata (
			id,
			schema_version,
			run_dir,
			index_path,
			mode,
			built_at,
			freshness,
			stale_reason,
			task_count
		) VALUES (1, ?, ?, ?, ?, ?, ?, ?, ?)
	`, SchemaVersion, "/runs/example", store.Path, "trace-backed", "2024-01-02T03:04:05Z", "fresh", "", 1)
	if err != nil {
		t.Fatalf("insert metadata row: %v", err)
	}
	_, err = store.DB.Exec(`
		INSERT INTO tasks (
			row_order,
			id,
			status,
			process,
			name,
			tag,
			workdir,
			exit_code,
			duration,
			realtime,
			cpus,
			memory,
			error_summary
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
	`, int64(1), "aa/bbbbb", "COMPLETED", "ALIGN", "sample-2", "normal", "/run/work/aa/bbbbb", 0, "2m", "120s", "1", "2 GB", "")
	if err != nil {
		t.Fatalf("insert task row: %v", err)
	}

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("second InitializeSchema() returned error: %v", err)
	}

	var taskCount int
	if err := store.DB.QueryRow("SELECT COUNT(*) FROM tasks").Scan(&taskCount); err != nil {
		t.Fatalf("count tasks: %v", err)
	}
	if taskCount != 1 {
		t.Fatalf("task row count after second InitializeSchema() = %d, want 1", taskCount)
	}

	var runDir string
	if err := store.DB.QueryRow("SELECT run_dir FROM index_metadata WHERE id = 1").Scan(&runDir); err != nil {
		t.Fatalf("query metadata row: %v", err)
	}
	if runDir != "/runs/example" {
		t.Fatalf("metadata run_dir after second InitializeSchema() = %q, want %q", runDir, "/runs/example")
	}
}

func TestInitializeSchemaRejectsNilStoreOrDatabase(t *testing.T) {
	if err := InitializeSchema(context.Background(), nil); err == nil {
		t.Fatalf("InitializeSchema(nil) returned nil error")
	}
	if err := InitializeSchema(context.Background(), &Store{}); err == nil {
		t.Fatalf("InitializeSchema(store with nil DB) returned nil error")
	}
}

func TestReadMetadataReturnsNilWhenMetadataRowMissing(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	metadata, err := ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata() returned error: %v", err)
	}
	if metadata != nil {
		t.Fatalf("ReadMetadata() = %#v, want nil metadata for an empty index", metadata)
	}
}

func TestReadMetadataLoadsMetadataAndFingerprintsByKind(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	builtAt := time.Date(2024, time.January, 2, 3, 4, 5, 6789, time.UTC)
	traceModTime := time.Date(2024, time.February, 3, 4, 5, 6, 7890, time.UTC)
	logModTime := time.Date(2024, time.March, 4, 5, 6, 7, 8901, time.UTC)
	_, err := store.DB.Exec(`
		INSERT INTO index_metadata (
			id,
			schema_version,
			run_dir,
			index_path,
			mode,
			built_at,
			freshness,
			stale_reason,
			task_count
		) VALUES (1, ?, ?, ?, ?, ?, ?, ?, ?)
	`, SchemaVersion, "/runs/example", store.Path, string(domain.IndexModeTraceBacked), builtAt.Format(time.RFC3339Nano), string(domain.IndexFreshnessFresh), "", 42)
	if err != nil {
		t.Fatalf("insert metadata row: %v", err)
	}
	_, err = store.DB.Exec(`
		INSERT INTO source_fingerprints (kind, path, mod_time, size)
		VALUES (?, ?, ?, ?), (?, ?, ?, ?)
	`,
		string(domain.SourceKindLog), "/runs/example/.nextflow.log", logModTime.Format(time.RFC3339Nano), int64(1234),
		string(domain.SourceKindTrace), "/runs/example/trace.txt", traceModTime.Format(time.RFC3339Nano), int64(5678),
	)
	if err != nil {
		t.Fatalf("insert source fingerprints: %v", err)
	}

	metadata, err := ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata() returned error: %v", err)
	}
	if metadata == nil {
		t.Fatalf("ReadMetadata() = nil, want metadata")
	}
	if metadata.SchemaVersion != SchemaVersion {
		t.Fatalf("SchemaVersion = %d, want %d", metadata.SchemaVersion, SchemaVersion)
	}
	if metadata.RunDir != "/runs/example" {
		t.Fatalf("RunDir = %q, want %q", metadata.RunDir, "/runs/example")
	}
	if metadata.IndexPath != store.Path {
		t.Fatalf("IndexPath = %q, want %q", metadata.IndexPath, store.Path)
	}
	if metadata.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("Mode = %q, want %q", metadata.Mode, domain.IndexModeTraceBacked)
	}
	if !metadata.BuiltAt.Equal(builtAt) {
		t.Fatalf("BuiltAt = %s, want %s", metadata.BuiltAt.Format(time.RFC3339Nano), builtAt.Format(time.RFC3339Nano))
	}
	if metadata.Freshness != domain.IndexFreshnessFresh {
		t.Fatalf("Freshness = %q, want %q", metadata.Freshness, domain.IndexFreshnessFresh)
	}
	if metadata.StaleReason != "" {
		t.Fatalf("StaleReason = %q, want empty", metadata.StaleReason)
	}
	if metadata.TaskCount != 42 {
		t.Fatalf("TaskCount = %d, want 42", metadata.TaskCount)
	}
	assertSourceFingerprint(t, "Trace", metadata.Trace, domain.SourceFingerprint{
		Kind:    domain.SourceKindTrace,
		Path:    "/runs/example/trace.txt",
		ModTime: traceModTime,
		Size:    5678,
	})
	assertSourceFingerprint(t, "Log", metadata.Log, domain.SourceFingerprint{
		Kind:    domain.SourceKindLog,
		Path:    "/runs/example/.nextflow.log",
		ModTime: logModTime,
		Size:    1234,
	})
}

func TestReadMetadataRejectsNilStoreOrDatabase(t *testing.T) {
	if _, err := ReadMetadata(context.Background(), nil); err == nil {
		t.Fatalf("ReadMetadata(nil) returned nil error")
	}
	if _, err := ReadMetadata(context.Background(), &Store{}); err == nil {
		t.Fatalf("ReadMetadata(store with nil DB) returned nil error")
	}
}

func TestWriteMetadataRoundTripsMetadataAndFingerprintsWithReadMetadata(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	builtAt := time.Date(2024, time.April, 5, 6, 7, 8, 9012, time.UTC)
	traceModTime := time.Date(2024, time.May, 6, 7, 8, 9, 1234, time.UTC)
	logModTime := time.Date(2024, time.June, 7, 8, 9, 10, 2345, time.UTC)
	want := domain.IndexMetadata{
		SchemaVersion: SchemaVersion,
		RunDir:        "/runs/round-trip",
		IndexPath:     store.Path,
		Mode:          domain.IndexModeTraceBacked,
		Trace: &domain.SourceFingerprint{
			Kind:    domain.SourceKindTrace,
			Path:    "/runs/round-trip/trace.txt",
			ModTime: traceModTime,
			Size:    5678,
		},
		Log: &domain.SourceFingerprint{
			Kind:    domain.SourceKindLog,
			Path:    "/runs/round-trip/.nextflow.log",
			ModTime: logModTime,
			Size:    1234,
		},
		BuiltAt:     builtAt,
		Freshness:   domain.IndexFreshnessFresh,
		StaleReason: "",
		TaskCount:   42,
	}

	if err := WriteMetadata(context.Background(), store, want); err != nil {
		t.Fatalf("WriteMetadata() returned error: %v", err)
	}

	got, err := ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata() returned error: %v", err)
	}
	if got == nil {
		t.Fatalf("ReadMetadata() = nil, want metadata")
	}
	if got.SchemaVersion != want.SchemaVersion {
		t.Fatalf("SchemaVersion = %d, want %d", got.SchemaVersion, want.SchemaVersion)
	}
	if got.RunDir != want.RunDir {
		t.Fatalf("RunDir = %q, want %q", got.RunDir, want.RunDir)
	}
	if got.IndexPath != want.IndexPath {
		t.Fatalf("IndexPath = %q, want %q", got.IndexPath, want.IndexPath)
	}
	if got.Mode != want.Mode {
		t.Fatalf("Mode = %q, want %q", got.Mode, want.Mode)
	}
	if !got.BuiltAt.Equal(want.BuiltAt) {
		t.Fatalf("BuiltAt = %s, want %s", got.BuiltAt.Format(time.RFC3339Nano), want.BuiltAt.Format(time.RFC3339Nano))
	}
	if got.Freshness != want.Freshness {
		t.Fatalf("Freshness = %q, want %q", got.Freshness, want.Freshness)
	}
	if got.StaleReason != want.StaleReason {
		t.Fatalf("StaleReason = %q, want %q", got.StaleReason, want.StaleReason)
	}
	if got.TaskCount != want.TaskCount {
		t.Fatalf("TaskCount = %d, want %d", got.TaskCount, want.TaskCount)
	}
	assertSourceFingerprint(t, "Trace", got.Trace, *want.Trace)
	assertSourceFingerprint(t, "Log", got.Log, *want.Log)
}

func TestWriteMetadataReplacesSingletonMetadataAndFingerprints(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	firstBuiltAt := time.Date(2024, time.January, 2, 3, 4, 5, 0, time.UTC)
	firstTraceModTime := time.Date(2024, time.January, 3, 4, 5, 6, 0, time.UTC)
	firstLogModTime := time.Date(2024, time.January, 4, 5, 6, 7, 0, time.UTC)
	first := domain.IndexMetadata{
		SchemaVersion: SchemaVersion,
		RunDir:        "/runs/first",
		IndexPath:     store.Path,
		Mode:          domain.IndexModeTraceBacked,
		Trace: &domain.SourceFingerprint{
			Kind:    domain.SourceKindTrace,
			Path:    "/runs/first/trace.txt",
			ModTime: firstTraceModTime,
			Size:    100,
		},
		Log: &domain.SourceFingerprint{
			Kind:    domain.SourceKindLog,
			Path:    "/runs/first/.nextflow.log",
			ModTime: firstLogModTime,
			Size:    200,
		},
		BuiltAt:     firstBuiltAt,
		Freshness:   domain.IndexFreshnessFresh,
		StaleReason: "",
		TaskCount:   10,
	}
	if err := WriteMetadata(context.Background(), store, first); err != nil {
		t.Fatalf("first WriteMetadata() returned error: %v", err)
	}

	secondBuiltAt := time.Date(2024, time.February, 2, 3, 4, 5, 0, time.UTC)
	secondLogModTime := time.Date(2024, time.February, 4, 5, 6, 7, 0, time.UTC)
	second := domain.IndexMetadata{
		SchemaVersion: SchemaVersion,
		RunDir:        "/runs/second",
		IndexPath:     store.Path,
		Mode:          domain.IndexModeLogOnly,
		Trace:         nil,
		Log: &domain.SourceFingerprint{
			Kind:    domain.SourceKindLog,
			Path:    "/runs/second/.nextflow.log",
			ModTime: secondLogModTime,
			Size:    300,
		},
		BuiltAt:     secondBuiltAt,
		Freshness:   domain.IndexFreshnessStale,
		StaleReason: "trace artifact removed",
		TaskCount:   0,
	}
	if err := WriteMetadata(context.Background(), store, second); err != nil {
		t.Fatalf("second WriteMetadata() returned error: %v", err)
	}

	got, err := ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata() returned error: %v", err)
	}
	if got == nil {
		t.Fatalf("ReadMetadata() = nil, want metadata")
	}
	if got.RunDir != second.RunDir {
		t.Fatalf("RunDir = %q, want %q", got.RunDir, second.RunDir)
	}
	if got.Mode != second.Mode {
		t.Fatalf("Mode = %q, want %q", got.Mode, second.Mode)
	}
	if got.Trace != nil {
		t.Fatalf("Trace = %#v, want nil after replacement without a trace fingerprint", got.Trace)
	}
	assertSourceFingerprint(t, "Log", got.Log, *second.Log)
	if !got.BuiltAt.Equal(second.BuiltAt) {
		t.Fatalf("BuiltAt = %s, want %s", got.BuiltAt.Format(time.RFC3339Nano), second.BuiltAt.Format(time.RFC3339Nano))
	}
	if got.Freshness != second.Freshness {
		t.Fatalf("Freshness = %q, want %q", got.Freshness, second.Freshness)
	}
	if got.StaleReason != second.StaleReason {
		t.Fatalf("StaleReason = %q, want %q", got.StaleReason, second.StaleReason)
	}
	if got.TaskCount != second.TaskCount {
		t.Fatalf("TaskCount = %d, want %d", got.TaskCount, second.TaskCount)
	}

	var metadataRows int
	if err := store.DB.QueryRow("SELECT COUNT(*) FROM index_metadata").Scan(&metadataRows); err != nil {
		t.Fatalf("count metadata rows: %v", err)
	}
	if metadataRows != 1 {
		t.Fatalf("metadata row count = %d, want 1", metadataRows)
	}

	var fingerprintRows int
	if err := store.DB.QueryRow("SELECT COUNT(*) FROM source_fingerprints WHERE kind IN (?, ?)", string(domain.SourceKindTrace), string(domain.SourceKindLog)).Scan(&fingerprintRows); err != nil {
		t.Fatalf("count source fingerprints: %v", err)
	}
	if fingerprintRows != 1 {
		t.Fatalf("trace/log fingerprint row count = %d, want 1", fingerprintRows)
	}
}

func TestWriteMetadataRollsBackMetadataWhenFingerprintWriteFails(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	originalBuiltAt := time.Date(2024, time.January, 2, 3, 4, 5, 0, time.UTC)
	originalTraceModTime := time.Date(2024, time.January, 3, 4, 5, 6, 0, time.UTC)
	_, err := store.DB.Exec(`
		INSERT INTO index_metadata (
			id,
			schema_version,
			run_dir,
			index_path,
			mode,
			built_at,
			freshness,
			stale_reason,
			task_count
		) VALUES (1, ?, ?, ?, ?, ?, ?, ?, ?)
	`, SchemaVersion, "/runs/original", store.Path, string(domain.IndexModeTraceBacked), originalBuiltAt.Format(time.RFC3339Nano), string(domain.IndexFreshnessFresh), "", 12)
	if err != nil {
		t.Fatalf("insert original metadata row: %v", err)
	}
	_, err = store.DB.Exec(`
		INSERT INTO source_fingerprints (kind, path, mod_time, size)
		VALUES (?, ?, ?, ?)
	`, string(domain.SourceKindTrace), "/runs/original/trace.txt", originalTraceModTime.Format(time.RFC3339Nano), int64(123))
	if err != nil {
		t.Fatalf("insert original source fingerprint: %v", err)
	}
	_, err = store.DB.Exec(`
		CREATE TRIGGER fail_source_fingerprint_insert
		BEFORE INSERT ON source_fingerprints
		BEGIN
			SELECT RAISE(ABORT, 'source fingerprint insert blocked');
		END
	`)
	if err != nil {
		t.Fatalf("create failing source fingerprint trigger: %v", err)
	}

	newBuiltAt := time.Date(2024, time.March, 2, 3, 4, 5, 0, time.UTC)
	newTraceModTime := time.Date(2024, time.March, 3, 4, 5, 6, 0, time.UTC)
	newMetadata := domain.IndexMetadata{
		SchemaVersion: SchemaVersion,
		RunDir:        "/runs/new",
		IndexPath:     store.Path,
		Mode:          domain.IndexModeTraceBacked,
		Trace: &domain.SourceFingerprint{
			Kind:    domain.SourceKindTrace,
			Path:    "/runs/new/trace.txt",
			ModTime: newTraceModTime,
			Size:    456,
		},
		BuiltAt:   newBuiltAt,
		Freshness: domain.IndexFreshnessFresh,
		TaskCount: 99,
	}
	if err := WriteMetadata(context.Background(), store, newMetadata); err == nil {
		t.Fatalf("WriteMetadata() returned nil error with failing source fingerprint trigger")
	}

	got, err := ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata() returned error after failed WriteMetadata(): %v", err)
	}
	if got == nil {
		t.Fatalf("ReadMetadata() = nil, want original metadata after failed WriteMetadata()")
	}
	if got.RunDir != "/runs/original" {
		t.Fatalf("RunDir after failed WriteMetadata() = %q, want %q", got.RunDir, "/runs/original")
	}
	if got.TaskCount != 12 {
		t.Fatalf("TaskCount after failed WriteMetadata() = %d, want 12", got.TaskCount)
	}
	assertSourceFingerprint(t, "Trace after failed WriteMetadata()", got.Trace, domain.SourceFingerprint{
		Kind:    domain.SourceKindTrace,
		Path:    "/runs/original/trace.txt",
		ModTime: originalTraceModTime,
		Size:    123,
	})
}

func TestWriteMetadataRejectsNilStoreOrDatabase(t *testing.T) {
	if err := WriteMetadata(context.Background(), nil, domain.IndexMetadata{}); err == nil {
		t.Fatalf("WriteMetadata(nil) returned nil error")
	}
	if err := WriteMetadata(context.Background(), &Store{}, domain.IndexMetadata{}); err == nil {
		t.Fatalf("WriteMetadata(store with nil DB) returned nil error")
	}
}

func TestCheckFreshnessReportsUnsupportedArtifacts(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	freshness, reason, err := CheckFreshness(context.Background(), store, domain.ArtifactSet{
		RunDir: domain.RunDir{Path: "/runs/empty"},
		Mode:   domain.IndexModeUnsupported,
	})
	if err != nil {
		t.Fatalf("CheckFreshness(unsupported artifacts) returned error: %v", err)
	}
	if freshness != domain.IndexFreshnessUnsupported {
		t.Fatalf("freshness = %q, want %q", freshness, domain.IndexFreshnessUnsupported)
	}
	if reason != "no supported artifacts" {
		t.Fatalf("reason = %q, want %q", reason, "no supported artifacts")
	}
}

func TestCheckFreshnessReportsMissingMetadata(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	traceModTime := time.Date(2024, time.January, 2, 3, 4, 5, 0, time.UTC)
	freshness, reason, err := CheckFreshness(context.Background(), store, domain.ArtifactSet{
		RunDir: domain.RunDir{Path: "/runs/missing"},
		Mode:   domain.IndexModeTraceBacked,
		Trace: &domain.SourceFingerprint{
			Kind:    domain.SourceKindTrace,
			Path:    "/runs/missing/trace.txt",
			ModTime: traceModTime,
			Size:    100,
		},
	})
	if err != nil {
		t.Fatalf("CheckFreshness(missing metadata) returned error: %v", err)
	}
	if freshness != domain.IndexFreshnessMissing {
		t.Fatalf("freshness = %q, want %q", freshness, domain.IndexFreshnessMissing)
	}
	if reason != "index metadata missing" {
		t.Fatalf("reason = %q, want %q", reason, "index metadata missing")
	}
}

func TestCheckFreshnessReportsFreshWhenModeAndFingerprintsMatch(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	builtAt := time.Date(2024, time.April, 5, 6, 7, 8, 0, time.UTC)
	traceModTime := time.Date(2024, time.April, 5, 7, 8, 9, 0, time.UTC)
	logModTime := time.Date(2024, time.April, 5, 8, 9, 10, 0, time.UTC)
	traceSource := domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: "/runs/fresh/trace.txt", ModTime: traceModTime, Size: 123}
	logSource := domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/fresh/.nextflow.log", ModTime: logModTime, Size: 456}

	if err := WriteMetadata(context.Background(), store, domain.IndexMetadata{
		SchemaVersion: SchemaVersion,
		RunDir:        "/runs/fresh",
		IndexPath:     store.Path,
		Mode:          domain.IndexModeTraceBacked,
		Trace:         &traceSource,
		Log:           &logSource,
		BuiltAt:       builtAt,
		Freshness:     domain.IndexFreshnessFresh,
		TaskCount:     12,
	}); err != nil {
		t.Fatalf("WriteMetadata() returned error: %v", err)
	}

	freshness, reason, err := CheckFreshness(context.Background(), store, domain.ArtifactSet{
		RunDir: domain.RunDir{Path: "/runs/fresh"},
		Mode:   domain.IndexModeTraceBacked,
		Trace:  &traceSource,
		Log:    &logSource,
	})
	if err != nil {
		t.Fatalf("CheckFreshness(fresh artifacts) returned error: %v", err)
	}
	if freshness != domain.IndexFreshnessFresh {
		t.Fatalf("freshness = %q, want %q", freshness, domain.IndexFreshnessFresh)
	}
	if reason != "" {
		t.Fatalf("reason = %q, want empty", reason)
	}
}

func TestCheckFreshnessReportsStaleWhenModeChanges(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	builtAt := time.Date(2024, time.May, 6, 7, 8, 9, 0, time.UTC)
	traceModTime := time.Date(2024, time.May, 6, 8, 9, 10, 0, time.UTC)
	logModTime := time.Date(2024, time.May, 6, 9, 10, 11, 0, time.UTC)
	traceSource := domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: "/runs/mode/trace.txt", ModTime: traceModTime, Size: 123}
	logSource := domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/mode/.nextflow.log", ModTime: logModTime, Size: 456}

	if err := WriteMetadata(context.Background(), store, domain.IndexMetadata{
		SchemaVersion: SchemaVersion,
		RunDir:        "/runs/mode",
		IndexPath:     store.Path,
		Mode:          domain.IndexModeTraceBacked,
		Trace:         &traceSource,
		Log:           &logSource,
		BuiltAt:       builtAt,
		Freshness:     domain.IndexFreshnessFresh,
		TaskCount:     12,
	}); err != nil {
		t.Fatalf("WriteMetadata() returned error: %v", err)
	}

	freshness, reason, err := CheckFreshness(context.Background(), store, domain.ArtifactSet{
		RunDir: domain.RunDir{Path: "/runs/mode"},
		Mode:   domain.IndexModeLogOnly,
		Log:    &logSource,
	})
	if err != nil {
		t.Fatalf("CheckFreshness(mode changed) returned error: %v", err)
	}
	if freshness != domain.IndexFreshnessStale {
		t.Fatalf("freshness = %q, want %q", freshness, domain.IndexFreshnessStale)
	}
	if reason != "index mode changed" {
		t.Fatalf("reason = %q, want %q", reason, "index mode changed")
	}
}

func TestCheckFreshnessReportsStaleWhenSelectedSourceFingerprintChanges(t *testing.T) {
	baseTraceModTime := time.Date(2024, time.June, 7, 8, 9, 10, 0, time.UTC)
	baseLogModTime := time.Date(2024, time.June, 7, 9, 10, 11, 0, time.UTC)
	baseTrace := domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: "/runs/sources/trace.txt", ModTime: baseTraceModTime, Size: 100}
	baseLog := domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/sources/.nextflow.log", ModTime: baseLogModTime, Size: 200}

	tracePathChanged := baseTrace
	tracePathChanged.Path = "/runs/sources/new-trace.txt"
	traceModTimeChanged := baseTrace
	traceModTimeChanged.ModTime = baseTraceModTime.Add(time.Second)
	traceSizeChanged := baseTrace
	traceSizeChanged.Size = 101
	logPathChanged := baseLog
	logPathChanged.Path = "/runs/sources/.nextflow.new.log"
	logModTimeChanged := baseLog
	logModTimeChanged.ModTime = baseLogModTime.Add(time.Second)
	logSizeChanged := baseLog
	logSizeChanged.Size = 201

	cases := []struct {
		name       string
		trace      domain.SourceFingerprint
		log        domain.SourceFingerprint
		wantReason string
	}{
		{name: "trace path changed", trace: tracePathChanged, log: baseLog, wantReason: "selected trace changed"},
		{name: "trace mtime changed", trace: traceModTimeChanged, log: baseLog, wantReason: "selected trace changed"},
		{name: "trace size changed", trace: traceSizeChanged, log: baseLog, wantReason: "selected trace changed"},
		{name: "log path changed", trace: baseTrace, log: logPathChanged, wantReason: "selected log changed"},
		{name: "log mtime changed", trace: baseTrace, log: logModTimeChanged, wantReason: "selected log changed"},
		{name: "log size changed", trace: baseTrace, log: logSizeChanged, wantReason: "selected log changed"},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			store, cleanup := openTempSQLiteStore(t)
			defer cleanup()

			if err := InitializeSchema(context.Background(), store); err != nil {
				t.Fatalf("InitializeSchema() returned error: %v", err)
			}

			builtAt := time.Date(2024, time.June, 7, 7, 8, 9, 0, time.UTC)
			if err := WriteMetadata(context.Background(), store, domain.IndexMetadata{
				SchemaVersion: SchemaVersion,
				RunDir:        "/runs/sources",
				IndexPath:     store.Path,
				Mode:          domain.IndexModeTraceBacked,
				Trace:         &baseTrace,
				Log:           &baseLog,
				BuiltAt:       builtAt,
				Freshness:     domain.IndexFreshnessFresh,
				TaskCount:     12,
			}); err != nil {
				t.Fatalf("WriteMetadata() returned error: %v", err)
			}

			traceSource := tc.trace
			logSource := tc.log
			freshness, reason, err := CheckFreshness(context.Background(), store, domain.ArtifactSet{
				RunDir: domain.RunDir{Path: "/runs/sources"},
				Mode:   domain.IndexModeTraceBacked,
				Trace:  &traceSource,
				Log:    &logSource,
			})
			if err != nil {
				t.Fatalf("CheckFreshness(source changed) returned error: %v", err)
			}
			if freshness != domain.IndexFreshnessStale {
				t.Fatalf("freshness = %q, want %q", freshness, domain.IndexFreshnessStale)
			}
			if reason != tc.wantReason {
				t.Fatalf("reason = %q, want %q", reason, tc.wantReason)
			}
		})
	}
}

func TestCheckFreshnessPropagatesReadMetadataErrors(t *testing.T) {
	freshness, reason, err := CheckFreshness(context.Background(), nil, domain.ArtifactSet{
		Mode: domain.IndexModeTraceBacked,
	})
	if err == nil {
		t.Fatalf("CheckFreshness(nil store) returned nil error")
	}
	if freshness != domain.IndexFreshnessUnknown {
		t.Fatalf("freshness = %q, want %q", freshness, domain.IndexFreshnessUnknown)
	}
	if reason != "" {
		t.Fatalf("reason = %q, want empty on error", reason)
	}
	if !strings.Contains(err.Error(), "read metadata") {
		t.Fatalf("error = %q, want it to mention read metadata", err.Error())
	}
}

func TestRefreshMetadataCreatesSchemaAndFreshTraceBackedMetadata(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	runDir := domain.RunDir{Path: t.TempDir()}
	traceModTime := time.Date(2024, time.January, 2, 3, 4, 5, 0, time.UTC)
	logModTime := time.Date(2024, time.January, 2, 4, 5, 6, 0, time.UTC)
	traceSource := domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: filepath.Join(runDir.Path, "trace.txt"), ModTime: traceModTime, Size: 1234}
	logSource := domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: filepath.Join(runDir.Path, ".nextflow.log"), ModTime: logModTime, Size: 5678}
	artifacts := domain.ArtifactSet{
		RunDir: runDir,
		Mode:   domain.IndexModeTraceBacked,
		Trace:  &traceSource,
		Log:    &logSource,
	}

	startedAt := time.Now()
	metadata, err := RefreshMetadata(context.Background(), store, runDir, artifacts)
	finishedAt := time.Now()
	if err != nil {
		t.Fatalf("RefreshMetadata(trace-backed) returned error: %v", err)
	}

	if metadata.SchemaVersion != SchemaVersion {
		t.Fatalf("SchemaVersion = %d, want %d", metadata.SchemaVersion, SchemaVersion)
	}
	if metadata.RunDir != runDir.Path {
		t.Fatalf("RunDir = %q, want %q", metadata.RunDir, runDir.Path)
	}
	if metadata.IndexPath != store.Path {
		t.Fatalf("IndexPath = %q, want %q", metadata.IndexPath, store.Path)
	}
	if metadata.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("Mode = %q, want %q", metadata.Mode, domain.IndexModeTraceBacked)
	}
	if metadata.Freshness != domain.IndexFreshnessFresh {
		t.Fatalf("Freshness = %q, want %q", metadata.Freshness, domain.IndexFreshnessFresh)
	}
	if metadata.StaleReason != "" {
		t.Fatalf("StaleReason = %q, want empty for a fresh metadata rebuild", metadata.StaleReason)
	}
	if metadata.TaskCount != 0 {
		t.Fatalf("TaskCount = %d, want 0 for metadata-only refresh", metadata.TaskCount)
	}
	if metadata.BuiltAt.IsZero() || metadata.BuiltAt.Before(startedAt) || metadata.BuiltAt.After(finishedAt) {
		t.Fatalf("BuiltAt = %s, want non-zero timestamp between %s and %s", metadata.BuiltAt.Format(time.RFC3339Nano), startedAt.Format(time.RFC3339Nano), finishedAt.Format(time.RFC3339Nano))
	}
	assertSourceFingerprint(t, "Trace", metadata.Trace, traceSource)
	assertSourceFingerprint(t, "Log", metadata.Log, logSource)

	var version int
	if err := store.DB.QueryRow("PRAGMA user_version").Scan(&version); err != nil {
		t.Fatalf("query user_version after RefreshMetadata(): %v", err)
	}
	if version != SchemaVersion {
		t.Fatalf("schema user_version = %d, want %d", version, SchemaVersion)
	}
	var taskRows int
	if err := store.DB.QueryRow("SELECT COUNT(*) FROM tasks").Scan(&taskRows); err != nil {
		t.Fatalf("count task rows after metadata-only refresh: %v", err)
	}
	if taskRows != 0 {
		t.Fatalf("task row count = %d, want 0 for metadata-only refresh", taskRows)
	}

	persisted, err := ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata() after RefreshMetadata() returned error: %v", err)
	}
	if persisted == nil {
		t.Fatalf("ReadMetadata() after RefreshMetadata() = nil, want metadata")
	}
	if !reflect.DeepEqual(*persisted, metadata) {
		t.Fatalf("persisted metadata = %#v, want returned metadata %#v", *persisted, metadata)
	}
}

func TestRefreshMetadataStoresFreshLogOnlyMetadata(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	runDir := domain.RunDir{Path: t.TempDir()}
	logModTime := time.Date(2024, time.February, 3, 4, 5, 6, 0, time.UTC)
	logSource := domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: filepath.Join(runDir.Path, ".nextflow.log"), ModTime: logModTime, Size: 2468}

	metadata, err := RefreshMetadata(context.Background(), store, runDir, domain.ArtifactSet{
		RunDir: runDir,
		Mode:   domain.IndexModeLogOnly,
		Log:    &logSource,
	})
	if err != nil {
		t.Fatalf("RefreshMetadata(log-only) returned error: %v", err)
	}

	if metadata.Mode != domain.IndexModeLogOnly {
		t.Fatalf("Mode = %q, want %q", metadata.Mode, domain.IndexModeLogOnly)
	}
	if metadata.Freshness != domain.IndexFreshnessFresh {
		t.Fatalf("Freshness = %q, want %q for fresh log-only metadata", metadata.Freshness, domain.IndexFreshnessFresh)
	}
	if metadata.StaleReason != "" {
		t.Fatalf("StaleReason = %q, want empty for fresh log-only metadata", metadata.StaleReason)
	}
	if metadata.Trace != nil {
		t.Fatalf("Trace = %#v, want nil for log-only metadata", metadata.Trace)
	}
	assertSourceFingerprint(t, "Log", metadata.Log, logSource)
	if metadata.TaskCount != 0 {
		t.Fatalf("TaskCount = %d, want 0 for metadata-only log-only refresh", metadata.TaskCount)
	}

	freshness, reason, err := CheckFreshness(context.Background(), store, domain.ArtifactSet{RunDir: runDir, Mode: domain.IndexModeLogOnly, Log: &logSource})
	if err != nil {
		t.Fatalf("CheckFreshness() after log-only RefreshMetadata() returned error: %v", err)
	}
	if freshness != domain.IndexFreshnessFresh || reason != "" {
		t.Fatalf("CheckFreshness() = (%q, %q), want (%q, empty)", freshness, reason, domain.IndexFreshnessFresh)
	}
}

func TestRefreshMetadataMarksUnsupportedArtifactsUnsupported(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	runDir := domain.RunDir{Path: t.TempDir()}
	metadata, err := RefreshMetadata(context.Background(), store, runDir, domain.ArtifactSet{
		RunDir: runDir,
		Mode:   domain.IndexModeUnsupported,
	})
	if err != nil {
		t.Fatalf("RefreshMetadata(unsupported) returned error: %v", err)
	}

	if metadata.Mode != domain.IndexModeUnsupported {
		t.Fatalf("Mode = %q, want %q", metadata.Mode, domain.IndexModeUnsupported)
	}
	if metadata.Freshness != domain.IndexFreshnessUnsupported {
		t.Fatalf("Freshness = %q, want %q", metadata.Freshness, domain.IndexFreshnessUnsupported)
	}
	if metadata.StaleReason != "no supported artifacts" {
		t.Fatalf("StaleReason = %q, want %q", metadata.StaleReason, "no supported artifacts")
	}
	if metadata.Trace != nil || metadata.Log != nil {
		t.Fatalf("Trace/Log = %#v/%#v, want nil fingerprints for unsupported artifacts", metadata.Trace, metadata.Log)
	}
	if metadata.TaskCount != 0 {
		t.Fatalf("TaskCount = %d, want 0 for unsupported metadata", metadata.TaskCount)
	}
	if metadata.RunDir != runDir.Path {
		t.Fatalf("RunDir = %q, want %q", metadata.RunDir, runDir.Path)
	}
	if metadata.IndexPath != store.Path {
		t.Fatalf("IndexPath = %q, want %q", metadata.IndexPath, store.Path)
	}

	persisted, err := ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata() after unsupported RefreshMetadata() returned error: %v", err)
	}
	if persisted == nil {
		t.Fatalf("ReadMetadata() after unsupported RefreshMetadata() = nil, want metadata")
	}
	if persisted.Freshness != domain.IndexFreshnessUnsupported || persisted.StaleReason != "no supported artifacts" {
		t.Fatalf("persisted freshness/reason = %q/%q, want %q/%q", persisted.Freshness, persisted.StaleReason, domain.IndexFreshnessUnsupported, "no supported artifacts")
	}
}

func TestRefreshMetadataRejectsNilContextStoreOrDatabase(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if _, err := RefreshMetadata(nil, store, domain.RunDir{Path: t.TempDir()}, domain.ArtifactSet{Mode: domain.IndexModeTraceBacked}); err == nil {
		t.Fatalf("RefreshMetadata(nil context) returned nil error")
	} else if !strings.Contains(err.Error(), "nil context") {
		t.Fatalf("RefreshMetadata(nil context) error = %v, want nil context error", err)
	}
	if _, err := RefreshMetadata(context.Background(), nil, domain.RunDir{Path: t.TempDir()}, domain.ArtifactSet{Mode: domain.IndexModeTraceBacked}); err == nil {
		t.Fatalf("RefreshMetadata(nil store) returned nil error")
	}
	if _, err := RefreshMetadata(context.Background(), &Store{}, domain.RunDir{Path: t.TempDir()}, domain.ArtifactSet{Mode: domain.IndexModeTraceBacked}); err == nil {
		t.Fatalf("RefreshMetadata(store with nil DB) returned nil error")
	}
}

func TestRebuildTraceIndexParsesSelectedTraceStoresTasksAndMetadata(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	runDir := domain.RunDir{Path: t.TempDir()}
	firstWorkdir := filepath.Join(runDir.Path, "work", "AB", "C123DEF")
	secondWorkdir := filepath.Join(runDir.Path, "work", "DE", "F456")
	tracePath := filepath.Join(runDir.Path, "trace.csv")
	traceContent := strings.Join([]string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AB/C123DEF,COMPLETED,ALIGN,ALIGN (sample-1),sample-1," + firstWorkdir + ",0,1m,60s,2,4 GB",
		"DE/F456,FAILED,QUANT,QUANT (sample-2),sample-2," + secondWorkdir + ",137,2m,120s,4,8 GB",
		"",
	}, "\n")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		t.Fatalf("write temp trace: %v", err)
	}
	traceInfo, err := os.Stat(tracePath)
	if err != nil {
		t.Fatalf("stat temp trace: %v", err)
	}
	traceSource := domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: tracePath, ModTime: traceInfo.ModTime(), Size: traceInfo.Size()}

	logPath := filepath.Join(runDir.Path, ".nextflow.log")
	if err := os.WriteFile(logPath, []byte("log remains source metadata only\n"), 0o644); err != nil {
		t.Fatalf("write temp log: %v", err)
	}
	logInfo, err := os.Stat(logPath)
	if err != nil {
		t.Fatalf("stat temp log: %v", err)
	}
	logSource := domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath, ModTime: logInfo.ModTime(), Size: logInfo.Size()}
	artifacts := domain.ArtifactSet{
		RunDir: runDir,
		Mode:   domain.IndexModeTraceBacked,
		Trace:  &traceSource,
		Log:    &logSource,
	}

	startedAt := time.Now()
	metadata, err := RebuildTraceIndex(context.Background(), store, runDir, artifacts)
	finishedAt := time.Now()
	if err != nil {
		t.Fatalf("RebuildTraceIndex(trace-backed) returned error: %v", err)
	}

	if metadata.SchemaVersion != SchemaVersion {
		t.Fatalf("SchemaVersion = %d, want %d", metadata.SchemaVersion, SchemaVersion)
	}
	if metadata.RunDir != runDir.Path {
		t.Fatalf("RunDir = %q, want %q", metadata.RunDir, runDir.Path)
	}
	if metadata.IndexPath != store.Path {
		t.Fatalf("IndexPath = %q, want %q", metadata.IndexPath, store.Path)
	}
	if metadata.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("Mode = %q, want %q", metadata.Mode, domain.IndexModeTraceBacked)
	}
	if metadata.Freshness != domain.IndexFreshnessFresh {
		t.Fatalf("Freshness = %q, want %q", metadata.Freshness, domain.IndexFreshnessFresh)
	}
	if metadata.StaleReason != "" {
		t.Fatalf("StaleReason = %q, want empty", metadata.StaleReason)
	}
	if metadata.TaskCount != 2 {
		t.Fatalf("TaskCount = %d, want 2 parsed tasks", metadata.TaskCount)
	}
	if metadata.BuiltAt.IsZero() || metadata.BuiltAt.Before(startedAt) || metadata.BuiltAt.After(finishedAt) {
		t.Fatalf("BuiltAt = %s, want non-zero timestamp between %s and %s", metadata.BuiltAt.Format(time.RFC3339Nano), startedAt.Format(time.RFC3339Nano), finishedAt.Format(time.RFC3339Nano))
	}
	assertSourceFingerprint(t, "Trace", metadata.Trace, traceSource)
	assertSourceFingerprint(t, "Log", metadata.Log, logSource)

	firstExit := 0
	secondExit := 137
	wantTasks := []domain.Task{
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
	gotTasks, err := QueryTasks(context.Background(), store, domain.TaskQuery{})
	if err != nil {
		t.Fatalf("QueryTasks() after RebuildTraceIndex() returned error: %v", err)
	}
	if !reflect.DeepEqual(gotTasks, wantTasks) {
		t.Fatalf("tasks after RebuildTraceIndex() = %#v, want %#v", gotTasks, wantTasks)
	}

	persisted, err := ReadMetadata(context.Background(), store)
	if err != nil {
		t.Fatalf("ReadMetadata() after RebuildTraceIndex() returned error: %v", err)
	}
	if persisted == nil {
		t.Fatalf("ReadMetadata() after RebuildTraceIndex() = nil, want metadata")
	}
	if !reflect.DeepEqual(*persisted, metadata) {
		t.Fatalf("persisted metadata = %#v, want returned metadata %#v", *persisted, metadata)
	}

	freshness, reason, err := CheckFreshness(context.Background(), store, artifacts)
	if err != nil {
		t.Fatalf("CheckFreshness() after RebuildTraceIndex() returned error: %v", err)
	}
	if freshness != domain.IndexFreshnessFresh || reason != "" {
		t.Fatalf("CheckFreshness() = (%q, %q), want (%q, empty)", freshness, reason, domain.IndexFreshnessFresh)
	}
}

func TestRebuildTraceIndexReplacesPriorRowsWithEmptyTrace(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}
	if err := InsertTasks(context.Background(), store, []domain.Task{{RowOrder: 1, ID: "aa/bbbbb", Status: domain.TaskStatusFailed}}); err != nil {
		t.Fatalf("InsertTasks(seed) returned error: %v", err)
	}

	runDir := domain.RunDir{Path: t.TempDir()}
	tracePath := filepath.Join(runDir.Path, "trace.txt")
	if err := os.WriteFile(tracePath, []byte("hash\tstatus\n"), 0o644); err != nil {
		t.Fatalf("write empty trace: %v", err)
	}
	traceInfo, err := os.Stat(tracePath)
	if err != nil {
		t.Fatalf("stat empty trace: %v", err)
	}
	traceSource := domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: tracePath, ModTime: traceInfo.ModTime(), Size: traceInfo.Size()}

	metadata, err := RebuildTraceIndex(context.Background(), store, runDir, domain.ArtifactSet{
		RunDir: runDir,
		Mode:   domain.IndexModeTraceBacked,
		Trace:  &traceSource,
	})
	if err != nil {
		t.Fatalf("RebuildTraceIndex(empty trace) returned error: %v", err)
	}
	if metadata.TaskCount != 0 {
		t.Fatalf("TaskCount = %d, want 0 for header-only trace", metadata.TaskCount)
	}
	if metadata.Freshness != domain.IndexFreshnessFresh || metadata.StaleReason != "" {
		t.Fatalf("Freshness/StaleReason = %q/%q, want fresh with empty reason", metadata.Freshness, metadata.StaleReason)
	}
	gotTasks, err := QueryTasks(context.Background(), store, domain.TaskQuery{})
	if err != nil {
		t.Fatalf("QueryTasks() after empty rebuild returned error: %v", err)
	}
	if len(gotTasks) != 0 {
		t.Fatalf("tasks after empty rebuild = %#v, want empty", gotTasks)
	}
}

func TestRebuildTraceIndexRejectsNonTraceBackedArtifacts(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	runDir := domain.RunDir{Path: t.TempDir()}
	logSource := domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: filepath.Join(runDir.Path, ".nextflow.log")}
	tests := []struct {
		name      string
		artifacts domain.ArtifactSet
		want      string
	}{
		{
			name:      "log-only",
			artifacts: domain.ArtifactSet{RunDir: runDir, Mode: domain.IndexModeLogOnly, Log: &logSource},
			want:      "trace-backed",
		},
		{
			name:      "unsupported",
			artifacts: domain.ArtifactSet{RunDir: runDir, Mode: domain.IndexModeUnsupported},
			want:      "trace-backed",
		},
		{
			name:      "missing trace source",
			artifacts: domain.ArtifactSet{RunDir: runDir, Mode: domain.IndexModeTraceBacked},
			want:      "trace source",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			metadata, err := RebuildTraceIndex(context.Background(), store, runDir, tt.artifacts)
			if err == nil {
				t.Fatalf("RebuildTraceIndex(%s) returned nil error and metadata %#v", tt.name, metadata)
			}
			if !strings.Contains(err.Error(), tt.want) {
				t.Fatalf("error = %q, want it to mention %q", err.Error(), tt.want)
			}
		})
	}
}

func TestEnsureFreshIndexRebuildsMissingTraceBackedIndex(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	workdir := filepath.Join(runDir.Path, "work", "AB", "C123DEF")
	traceSource := writeTraceSource(t, runDir, "trace.csv", []string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AB/C123DEF,FAILED,ALIGN,ALIGN (sample-1),sample-1," + workdir + ",137,1m,60s,2,4 GB",
		"",
	})
	artifacts := domain.ArtifactSet{
		RunDir: runDir,
		Mode:   domain.IndexModeTraceBacked,
		Trace:  &traceSource,
	}

	if _, err := os.Stat(run.IndexPath(runDir)); !os.IsNotExist(err) {
		t.Fatalf("test setup index path stat error = %v, want missing index", err)
	}

	startedAt := time.Now()
	store, metadata, err := EnsureFreshIndex(context.Background(), runDir, artifacts)
	finishedAt := time.Now()
	if err != nil {
		t.Fatalf("EnsureFreshIndex(missing trace-backed index) returned error: %v", err)
	}
	if store == nil {
		t.Fatalf("EnsureFreshIndex() store = nil, want open store")
	}
	defer store.Close()

	if store.Path != run.IndexPath(runDir) {
		t.Fatalf("Store.Path = %q, want %q", store.Path, run.IndexPath(runDir))
	}
	if metadata.SchemaVersion != SchemaVersion {
		t.Fatalf("SchemaVersion = %d, want %d", metadata.SchemaVersion, SchemaVersion)
	}
	if metadata.RunDir != runDir.Path {
		t.Fatalf("RunDir = %q, want %q", metadata.RunDir, runDir.Path)
	}
	if metadata.IndexPath != run.IndexPath(runDir) {
		t.Fatalf("IndexPath = %q, want %q", metadata.IndexPath, run.IndexPath(runDir))
	}
	if metadata.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("Mode = %q, want %q", metadata.Mode, domain.IndexModeTraceBacked)
	}
	if metadata.Freshness != domain.IndexFreshnessFresh || metadata.StaleReason != "" {
		t.Fatalf("Freshness/StaleReason = %q/%q, want fresh with empty reason", metadata.Freshness, metadata.StaleReason)
	}
	if metadata.TaskCount != 1 {
		t.Fatalf("TaskCount = %d, want 1", metadata.TaskCount)
	}
	if metadata.BuiltAt.IsZero() || metadata.BuiltAt.Before(startedAt) || metadata.BuiltAt.After(finishedAt) {
		t.Fatalf("BuiltAt = %s, want non-zero timestamp between %s and %s", metadata.BuiltAt.Format(time.RFC3339Nano), startedAt.Format(time.RFC3339Nano), finishedAt.Format(time.RFC3339Nano))
	}
	assertSourceFingerprint(t, "Trace", metadata.Trace, traceSource)

	gotTasks, err := QueryTasks(context.Background(), store, domain.TaskQuery{})
	if err != nil {
		t.Fatalf("QueryTasks() after EnsureFreshIndex() returned error: %v", err)
	}
	exitCode := 137
	wantTasks := []domain.Task{{
		RowOrder: 1,
		ID:       "ab/c123def",
		Status:   domain.TaskStatusFailed,
		Process:  "ALIGN",
		Name:     "ALIGN (sample-1)",
		Tag:      "sample-1",
		Workdir:  filepath.Clean(workdir),
		Exit:     &exitCode,
		Duration: "1m",
		Realtime: "60s",
		CPUs:     "2",
		Memory:   "4 GB",
	}}
	if !reflect.DeepEqual(gotTasks, wantTasks) {
		t.Fatalf("tasks after EnsureFreshIndex() = %#v, want %#v", gotTasks, wantTasks)
	}
}

func TestEnsureFreshIndexUsesFreshMetadataWithoutReparsingTrace(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	workdir := filepath.Join(runDir.Path, "work", "AA", "111111")
	traceSource := writeTraceSource(t, runDir, "trace.csv", []string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AA/111111,COMPLETED,ALIGN,ALIGN (sample-1),sample-1," + workdir + ",0,1m,60s,2,4 GB",
		"",
	})
	artifacts := domain.ArtifactSet{RunDir: runDir, Mode: domain.IndexModeTraceBacked, Trace: &traceSource}

	seedStore, err := OpenStore(context.Background(), runDir)
	if err != nil {
		t.Fatalf("OpenStore(seed) returned error: %v", err)
	}
	seedMetadata, err := RebuildTraceIndex(context.Background(), seedStore, runDir, artifacts)
	if err != nil {
		_ = seedStore.Close()
		t.Fatalf("RebuildTraceIndex(seed) returned error: %v", err)
	}
	if err := seedStore.Close(); err != nil {
		t.Fatalf("close seed store: %v", err)
	}

	if err := os.WriteFile(traceSource.Path, []byte("status\nFAILED\n"), 0o644); err != nil {
		t.Fatalf("overwrite selected trace with invalid content: %v", err)
	}

	store, metadata, err := EnsureFreshIndex(context.Background(), runDir, artifacts)
	if err != nil {
		t.Fatalf("EnsureFreshIndex(fresh metadata) returned error: %v", err)
	}
	if store == nil {
		t.Fatalf("EnsureFreshIndex() store = nil, want open store")
	}
	defer store.Close()

	if !reflect.DeepEqual(metadata, seedMetadata) {
		t.Fatalf("metadata = %#v, want existing fresh metadata %#v", metadata, seedMetadata)
	}
	gotTasks, err := QueryTasks(context.Background(), store, domain.TaskQuery{})
	if err != nil {
		t.Fatalf("QueryTasks() after fresh EnsureFreshIndex() returned error: %v", err)
	}
	if len(gotTasks) != 1 || gotTasks[0].ID != "aa/111111" || gotTasks[0].Status != domain.TaskStatusCompleted {
		t.Fatalf("tasks after fresh EnsureFreshIndex() = %#v, want original completed task", gotTasks)
	}
}

func TestEnsureFreshIndexRebuildsStaleTraceBackedIndex(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	firstWorkdir := filepath.Join(runDir.Path, "work", "AA", "111111")
	firstTraceSource := writeTraceSource(t, runDir, "trace.csv", []string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AA/111111,COMPLETED,ALIGN,ALIGN (sample-1),sample-1," + firstWorkdir + ",0,1m,60s,2,4 GB",
		"",
	})
	artifacts := domain.ArtifactSet{RunDir: runDir, Mode: domain.IndexModeTraceBacked, Trace: &firstTraceSource}

	seedStore, err := OpenStore(context.Background(), runDir)
	if err != nil {
		t.Fatalf("OpenStore(seed) returned error: %v", err)
	}
	oldMetadata, err := RebuildTraceIndex(context.Background(), seedStore, runDir, artifacts)
	if err != nil {
		_ = seedStore.Close()
		t.Fatalf("RebuildTraceIndex(seed) returned error: %v", err)
	}
	if err := seedStore.Close(); err != nil {
		t.Fatalf("close seed store: %v", err)
	}

	secondWorkdir := filepath.Join(runDir.Path, "work", "BB", "222222")
	updatedTraceSource := writeTraceSource(t, runDir, "trace.csv", []string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"BB/222222,FAILED,QUANTIFY,QUANTIFY (sample-2),sample-2," + secondWorkdir + ",137,2m,120s,4,8 GB",
		"",
	})
	if updatedTraceSource.Size == firstTraceSource.Size && updatedTraceSource.ModTime.Equal(firstTraceSource.ModTime) {
		t.Fatalf("test setup updated trace fingerprint = original fingerprint %#v", updatedTraceSource)
	}
	artifacts.Trace = &updatedTraceSource

	store, metadata, err := EnsureFreshIndex(context.Background(), runDir, artifacts)
	if err != nil {
		t.Fatalf("EnsureFreshIndex(stale trace-backed index) returned error: %v", err)
	}
	if store == nil {
		t.Fatalf("EnsureFreshIndex() store = nil, want open store")
	}
	defer store.Close()

	if metadata.TaskCount != 1 {
		t.Fatalf("TaskCount = %d, want 1 rebuilt task", metadata.TaskCount)
	}
	if metadata.Freshness != domain.IndexFreshnessFresh || metadata.StaleReason != "" {
		t.Fatalf("Freshness/StaleReason = %q/%q, want fresh with empty reason", metadata.Freshness, metadata.StaleReason)
	}
	assertSourceFingerprint(t, "Trace", metadata.Trace, updatedTraceSource)
	if reflect.DeepEqual(metadata.Trace, oldMetadata.Trace) {
		t.Fatalf("Trace metadata was not updated: got %#v, old %#v", metadata.Trace, oldMetadata.Trace)
	}

	gotTasks, err := QueryTasks(context.Background(), store, domain.TaskQuery{})
	if err != nil {
		t.Fatalf("QueryTasks() after stale EnsureFreshIndex() returned error: %v", err)
	}
	if len(gotTasks) != 1 || gotTasks[0].ID != "bb/222222" || gotTasks[0].Status != domain.TaskStatusFailed {
		t.Fatalf("tasks after stale EnsureFreshIndex() = %#v, want rebuilt failed task only", gotTasks)
	}
}

func TestEnsureFreshIndexRejectsNonTraceBackedArtifactsWithoutOpeningStore(t *testing.T) {
	tests := []struct {
		name    string
		mode    domain.IndexMode
		withLog bool
	}{
		{name: "log-only", mode: domain.IndexModeLogOnly, withLog: true},
		{name: "unsupported", mode: domain.IndexModeUnsupported},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			runDir := domain.RunDir{Path: t.TempDir()}
			artifacts := domain.ArtifactSet{RunDir: runDir, Mode: tt.mode}
			if tt.withLog {
				logSource := writeSourceFingerprint(t, domain.SourceKindLog, filepath.Join(runDir.Path, ".nextflow.log"), "log-only failure evidence\n")
				artifacts.Log = &logSource
			}

			store, metadata, err := EnsureFreshIndex(context.Background(), runDir, artifacts)
			if err == nil {
				t.Fatalf("EnsureFreshIndex(%s) returned nil error with store %#v and metadata %#v", tt.name, store, metadata)
			}
			if store != nil {
				_ = store.Close()
				t.Fatalf("EnsureFreshIndex(%s) store = %#v, want nil on error", tt.name, store)
			}
			if !reflect.DeepEqual(metadata, domain.IndexMetadata{}) {
				t.Fatalf("EnsureFreshIndex(%s) metadata = %#v, want zero metadata on error", tt.name, metadata)
			}
			if !strings.Contains(err.Error(), "trace-backed") {
				t.Fatalf("EnsureFreshIndex(%s) error = %q, want it to mention trace-backed", tt.name, err.Error())
			}
			if _, statErr := os.Stat(filepath.Join(runDir.Path, run.IndexDirName)); !os.IsNotExist(statErr) {
				t.Fatalf("EnsureFreshIndex(%s) cache directory stat error = %v, want cache directory not created", tt.name, statErr)
			}
		})
	}
}

func TestIndexDiagnosticsReportsFreshExistingIndex(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	workdir := filepath.Join(runDir.Path, "work", "AA", "111111")
	traceSource := writeTraceSource(t, runDir, "trace.csv", []string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AA/111111,COMPLETED,ALIGN,ALIGN (sample-1),sample-1," + workdir + ",0,1m,60s,2,4 GB",
		"",
	})
	artifacts := domain.ArtifactSet{RunDir: runDir, Mode: domain.IndexModeTraceBacked, Trace: &traceSource}

	seedStore, err := OpenStore(context.Background(), runDir)
	if err != nil {
		t.Fatalf("OpenStore(seed) returned error: %v", err)
	}
	seedMetadata, err := RebuildTraceIndex(context.Background(), seedStore, runDir, artifacts)
	if err != nil {
		_ = seedStore.Close()
		t.Fatalf("RebuildTraceIndex(seed) returned error: %v", err)
	}
	if err := seedStore.Close(); err != nil {
		t.Fatalf("close seed store: %v", err)
	}

	diagnostics, err := IndexDiagnostics(context.Background(), runDir, artifacts)
	if err != nil {
		t.Fatalf("IndexDiagnostics(fresh index) returned error: %v", err)
	}
	if diagnostics.RunDir != runDir {
		t.Fatalf("RunDir = %#v, want %#v", diagnostics.RunDir, runDir)
	}
	if diagnostics.Artifacts.Mode != domain.IndexModeTraceBacked {
		t.Fatalf("Artifacts.Mode = %q, want %q", diagnostics.Artifacts.Mode, domain.IndexModeTraceBacked)
	}
	if diagnostics.Metadata == nil {
		t.Fatalf("Metadata = nil, want existing index metadata")
	}
	if diagnostics.Metadata.IndexPath != run.IndexPath(runDir) {
		t.Fatalf("IndexPath = %q, want %q", diagnostics.Metadata.IndexPath, run.IndexPath(runDir))
	}
	if diagnostics.Metadata.Freshness != domain.IndexFreshnessFresh || diagnostics.Metadata.StaleReason != "" {
		t.Fatalf("Freshness/StaleReason = %q/%q, want fresh with empty reason", diagnostics.Metadata.Freshness, diagnostics.Metadata.StaleReason)
	}
	if diagnostics.Metadata.TaskCount != seedMetadata.TaskCount {
		t.Fatalf("TaskCount = %d, want %d", diagnostics.Metadata.TaskCount, seedMetadata.TaskCount)
	}
	if !diagnostics.Metadata.BuiltAt.Equal(seedMetadata.BuiltAt) {
		t.Fatalf("BuiltAt = %s, want existing build time %s", diagnostics.Metadata.BuiltAt.Format(time.RFC3339Nano), seedMetadata.BuiltAt.Format(time.RFC3339Nano))
	}
	if len(diagnostics.Diagnostics) != 0 {
		t.Fatalf("Diagnostics = %#v, want none for fresh index", diagnostics.Diagnostics)
	}
}

func TestIndexDiagnosticsReportsStaleExistingIndexWithoutRebuilding(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	firstWorkdir := filepath.Join(runDir.Path, "work", "AA", "111111")
	firstTraceSource := writeTraceSource(t, runDir, "trace.csv", []string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"AA/111111,COMPLETED,ALIGN,ALIGN (sample-1),sample-1," + firstWorkdir + ",0,1m,60s,2,4 GB",
		"",
	})
	artifacts := domain.ArtifactSet{RunDir: runDir, Mode: domain.IndexModeTraceBacked, Trace: &firstTraceSource}

	seedStore, err := OpenStore(context.Background(), runDir)
	if err != nil {
		t.Fatalf("OpenStore(seed) returned error: %v", err)
	}
	seedMetadata, err := RebuildTraceIndex(context.Background(), seedStore, runDir, artifacts)
	if err != nil {
		_ = seedStore.Close()
		t.Fatalf("RebuildTraceIndex(seed) returned error: %v", err)
	}
	if err := seedStore.Close(); err != nil {
		t.Fatalf("close seed store: %v", err)
	}

	secondWorkdir := filepath.Join(runDir.Path, "work", "BB", "222222")
	updatedTraceSource := writeTraceSource(t, runDir, "trace.csv", []string{
		"hash,status,process,name,tag,workdir,exit,duration,realtime,cpus,memory",
		"BB/222222,FAILED,QUANTIFY,QUANTIFY (sample-2),sample-2," + secondWorkdir + ",137,2m,120s,4,8 GB",
		"",
	})
	if updatedTraceSource.Size == firstTraceSource.Size && updatedTraceSource.ModTime.Equal(firstTraceSource.ModTime) {
		t.Fatalf("test setup updated trace fingerprint = original fingerprint %#v", updatedTraceSource)
	}
	artifacts.Trace = &updatedTraceSource

	diagnostics, err := IndexDiagnostics(context.Background(), runDir, artifacts)
	if err != nil {
		t.Fatalf("IndexDiagnostics(stale index) returned error: %v", err)
	}
	if diagnostics.Metadata == nil {
		t.Fatalf("Metadata = nil, want stale index metadata")
	}
	if diagnostics.Metadata.Freshness != domain.IndexFreshnessStale {
		t.Fatalf("Freshness = %q, want %q", diagnostics.Metadata.Freshness, domain.IndexFreshnessStale)
	}
	if diagnostics.Metadata.StaleReason != "selected trace changed" {
		t.Fatalf("StaleReason = %q, want selected trace changed", diagnostics.Metadata.StaleReason)
	}
	if diagnostics.Metadata.TaskCount != seedMetadata.TaskCount {
		t.Fatalf("TaskCount = %d, want existing task count %d without rebuild", diagnostics.Metadata.TaskCount, seedMetadata.TaskCount)
	}
	if !diagnostics.Metadata.BuiltAt.Equal(seedMetadata.BuiltAt) {
		t.Fatalf("BuiltAt = %s, want existing build time %s without rebuild", diagnostics.Metadata.BuiltAt.Format(time.RFC3339Nano), seedMetadata.BuiltAt.Format(time.RFC3339Nano))
	}
	if len(diagnostics.Diagnostics) != 1 {
		t.Fatalf("Diagnostics = %#v, want one stale warning", diagnostics.Diagnostics)
	}
	if diagnostics.Diagnostics[0].Severity != domain.DiagnosticWarning || diagnostics.Diagnostics[0].Code != "index_stale" {
		t.Fatalf("stale diagnostic = %#v, want warning index_stale", diagnostics.Diagnostics[0])
	}

	verifyStore, err := OpenStore(context.Background(), runDir)
	if err != nil {
		t.Fatalf("OpenStore(verify) returned error: %v", err)
	}
	defer verifyStore.Close()
	gotTasks, err := QueryTasks(context.Background(), verifyStore, domain.TaskQuery{})
	if err != nil {
		t.Fatalf("QueryTasks() after IndexDiagnostics() returned error: %v", err)
	}
	if len(gotTasks) != 1 || gotTasks[0].ID != "aa/111111" || gotTasks[0].Status != domain.TaskStatusCompleted {
		t.Fatalf("tasks after IndexDiagnostics() = %#v, want original completed task without rebuild", gotTasks)
	}
}

func TestIndexDiagnosticsReportsMissingIndexWithoutCreatingCacheDirectory(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	traceSource := writeTraceSource(t, runDir, "trace.csv", []string{
		"hash,status,process",
		"AA/111111,COMPLETED,ALIGN",
		"",
	})
	artifacts := domain.ArtifactSet{RunDir: runDir, Mode: domain.IndexModeTraceBacked, Trace: &traceSource}

	diagnostics, err := IndexDiagnostics(context.Background(), runDir, artifacts)
	if err != nil {
		t.Fatalf("IndexDiagnostics(missing index) returned error: %v", err)
	}
	if diagnostics.Metadata == nil {
		t.Fatalf("Metadata = nil, want missing-index metadata for diagnostics")
	}
	if diagnostics.Metadata.IndexPath != run.IndexPath(runDir) {
		t.Fatalf("IndexPath = %q, want %q", diagnostics.Metadata.IndexPath, run.IndexPath(runDir))
	}
	if diagnostics.Metadata.Freshness != domain.IndexFreshnessMissing {
		t.Fatalf("Freshness = %q, want %q", diagnostics.Metadata.Freshness, domain.IndexFreshnessMissing)
	}
	if diagnostics.Metadata.StaleReason != "index metadata missing" {
		t.Fatalf("StaleReason = %q, want index metadata missing", diagnostics.Metadata.StaleReason)
	}
	assertSourceFingerprint(t, "Trace", diagnostics.Metadata.Trace, traceSource)
	if len(diagnostics.Diagnostics) != 1 {
		t.Fatalf("Diagnostics = %#v, want one missing-index warning", diagnostics.Diagnostics)
	}
	if diagnostics.Diagnostics[0].Severity != domain.DiagnosticWarning || diagnostics.Diagnostics[0].Code != "index_missing" {
		t.Fatalf("missing diagnostic = %#v, want warning index_missing", diagnostics.Diagnostics[0])
	}
	if _, statErr := os.Stat(filepath.Join(runDir.Path, run.IndexDirName)); !os.IsNotExist(statErr) {
		t.Fatalf("cache directory stat error = %v, want cache directory not created", statErr)
	}
}

func TestIndexDiagnosticsReportsUnsupportedArtifactsWithoutCreatingCacheDirectory(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	artifacts := domain.ArtifactSet{
		RunDir: runDir,
		Mode:   domain.IndexModeUnsupported,
		Diagnostics: []domain.Diagnostic{
			{Severity: domain.DiagnosticError, Code: "unsupported_artifacts", Message: "No supported artifacts"},
		},
	}

	diagnostics, err := IndexDiagnostics(context.Background(), runDir, artifacts)
	if err != nil {
		t.Fatalf("IndexDiagnostics(unsupported artifacts) returned error: %v", err)
	}
	if diagnostics.RunDir != runDir {
		t.Fatalf("RunDir = %#v, want %#v", diagnostics.RunDir, runDir)
	}
	if diagnostics.Artifacts.Mode != domain.IndexModeUnsupported {
		t.Fatalf("Artifacts.Mode = %q, want %q", diagnostics.Artifacts.Mode, domain.IndexModeUnsupported)
	}
	if len(diagnostics.Artifacts.Diagnostics) != 1 || diagnostics.Artifacts.Diagnostics[0].Code != "unsupported_artifacts" {
		t.Fatalf("Artifacts.Diagnostics = %#v, want input unsupported diagnostic", diagnostics.Artifacts.Diagnostics)
	}
	if diagnostics.Metadata == nil {
		t.Fatalf("Metadata = nil, want unsupported metadata for diagnostics")
	}
	if diagnostics.Metadata.IndexPath != "" {
		t.Fatalf("IndexPath = %q, want empty path so unsupported diagnostics do not imply an existing index", diagnostics.Metadata.IndexPath)
	}
	if diagnostics.Metadata.Freshness != domain.IndexFreshnessUnsupported {
		t.Fatalf("Freshness = %q, want %q", diagnostics.Metadata.Freshness, domain.IndexFreshnessUnsupported)
	}
	if diagnostics.Metadata.StaleReason != "no supported artifacts" {
		t.Fatalf("StaleReason = %q, want no supported artifacts", diagnostics.Metadata.StaleReason)
	}
	if diagnostics.Metadata.TaskCount != 0 {
		t.Fatalf("TaskCount = %d, want 0 for unsupported artifacts", diagnostics.Metadata.TaskCount)
	}
	if len(diagnostics.Diagnostics) != 0 {
		t.Fatalf("Diagnostics = %#v, want no duplicate index diagnostics for unsupported artifacts", diagnostics.Diagnostics)
	}
	if _, statErr := os.Stat(filepath.Join(runDir.Path, run.IndexDirName)); !os.IsNotExist(statErr) {
		t.Fatalf("cache directory stat error = %v, want cache directory not created", statErr)
	}
}

func TestIndexDiagnosticsRejectsNilContextWithoutCreatingCacheDirectory(t *testing.T) {
	runDir := domain.RunDir{Path: t.TempDir()}
	traceSource := writeTraceSource(t, runDir, "trace.csv", []string{"hash,status", "AA/111111,COMPLETED"})

	diagnostics, err := IndexDiagnostics(nil, runDir, domain.ArtifactSet{RunDir: runDir, Mode: domain.IndexModeTraceBacked, Trace: &traceSource})
	if err == nil {
		t.Fatalf("IndexDiagnostics(nil context) returned nil error and diagnostics %#v", diagnostics)
	}
	if !strings.Contains(err.Error(), "nil context") {
		t.Fatalf("error = %q, want nil context message", err.Error())
	}
	if _, statErr := os.Stat(filepath.Join(runDir.Path, run.IndexDirName)); !os.IsNotExist(statErr) {
		t.Fatalf("cache directory stat error = %v, want cache directory not created", statErr)
	}
}

func TestInsertTasksStoresNormalizedFieldsInSourceOrder(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	exitCode := 137
	tasks := []domain.Task{
		{
			RowOrder:     2,
			ID:           "bb/222222",
			Status:       domain.TaskStatusFailed,
			Process:      "CALL_VARIANTS",
			Name:         "sample-2",
			Tag:          "tumor",
			Workdir:      "/runs/example/work/bb/222222",
			Exit:         &exitCode,
			Duration:     "1h 2m",
			Realtime:     "62m",
			CPUs:         "8",
			Memory:       "16 GB",
			ErrorSummary: "process failed: command exited with 137",
		},
		{
			RowOrder:     1,
			ID:           "aa/111111",
			Status:       domain.TaskStatusCompleted,
			Process:      "ALIGN",
			Name:         "sample-1",
			Tag:          "normal",
			Workdir:      "/runs/example/work/aa/111111",
			Exit:         nil,
			Duration:     "5m",
			Realtime:     "300s",
			CPUs:         "2",
			Memory:       "4 GB",
			ErrorSummary: "",
		},
	}

	if err := InsertTasks(context.Background(), store, tasks); err != nil {
		t.Fatalf("InsertTasks() returned error: %v", err)
	}

	rows, err := store.DB.Query(`
		SELECT
			row_order,
			id,
			status,
			process,
			name,
			tag,
			workdir,
			exit_code,
			duration,
			realtime,
			cpus,
			memory,
			error_summary
		FROM tasks
		ORDER BY row_order
	`)
	if err != nil {
		t.Fatalf("query tasks: %v", err)
	}
	defer rows.Close()

	type taskRow struct {
		rowOrder     int64
		id           string
		status       string
		process      string
		name         string
		tag          string
		workdir      string
		exit         sql.NullInt64
		duration     string
		realtime     string
		cpus         string
		memory       string
		errorSummary string
	}
	var got []taskRow
	for rows.Next() {
		var row taskRow
		if err := rows.Scan(
			&row.rowOrder,
			&row.id,
			&row.status,
			&row.process,
			&row.name,
			&row.tag,
			&row.workdir,
			&row.exit,
			&row.duration,
			&row.realtime,
			&row.cpus,
			&row.memory,
			&row.errorSummary,
		); err != nil {
			t.Fatalf("scan task row: %v", err)
		}
		got = append(got, row)
	}
	if err := rows.Err(); err != nil {
		t.Fatalf("iterate task rows: %v", err)
	}
	if len(got) != 2 {
		t.Fatalf("task row count = %d, want 2", len(got))
	}

	first := got[0]
	if first.rowOrder != 1 || first.id != "aa/111111" || first.status != string(domain.TaskStatusCompleted) || first.process != "ALIGN" || first.name != "sample-1" || first.tag != "normal" || first.workdir != "/runs/example/work/aa/111111" || first.duration != "5m" || first.realtime != "300s" || first.cpus != "2" || first.memory != "4 GB" || first.errorSummary != "" {
		t.Fatalf("first task row = %#v, want normalized completed row in source order", first)
	}
	if first.exit.Valid {
		t.Fatalf("first task exit_code.Valid = true, want false for nil task exit")
	}

	second := got[1]
	if second.rowOrder != 2 || second.id != "bb/222222" || second.status != string(domain.TaskStatusFailed) || second.process != "CALL_VARIANTS" || second.name != "sample-2" || second.tag != "tumor" || second.workdir != "/runs/example/work/bb/222222" || second.duration != "1h 2m" || second.realtime != "62m" || second.cpus != "8" || second.memory != "16 GB" || second.errorSummary != "process failed: command exited with 137" {
		t.Fatalf("second task row = %#v, want normalized failed row in source order", second)
	}
	if !second.exit.Valid || second.exit.Int64 != int64(exitCode) {
		t.Fatalf("second task exit_code = (%d, valid=%t), want (%d, valid=true)", second.exit.Int64, second.exit.Valid, exitCode)
	}
}

func TestInsertTasksReplacesPriorRowsForRebuild(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	initialExit := 1
	if err := InsertTasks(context.Background(), store, []domain.Task{
		{
			RowOrder: 1,
			ID:       "old/task",
			Status:   domain.TaskStatusFailed,
			Process:  "OLD_PROCESS",
			Exit:     &initialExit,
		},
	}); err != nil {
		t.Fatalf("initial InsertTasks() returned error: %v", err)
	}

	if err := InsertTasks(context.Background(), store, []domain.Task{
		{
			RowOrder: 10,
			ID:       "new/task",
			Status:   domain.TaskStatusCached,
			Process:  "NEW_PROCESS",
		},
	}); err != nil {
		t.Fatalf("replacement InsertTasks() returned error: %v", err)
	}

	var count int
	if err := store.DB.QueryRow("SELECT COUNT(*) FROM tasks").Scan(&count); err != nil {
		t.Fatalf("count replacement tasks: %v", err)
	}
	if count != 1 {
		t.Fatalf("task row count after replacement = %d, want 1", count)
	}

	var oldCount int
	if err := store.DB.QueryRow("SELECT COUNT(*) FROM tasks WHERE id = ?", "old/task").Scan(&oldCount); err != nil {
		t.Fatalf("count old task rows: %v", err)
	}
	if oldCount != 0 {
		t.Fatalf("old task row count after replacement = %d, want 0", oldCount)
	}

	var rowOrder int64
	var process string
	if err := store.DB.QueryRow("SELECT row_order, process FROM tasks WHERE id = ?", "new/task").Scan(&rowOrder, &process); err != nil {
		t.Fatalf("query replacement task: %v", err)
	}
	if rowOrder != 10 || process != "NEW_PROCESS" {
		t.Fatalf("replacement row = (row_order=%d, process=%q), want (10, %q)", rowOrder, process, "NEW_PROCESS")
	}

	if err := InsertTasks(context.Background(), store, nil); err != nil {
		t.Fatalf("empty InsertTasks() returned error: %v", err)
	}
	if err := store.DB.QueryRow("SELECT COUNT(*) FROM tasks").Scan(&count); err != nil {
		t.Fatalf("count tasks after empty replacement: %v", err)
	}
	if count != 0 {
		t.Fatalf("task row count after empty replacement = %d, want 0", count)
	}
}

func TestInsertTasksRollsBackReplacementWhenInsertFails(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	_, err := store.DB.Exec(`
		INSERT INTO tasks (row_order, id, status, process)
		VALUES (?, ?, ?, ?)
	`, int64(1), "original/task", string(domain.TaskStatusCompleted), "ORIGINAL_PROCESS")
	if err != nil {
		t.Fatalf("insert original task row: %v", err)
	}
	_, err = store.DB.Exec(`
		CREATE TRIGGER fail_task_insert
		BEFORE INSERT ON tasks
		BEGIN
			SELECT RAISE(ABORT, 'task insert blocked');
		END
	`)
	if err != nil {
		t.Fatalf("create failing task insert trigger: %v", err)
	}

	if err := InsertTasks(context.Background(), store, []domain.Task{
		{
			RowOrder: 2,
			ID:       "replacement/task",
			Status:   domain.TaskStatusRunning,
			Process:  "REPLACEMENT_PROCESS",
		},
	}); err == nil {
		t.Fatalf("InsertTasks() returned nil error with failing task insert trigger")
	}

	var count int
	if err := store.DB.QueryRow("SELECT COUNT(*) FROM tasks").Scan(&count); err != nil {
		t.Fatalf("count tasks after failed replacement: %v", err)
	}
	if count != 1 {
		t.Fatalf("task row count after failed replacement = %d, want 1", count)
	}

	var process string
	if err := store.DB.QueryRow("SELECT process FROM tasks WHERE id = ?", "original/task").Scan(&process); err != nil {
		t.Fatalf("query original task after failed replacement: %v", err)
	}
	if process != "ORIGINAL_PROCESS" {
		t.Fatalf("original task process after failed replacement = %q, want %q", process, "ORIGINAL_PROCESS")
	}
}

func TestInsertTasksRejectsNilStoreOrDatabase(t *testing.T) {
	if err := InsertTasks(context.Background(), nil, nil); err == nil {
		t.Fatalf("InsertTasks(nil) returned nil error")
	}
	if err := InsertTasks(context.Background(), &Store{}, nil); err == nil {
		t.Fatalf("InsertTasks(store with nil DB) returned nil error")
	}
}

func TestQueryTasksReturnsRowsInSourceOrderAndPreservesNullableExit(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	exitCode := 137
	inserted := []domain.Task{
		{
			RowOrder:     2,
			ID:           "bb/222222",
			Status:       domain.TaskStatusFailed,
			Process:      "CALL_VARIANTS",
			Name:         "sample-2",
			Tag:          "tumor",
			Workdir:      "/runs/example/work/bb/222222",
			Exit:         &exitCode,
			Duration:     "1h 2m",
			Realtime:     "62m",
			CPUs:         "8",
			Memory:       "16 GB",
			ErrorSummary: "process failed: command exited with 137",
		},
		{
			RowOrder: 1,
			ID:       "aa/111111",
			Status:   domain.TaskStatusCompleted,
			Process:  "ALIGN",
			Name:     "sample-1",
			Tag:      "normal",
			Workdir:  "/runs/example/work/aa/111111",
			Exit:     nil,
			Duration: "5m",
			Realtime: "300s",
			CPUs:     "2",
			Memory:   "4 GB",
		},
	}
	if err := InsertTasks(context.Background(), store, inserted); err != nil {
		t.Fatalf("InsertTasks() returned error: %v", err)
	}

	got, err := QueryTasks(context.Background(), store, domain.TaskQuery{})
	if err != nil {
		t.Fatalf("QueryTasks() returned error: %v", err)
	}

	want := []domain.Task{inserted[1], inserted[0]}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("QueryTasks() = %#v, want %#v", got, want)
	}
	if got[0].Exit != nil {
		t.Fatalf("first queried task Exit = %#v, want nil", got[0].Exit)
	}
	if got[1].Exit == nil || *got[1].Exit != exitCode {
		t.Fatalf("second queried task Exit = %#v, want %d", got[1].Exit, exitCode)
	}
}

func TestQueryTasksAppliesV1FiltersCaseInsensitively(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	inserted := []domain.Task{
		{RowOrder: 3, ID: "tumor-align", Status: domain.TaskStatusCompleted, Process: "ALIGN", Name: "WES-03", Tag: "tumor"},
		{RowOrder: 1, ID: "normal-call", Status: domain.TaskStatusFailed, Process: "CALL_VARIANTS", Name: "WES-02", Tag: "normal"},
		{RowOrder: 2, ID: "tumor-call", Status: domain.TaskStatusFailed, Process: "CALL_VARIANTS", Name: "WES-01", Tag: "tumor"},
	}
	if err := InsertTasks(context.Background(), store, inserted); err != nil {
		t.Fatalf("InsertTasks() returned error: %v", err)
	}

	got, err := QueryTasks(context.Background(), store, domain.TaskQuery{
		ProcessSubstring: " call ",
		NameSubstring:    " wes-01 ",
		SampleSubstring:  " TUMOR ",
		StatusRaw:        " failed ",
	})
	if err != nil {
		t.Fatalf("QueryTasks() returned error: %v", err)
	}

	want := []domain.Task{inserted[2]}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("QueryTasks() = %#v, want %#v", got, want)
	}
}

func TestQueryTasksReturnsEmptySliceForEmptyTaskTable(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	got, err := QueryTasks(context.Background(), store, domain.TaskQuery{})
	if err != nil {
		t.Fatalf("QueryTasks() returned error: %v", err)
	}
	if got == nil {
		t.Fatalf("QueryTasks() = nil, want empty slice")
	}
	if len(got) != 0 {
		t.Fatalf("QueryTasks() length = %d, want 0: %#v", len(got), got)
	}
}

func TestQueryTasksRejectsInvalidQueryOrMissingStore(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	if _, err := QueryTasks(context.Background(), store, domain.TaskQuery{StatusRaw: "finished"}); err == nil {
		t.Fatalf("QueryTasks() returned nil error for unknown status")
	}
	if _, err := QueryTasks(context.Background(), nil, domain.TaskQuery{}); err == nil {
		t.Fatalf("QueryTasks(nil) returned nil error")
	}
	if _, err := QueryTasks(context.Background(), &Store{}, domain.TaskQuery{}); err == nil {
		t.Fatalf("QueryTasks(store with nil DB) returned nil error")
	}
}

func TestCountTasksByStatusAggregatesRowsInLexicalStatusOrder(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	tasks := []domain.Task{
		{RowOrder: 1, ID: "failed-1", Status: domain.TaskStatusFailed},
		{RowOrder: 2, ID: "completed-1", Status: domain.TaskStatusCompleted},
		{RowOrder: 3, ID: "failed-2", Status: domain.TaskStatusFailed},
		{RowOrder: 4, ID: "cached-1", Status: domain.TaskStatusCached},
		{RowOrder: 5, ID: "completed-2", Status: domain.TaskStatusCompleted},
	}
	if err := InsertTasks(context.Background(), store, tasks); err != nil {
		t.Fatalf("InsertTasks() returned error: %v", err)
	}

	got, err := CountTasksByStatus(context.Background(), store)
	if err != nil {
		t.Fatalf("CountTasksByStatus() returned error: %v", err)
	}

	want := []domain.StatusCount{
		{Status: domain.TaskStatusCached, Count: 1},
		{Status: domain.TaskStatusCompleted, Count: 2},
		{Status: domain.TaskStatusFailed, Count: 2},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("CountTasksByStatus() = %#v, want %#v", got, want)
	}
}

func TestCountTasksByStatusReturnsEmptySliceForEmptyTaskTable(t *testing.T) {
	store, cleanup := openTempSQLiteStore(t)
	defer cleanup()

	if err := InitializeSchema(context.Background(), store); err != nil {
		t.Fatalf("InitializeSchema() returned error: %v", err)
	}

	got, err := CountTasksByStatus(context.Background(), store)
	if err != nil {
		t.Fatalf("CountTasksByStatus() returned error: %v", err)
	}
	if got == nil {
		t.Fatalf("CountTasksByStatus() = nil, want empty slice")
	}
	if len(got) != 0 {
		t.Fatalf("CountTasksByStatus() length = %d, want 0: %#v", len(got), got)
	}
}

func TestCountTasksByStatusRejectsNilStoreOrDatabase(t *testing.T) {
	if _, err := CountTasksByStatus(context.Background(), nil); err == nil {
		t.Fatalf("CountTasksByStatus(nil) returned nil error")
	}
	if _, err := CountTasksByStatus(context.Background(), &Store{}); err == nil {
		t.Fatalf("CountTasksByStatus(store with nil DB) returned nil error")
	}
}

func writeTraceSource(t *testing.T, runDir domain.RunDir, name string, lines []string) domain.SourceFingerprint {
	t.Helper()

	return writeSourceFingerprint(t, domain.SourceKindTrace, filepath.Join(runDir.Path, name), strings.Join(lines, "\n"))
}

func writeSourceFingerprint(t *testing.T, kind domain.SourceKind, path string, content string) domain.SourceFingerprint {
	t.Helper()

	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		t.Fatalf("create parent directory for %s source: %v", kind, err)
	}
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s source %q: %v", kind, path, err)
	}
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("stat %s source %q: %v", kind, path, err)
	}

	return domain.SourceFingerprint{Kind: kind, Path: path, ModTime: info.ModTime(), Size: info.Size()}
}

func openTempSQLiteStore(t *testing.T) (*Store, func()) {
	t.Helper()

	path := filepath.Join(t.TempDir(), "index.sqlite")
	db, err := sql.Open("sqlite", path)
	if err != nil {
		t.Fatalf("open temp sqlite database: %v", err)
	}
	cleanup := func() {
		if err := db.Close(); err != nil {
			t.Fatalf("close temp sqlite database: %v", err)
		}
	}

	return &Store{Path: path, DB: db}, cleanup
}

func requireTableColumns(t *testing.T, db *sql.DB, table string, want []string) {
	t.Helper()

	rows, err := db.Query(fmt.Sprintf("PRAGMA table_info(%s)", table))
	if err != nil {
		t.Fatalf("table_info(%s): %v", table, err)
	}
	defer rows.Close()

	columns := map[string]bool{}
	for rows.Next() {
		var cid int
		var name string
		var dataType string
		var notNull int
		var defaultValue sql.NullString
		var primaryKey int
		if err := rows.Scan(&cid, &name, &dataType, &notNull, &defaultValue, &primaryKey); err != nil {
			t.Fatalf("scan table_info(%s): %v", table, err)
		}
		columns[name] = true
	}
	if err := rows.Err(); err != nil {
		t.Fatalf("iterate table_info(%s): %v", table, err)
	}
	if len(columns) == 0 {
		t.Fatalf("table %q does not exist or has no columns", table)
	}

	for _, column := range want {
		if !columns[column] {
			t.Fatalf("table %q columns = %s, missing %q", table, sortedColumnList(columns), column)
		}
	}
}

func assertSourceFingerprint(t *testing.T, label string, got *domain.SourceFingerprint, want domain.SourceFingerprint) {
	t.Helper()

	if got == nil {
		t.Fatalf("%s fingerprint = nil, want %#v", label, want)
	}
	if got.Kind != want.Kind {
		t.Fatalf("%s Kind = %q, want %q", label, got.Kind, want.Kind)
	}
	if got.Path != want.Path {
		t.Fatalf("%s Path = %q, want %q", label, got.Path, want.Path)
	}
	if !got.ModTime.Equal(want.ModTime) {
		t.Fatalf("%s ModTime = %s, want %s", label, got.ModTime.Format(time.RFC3339Nano), want.ModTime.Format(time.RFC3339Nano))
	}
	if got.Size != want.Size {
		t.Fatalf("%s Size = %d, want %d", label, got.Size, want.Size)
	}
}

func sortedColumnList(columns map[string]bool) string {
	names := make([]string, 0, len(columns))
	for name := range columns {
		names = append(names, name)
	}
	sort.Strings(names)
	return strings.Join(names, ", ")
}
