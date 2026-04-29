package index

import (
	"context"
	"database/sql"
	"fmt"
	"os"
	"path/filepath"
	"time"

	"github.com/mskilab-org/gosh/internal/domain"
	"github.com/mskilab-org/gosh/internal/run"
	"github.com/mskilab-org/gosh/internal/tasks"
	"github.com/mskilab-org/gosh/internal/trace"
	_ "modernc.org/sqlite"
)

const SchemaVersion = 1

type Store struct {
	Path string
	DB   *sql.DB
}

func OpenStore(ctx context.Context, runDir domain.RunDir) (*Store, error) {
	if ctx == nil {
		return nil, fmt.Errorf("open store: nil context")
	}
	if err := ctx.Err(); err != nil {
		return nil, fmt.Errorf("open store: %w", err)
	}
	if runDir.Path == "" {
		return nil, fmt.Errorf("open store: empty run dir")
	}

	runInfo, err := os.Stat(runDir.Path)
	if err != nil {
		return nil, fmt.Errorf("open store run dir %q: %w", runDir.Path, err)
	}
	if !runInfo.IsDir() {
		return nil, fmt.Errorf("open store run dir %q: not a directory", runDir.Path)
	}

	indexPath := run.IndexPath(runDir)
	indexDir := filepath.Dir(indexPath)
	if err := os.MkdirAll(indexDir, 0o755); err != nil {
		return nil, fmt.Errorf("create index directory %q: %w", indexDir, err)
	}
	if err := ctx.Err(); err != nil {
		return nil, fmt.Errorf("open store: %w", err)
	}

	db, err := sql.Open("sqlite", indexPath)
	if err != nil {
		return nil, fmt.Errorf("open sqlite index %q: %w", indexPath, err)
	}
	store := &Store{Path: indexPath, DB: db}
	if err := db.PingContext(ctx); err != nil {
		_ = db.Close()
		return nil, fmt.Errorf("ping sqlite index %q: %w", indexPath, err)
	}

	return store, nil
}

func (s *Store) Close() error {
	if s == nil || s.DB == nil {
		return nil
	}
	return s.DB.Close()
}

func validateStore(ctx context.Context, store *Store, operation string) error {
	if err := ctx.Err(); err != nil {
		return err
	}
	if store == nil {
		return fmt.Errorf("%s: nil store", operation)
	}
	if store.DB == nil {
		return fmt.Errorf("%s: nil database", operation)
	}
	return nil
}

func copySourceFingerprint(kind domain.SourceKind, fingerprint *domain.SourceFingerprint) *domain.SourceFingerprint {
	if fingerprint == nil {
		return nil
	}
	copied := *fingerprint
	copied.Kind = kind
	return &copied
}

func sourceFingerprintsMatch(persisted *domain.SourceFingerprint, selected *domain.SourceFingerprint) bool {
	if persisted == nil || selected == nil {
		return persisted == nil && selected == nil
	}
	return persisted.Path == selected.Path && persisted.ModTime.Equal(selected.ModTime) && persisted.Size == selected.Size
}

func effectiveIndexPath(store *Store, runDirPath string) string {
	if store != nil && store.Path != "" {
		return store.Path
	}
	if runDirPath == "" {
		return ""
	}
	return run.IndexPath(domain.RunDir{Path: runDirPath})
}

func metadataForArtifacts(runDirPath string, indexPath string, artifacts domain.ArtifactSet) domain.IndexMetadata {
	return domain.IndexMetadata{
		SchemaVersion: SchemaVersion,
		RunDir:        runDirPath,
		IndexPath:     indexPath,
		Mode:          artifacts.Mode,
		Trace:         copySourceFingerprint(domain.SourceKindTrace, artifacts.Trace),
		Log:           copySourceFingerprint(domain.SourceKindLog, artifacts.Log),
	}
}

func InitializeSchema(ctx context.Context, store *Store) error {
	if err := validateStore(ctx, store, "initialize schema"); err != nil {
		return err
	}

	tx, err := store.DB.BeginTx(ctx, nil)
	if err != nil {
		return fmt.Errorf("begin schema transaction: %w", err)
	}
	committed := false
	defer func() {
		if !committed {
			_ = tx.Rollback()
		}
	}()

	statements := []string{
		fmt.Sprintf("PRAGMA user_version = %d", SchemaVersion),
		`CREATE TABLE IF NOT EXISTS index_metadata (
			id INTEGER NOT NULL PRIMARY KEY CHECK (id = 1),
			schema_version INTEGER NOT NULL,
			run_dir TEXT NOT NULL,
			index_path TEXT NOT NULL,
			mode TEXT NOT NULL,
			built_at TEXT NOT NULL,
			freshness TEXT NOT NULL,
			stale_reason TEXT NOT NULL DEFAULT '',
			task_count INTEGER NOT NULL DEFAULT 0
		)`,
		`CREATE TABLE IF NOT EXISTS source_fingerprints (
			kind TEXT NOT NULL PRIMARY KEY,
			path TEXT NOT NULL,
			mod_time TEXT NOT NULL,
			size INTEGER NOT NULL
		)`,
		`CREATE TABLE IF NOT EXISTS tasks (
			row_order INTEGER NOT NULL,
			id TEXT NOT NULL PRIMARY KEY,
			status TEXT NOT NULL,
			process TEXT NOT NULL DEFAULT '',
			name TEXT NOT NULL DEFAULT '',
			tag TEXT NOT NULL DEFAULT '',
			workdir TEXT NOT NULL DEFAULT '',
			exit_code INTEGER,
			duration TEXT NOT NULL DEFAULT '',
			realtime TEXT NOT NULL DEFAULT '',
			cpus TEXT NOT NULL DEFAULT '',
			memory TEXT NOT NULL DEFAULT '',
			error_summary TEXT NOT NULL DEFAULT ''
		)`,
		`CREATE INDEX IF NOT EXISTS tasks_row_order_idx ON tasks (row_order)`,
		`CREATE INDEX IF NOT EXISTS tasks_status_row_order_idx ON tasks (status, row_order)`,
		`CREATE INDEX IF NOT EXISTS tasks_workdir_idx ON tasks (workdir)`,
	}

	for _, statement := range statements {
		if _, err := tx.ExecContext(ctx, statement); err != nil {
			return fmt.Errorf("apply schema statement: %w", err)
		}
	}

	if err := tx.Commit(); err != nil {
		return fmt.Errorf("commit schema transaction: %w", err)
	}
	committed = true
	return nil
}

func ReadMetadata(ctx context.Context, store *Store) (*domain.IndexMetadata, error) {
	if err := validateStore(ctx, store, "read metadata"); err != nil {
		return nil, err
	}

	var metadata domain.IndexMetadata
	var mode string
	var builtAt string
	var freshness string
	err := store.DB.QueryRowContext(ctx, `
		SELECT
			schema_version,
			run_dir,
			index_path,
			mode,
			built_at,
			freshness,
			stale_reason,
			task_count
		FROM index_metadata
		WHERE id = 1
	`).Scan(
		&metadata.SchemaVersion,
		&metadata.RunDir,
		&metadata.IndexPath,
		&mode,
		&builtAt,
		&freshness,
		&metadata.StaleReason,
		&metadata.TaskCount,
	)
	if err != nil {
		if err == sql.ErrNoRows {
			return nil, nil
		}
		return nil, fmt.Errorf("read metadata row: %w", err)
	}

	parsedBuiltAt, err := time.Parse(time.RFC3339Nano, builtAt)
	if err != nil {
		return nil, fmt.Errorf("parse metadata built_at %q: %w", builtAt, err)
	}
	metadata.Mode = domain.IndexMode(mode)
	metadata.BuiltAt = parsedBuiltAt
	metadata.Freshness = domain.IndexFreshness(freshness)

	rows, err := store.DB.QueryContext(ctx, `
		SELECT kind, path, mod_time, size
		FROM source_fingerprints
		WHERE kind IN (?, ?)
	`, string(domain.SourceKindTrace), string(domain.SourceKindLog))
	if err != nil {
		return nil, fmt.Errorf("read source fingerprints: %w", err)
	}
	defer rows.Close()

	for rows.Next() {
		var kind string
		var path string
		var modTime string
		var size int64
		if err := rows.Scan(&kind, &path, &modTime, &size); err != nil {
			return nil, fmt.Errorf("scan source fingerprint: %w", err)
		}
		parsedModTime, err := time.Parse(time.RFC3339Nano, modTime)
		if err != nil {
			return nil, fmt.Errorf("parse %s source mod_time %q: %w", kind, modTime, err)
		}

		fingerprint := &domain.SourceFingerprint{
			Kind:    domain.SourceKind(kind),
			Path:    path,
			ModTime: parsedModTime,
			Size:    size,
		}
		switch fingerprint.Kind {
		case domain.SourceKindTrace:
			metadata.Trace = fingerprint
		case domain.SourceKindLog:
			metadata.Log = fingerprint
		}
	}
	if err := rows.Err(); err != nil {
		return nil, fmt.Errorf("iterate source fingerprints: %w", err)
	}

	return &metadata, nil
}

func WriteMetadata(ctx context.Context, store *Store, metadata domain.IndexMetadata) error {
	if err := validateStore(ctx, store, "write metadata"); err != nil {
		return err
	}

	tx, err := store.DB.BeginTx(ctx, nil)
	if err != nil {
		return fmt.Errorf("begin metadata transaction: %w", err)
	}
	committed := false
	defer func() {
		if !committed {
			_ = tx.Rollback()
		}
	}()

	_, err = tx.ExecContext(ctx, `
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
		ON CONFLICT(id) DO UPDATE SET
			schema_version = excluded.schema_version,
			run_dir = excluded.run_dir,
			index_path = excluded.index_path,
			mode = excluded.mode,
			built_at = excluded.built_at,
			freshness = excluded.freshness,
			stale_reason = excluded.stale_reason,
			task_count = excluded.task_count
	`,
		metadata.SchemaVersion,
		metadata.RunDir,
		metadata.IndexPath,
		string(metadata.Mode),
		metadata.BuiltAt.Format(time.RFC3339Nano),
		string(metadata.Freshness),
		metadata.StaleReason,
		metadata.TaskCount,
	)
	if err != nil {
		return fmt.Errorf("write metadata row: %w", err)
	}

	_, err = tx.ExecContext(ctx, `
		DELETE FROM source_fingerprints
		WHERE kind IN (?, ?)
	`, string(domain.SourceKindTrace), string(domain.SourceKindLog))
	if err != nil {
		return fmt.Errorf("replace source fingerprints: %w", err)
	}

	writeFingerprint := func(kind domain.SourceKind, fingerprint *domain.SourceFingerprint) error {
		if fingerprint == nil {
			return nil
		}
		_, err := tx.ExecContext(ctx, `
			INSERT INTO source_fingerprints (kind, path, mod_time, size)
			VALUES (?, ?, ?, ?)
		`, string(kind), fingerprint.Path, fingerprint.ModTime.Format(time.RFC3339Nano), fingerprint.Size)
		if err != nil {
			return fmt.Errorf("write %s source fingerprint: %w", kind, err)
		}
		return nil
	}

	if err := writeFingerprint(domain.SourceKindTrace, metadata.Trace); err != nil {
		return err
	}
	if err := writeFingerprint(domain.SourceKindLog, metadata.Log); err != nil {
		return err
	}

	if err := tx.Commit(); err != nil {
		return fmt.Errorf("commit metadata transaction: %w", err)
	}
	committed = true
	return nil
}

func CheckFreshness(ctx context.Context, store *Store, artifacts domain.ArtifactSet) (domain.IndexFreshness, string, error) {
	if ctx == nil {
		return domain.IndexFreshnessUnknown, "", fmt.Errorf("check freshness: nil context")
	}
	if err := ctx.Err(); err != nil {
		return domain.IndexFreshnessUnknown, "", err
	}
	if artifacts.Mode == domain.IndexModeUnsupported {
		return domain.IndexFreshnessUnsupported, "no supported artifacts", nil
	}

	metadata, err := ReadMetadata(ctx, store)
	if err != nil {
		return domain.IndexFreshnessUnknown, "", err
	}
	if metadata == nil {
		return domain.IndexFreshnessMissing, "index metadata missing", nil
	}
	if metadata.Mode != artifacts.Mode {
		return domain.IndexFreshnessStale, "index mode changed", nil
	}

	if !sourceFingerprintsMatch(metadata.Trace, artifacts.Trace) {
		return domain.IndexFreshnessStale, "selected trace changed", nil
	}
	if !sourceFingerprintsMatch(metadata.Log, artifacts.Log) {
		return domain.IndexFreshnessStale, "selected log changed", nil
	}

	return domain.IndexFreshnessFresh, "", nil
}

func RefreshMetadata(ctx context.Context, store *Store, runDir domain.RunDir, artifacts domain.ArtifactSet) (domain.IndexMetadata, error) {
	if ctx == nil {
		return domain.IndexMetadata{}, fmt.Errorf("refresh metadata: nil context")
	}
	if err := ctx.Err(); err != nil {
		return domain.IndexMetadata{}, err
	}

	if err := InitializeSchema(ctx, store); err != nil {
		return domain.IndexMetadata{}, fmt.Errorf("refresh metadata initialize schema: %w", err)
	}

	metadataRunDir := runDir.Path
	if metadataRunDir == "" {
		metadataRunDir = artifacts.RunDir.Path
	}

	freshness := domain.IndexFreshnessFresh
	staleReason := ""
	if artifacts.Mode == domain.IndexModeUnsupported {
		freshness = domain.IndexFreshnessUnsupported
		staleReason = "no supported artifacts"
	}

	metadata := metadataForArtifacts(metadataRunDir, effectiveIndexPath(store, metadataRunDir), artifacts)
	metadata.BuiltAt = time.Now().UTC()
	metadata.Freshness = freshness
	metadata.StaleReason = staleReason

	if err := WriteMetadata(ctx, store, metadata); err != nil {
		return domain.IndexMetadata{}, fmt.Errorf("refresh metadata write metadata: %w", err)
	}

	return metadata, nil
}

func RebuildTraceIndex(ctx context.Context, store *Store, runDir domain.RunDir, artifacts domain.ArtifactSet) (domain.IndexMetadata, error) {
	if ctx == nil {
		return domain.IndexMetadata{}, fmt.Errorf("rebuild trace index: nil context")
	}
	if err := ctx.Err(); err != nil {
		return domain.IndexMetadata{}, fmt.Errorf("rebuild trace index: %w", err)
	}
	if store == nil {
		return domain.IndexMetadata{}, fmt.Errorf("rebuild trace index: nil store")
	}
	if store.DB == nil {
		return domain.IndexMetadata{}, fmt.Errorf("rebuild trace index: nil database")
	}
	if artifacts.Mode != domain.IndexModeTraceBacked {
		return domain.IndexMetadata{}, fmt.Errorf("rebuild trace index: requires trace-backed artifacts, got %q", artifacts.Mode)
	}
	if artifacts.Trace == nil {
		return domain.IndexMetadata{}, fmt.Errorf("rebuild trace index: missing trace source")
	}

	effectiveRunDir := runDir
	if effectiveRunDir.Path == "" {
		effectiveRunDir = artifacts.RunDir
	}

	if err := InitializeSchema(ctx, store); err != nil {
		return domain.IndexMetadata{}, fmt.Errorf("rebuild trace index initialize schema: %w", err)
	}

	parsedTasks, err := trace.ParseTrace(ctx, effectiveRunDir, *artifacts.Trace)
	if err != nil {
		return domain.IndexMetadata{}, fmt.Errorf("rebuild trace index parse trace: %w", err)
	}

	if err := InsertTasks(ctx, store, parsedTasks); err != nil {
		return domain.IndexMetadata{}, fmt.Errorf("rebuild trace index insert tasks: %w", err)
	}

	metadataRunDir := effectiveRunDir.Path
	metadata := metadataForArtifacts(metadataRunDir, effectiveIndexPath(store, metadataRunDir), artifacts)
	metadata.BuiltAt = time.Now().UTC()
	metadata.Freshness = domain.IndexFreshnessFresh
	metadata.TaskCount = len(parsedTasks)

	if err := WriteMetadata(ctx, store, metadata); err != nil {
		return domain.IndexMetadata{}, fmt.Errorf("rebuild trace index write metadata: %w", err)
	}

	return metadata, nil
}

func EnsureFreshIndex(ctx context.Context, runDir domain.RunDir, artifacts domain.ArtifactSet) (*Store, domain.IndexMetadata, error) {
	if ctx == nil {
		return nil, domain.IndexMetadata{}, fmt.Errorf("ensure fresh index: nil context")
	}
	if err := ctx.Err(); err != nil {
		return nil, domain.IndexMetadata{}, fmt.Errorf("ensure fresh index: %w", err)
	}
	if artifacts.Mode != domain.IndexModeTraceBacked {
		return nil, domain.IndexMetadata{}, fmt.Errorf("ensure fresh index: requires trace-backed artifacts, got %q", artifacts.Mode)
	}
	if artifacts.Trace == nil {
		return nil, domain.IndexMetadata{}, fmt.Errorf("ensure fresh index: missing trace source")
	}

	effectiveRunDir := runDir
	if effectiveRunDir.Path == "" {
		effectiveRunDir = artifacts.RunDir
	}

	store, err := OpenStore(ctx, effectiveRunDir)
	if err != nil {
		return nil, domain.IndexMetadata{}, fmt.Errorf("ensure fresh index open store: %w", err)
	}
	closeOnError := true
	defer func() {
		if closeOnError {
			_ = store.Close()
		}
	}()

	if err := InitializeSchema(ctx, store); err != nil {
		return nil, domain.IndexMetadata{}, fmt.Errorf("ensure fresh index initialize schema: %w", err)
	}

	freshness, reason, err := CheckFreshness(ctx, store, artifacts)
	if err != nil {
		return nil, domain.IndexMetadata{}, fmt.Errorf("ensure fresh index check freshness: %w", err)
	}

	switch freshness {
	case domain.IndexFreshnessFresh:
		metadata, err := ReadMetadata(ctx, store)
		if err != nil {
			return nil, domain.IndexMetadata{}, fmt.Errorf("ensure fresh index read metadata: %w", err)
		}
		if metadata == nil {
			return nil, domain.IndexMetadata{}, fmt.Errorf("ensure fresh index: fresh index metadata missing")
		}
		closeOnError = false
		return store, *metadata, nil
	case domain.IndexFreshnessMissing, domain.IndexFreshnessStale:
		metadata, err := RebuildTraceIndex(ctx, store, effectiveRunDir, artifacts)
		if err != nil {
			return nil, domain.IndexMetadata{}, fmt.Errorf("ensure fresh index rebuild trace index: %w", err)
		}
		closeOnError = false
		return store, metadata, nil
	default:
		if reason == "" {
			reason = "cannot build a fresh trace-backed task index"
		}
		return nil, domain.IndexMetadata{}, fmt.Errorf("ensure fresh index: freshness %q is not rebuildable: %s", freshness, reason)
	}
}

func InsertTasks(ctx context.Context, store *Store, tasks []domain.Task) error {
	if err := validateStore(ctx, store, "insert tasks"); err != nil {
		return err
	}

	tx, err := store.DB.BeginTx(ctx, nil)
	if err != nil {
		return fmt.Errorf("begin tasks transaction: %w", err)
	}
	committed := false
	defer func() {
		if !committed {
			_ = tx.Rollback()
		}
	}()

	if _, err := tx.ExecContext(ctx, `DELETE FROM tasks`); err != nil {
		return fmt.Errorf("replace task rows: %w", err)
	}

	stmt, err := tx.PrepareContext(ctx, `
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
	`)
	if err != nil {
		return fmt.Errorf("prepare task insert: %w", err)
	}
	defer stmt.Close()

	for index, task := range tasks {
		var exitCode any
		if task.Exit != nil {
			exitCode = *task.Exit
		}

		_, err := stmt.ExecContext(
			ctx,
			task.RowOrder,
			task.ID,
			string(task.Status),
			task.Process,
			task.Name,
			task.Tag,
			task.Workdir,
			exitCode,
			task.Duration,
			task.Realtime,
			task.CPUs,
			task.Memory,
			task.ErrorSummary,
		)
		if err != nil {
			return fmt.Errorf("insert task %q at index %d: %w", task.ID, index, err)
		}
	}

	if err := tx.Commit(); err != nil {
		return fmt.Errorf("commit tasks transaction: %w", err)
	}
	committed = true
	return nil
}

func QueryTasks(ctx context.Context, store *Store, query domain.TaskQuery) ([]domain.Task, error) {
	if err := validateStore(ctx, store, "query tasks"); err != nil {
		return nil, err
	}

	rows, err := store.DB.QueryContext(ctx, `
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
		return nil, fmt.Errorf("query tasks: %w", err)
	}
	defer rows.Close()

	taskList := make([]domain.Task, 0)
	for rows.Next() {
		var task domain.Task
		var status string
		var exitCode sql.NullInt64
		if err := rows.Scan(
			&task.RowOrder,
			&task.ID,
			&status,
			&task.Process,
			&task.Name,
			&task.Tag,
			&task.Workdir,
			&exitCode,
			&task.Duration,
			&task.Realtime,
			&task.CPUs,
			&task.Memory,
			&task.ErrorSummary,
		); err != nil {
			return nil, fmt.Errorf("scan task row: %w", err)
		}
		task.Status = domain.TaskStatus(status)
		if exitCode.Valid {
			value := int(exitCode.Int64)
			task.Exit = &value
		}
		taskList = append(taskList, task)
	}
	if err := rows.Err(); err != nil {
		return nil, fmt.Errorf("iterate task rows: %w", err)
	}

	filtered, err := tasks.ApplyTaskQuery(taskList, query)
	if err != nil {
		return nil, fmt.Errorf("apply task query: %w", err)
	}
	return filtered, nil
}

func CountTasksByStatus(ctx context.Context, store *Store) ([]domain.StatusCount, error) {
	if err := validateStore(ctx, store, "count tasks by status"); err != nil {
		return nil, err
	}

	rows, err := store.DB.QueryContext(ctx, `
		SELECT status, COUNT(*)
		FROM tasks
		GROUP BY status
		ORDER BY status
	`)
	if err != nil {
		return nil, fmt.Errorf("query task status counts: %w", err)
	}
	defer rows.Close()

	counts := make([]domain.StatusCount, 0)
	for rows.Next() {
		var status string
		var count int
		if err := rows.Scan(&status, &count); err != nil {
			return nil, fmt.Errorf("scan task status count: %w", err)
		}
		counts = append(counts, domain.StatusCount{
			Status: domain.TaskStatus(status),
			Count:  count,
		})
	}
	if err := rows.Err(); err != nil {
		return nil, fmt.Errorf("iterate task status counts: %w", err)
	}

	return counts, nil
}

func IndexDiagnostics(ctx context.Context, runDir domain.RunDir, artifacts domain.ArtifactSet) (domain.IndexDiagnostics, error) {
	if ctx == nil {
		return domain.IndexDiagnostics{}, fmt.Errorf("index diagnostics: nil context")
	}
	if err := ctx.Err(); err != nil {
		return domain.IndexDiagnostics{}, fmt.Errorf("index diagnostics: %w", err)
	}

	effectiveRunDir := runDir
	if effectiveRunDir.Path == "" {
		effectiveRunDir = artifacts.RunDir
	}
	if artifacts.RunDir.Path == "" {
		artifacts.RunDir = effectiveRunDir
	}

	diagnostics := domain.IndexDiagnostics{
		RunDir:    effectiveRunDir,
		Artifacts: artifacts,
	}

	metadataFor := func(freshness domain.IndexFreshness, reason string, indexPath string) *domain.IndexMetadata {
		metadata := metadataForArtifacts(effectiveRunDir.Path, indexPath, artifacts)
		metadata.Freshness = freshness
		metadata.StaleReason = reason
		return &metadata
	}

	indexDiagnostic := func(freshness domain.IndexFreshness, reason string) []domain.Diagnostic {
		switch freshness {
		case domain.IndexFreshnessMissing:
			if reason == "" {
				reason = "index metadata missing"
			}
			return []domain.Diagnostic{{
				Severity: domain.DiagnosticWarning,
				Code:     "index_missing",
				Message:  "index is missing",
				Detail:   "Run `gosh index --refresh` to create the index.",
			}}
		case domain.IndexFreshnessStale:
			return []domain.Diagnostic{{
				Severity: domain.DiagnosticWarning,
				Code:     "index_stale",
				Message:  "index is stale",
				Detail:   "Refresh with `gosh index --refresh`.",
			}}
		case domain.IndexFreshnessUnknown:
			if reason == "" {
				reason = "freshness could not be determined"
			}
			return []domain.Diagnostic{{
				Severity: domain.DiagnosticWarning,
				Code:     "index_freshness_unknown",
				Message:  "index freshness is unknown",
				Detail:   reason,
			}}
		default:
			return nil
		}
	}

	if artifacts.Mode == domain.IndexModeUnsupported {
		diagnostics.Metadata = metadataFor(domain.IndexFreshnessUnsupported, "no supported artifacts", "")
		return diagnostics, nil
	}

	if effectiveRunDir.Path == "" {
		return diagnostics, fmt.Errorf("index diagnostics: empty run dir")
	}

	indexPath := run.IndexPath(effectiveRunDir)
	if _, err := os.Stat(indexPath); err != nil {
		if os.IsNotExist(err) {
			reason := "index metadata missing"
			diagnostics.Metadata = metadataFor(domain.IndexFreshnessMissing, reason, indexPath)
			diagnostics.Diagnostics = indexDiagnostic(domain.IndexFreshnessMissing, reason)
			return diagnostics, nil
		}
		return diagnostics, fmt.Errorf("index diagnostics stat index %q: %w", indexPath, err)
	}

	store, err := OpenStore(ctx, effectiveRunDir)
	if err != nil {
		return diagnostics, fmt.Errorf("index diagnostics open store: %w", err)
	}
	defer store.Close()

	freshness, reason, err := CheckFreshness(ctx, store, artifacts)
	if err != nil {
		return diagnostics, fmt.Errorf("index diagnostics check freshness: %w", err)
	}

	metadata, err := ReadMetadata(ctx, store)
	if err != nil {
		return diagnostics, fmt.Errorf("index diagnostics read metadata: %w", err)
	}
	if metadata == nil {
		if reason == "" {
			reason = "index metadata missing"
		}
		diagnostics.Metadata = metadataFor(freshness, reason, store.Path)
		diagnostics.Diagnostics = indexDiagnostic(freshness, reason)
		return diagnostics, nil
	}

	metadataCopy := *metadata
	if metadataCopy.RunDir == "" {
		metadataCopy.RunDir = effectiveRunDir.Path
	}
	if metadataCopy.IndexPath == "" {
		metadataCopy.IndexPath = store.Path
	}
	metadataCopy.Freshness = freshness
	metadataCopy.StaleReason = reason
	if freshness == domain.IndexFreshnessFresh {
		metadataCopy.StaleReason = ""
	}
	diagnostics.Metadata = &metadataCopy
	diagnostics.Diagnostics = indexDiagnostic(freshness, metadataCopy.StaleReason)

	return diagnostics, nil
}
