package run

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"time"

	"github.com/mskilab-org/gosh/internal/domain"
)

const (
	DefaultRunDir = "."
	IndexDirName  = ".gosh"
	IndexFileName = "index.sqlite"
)

var TracePatterns = []string{"trace*.txt", "trace*.csv", "trace*.tsv"}
var LogPatterns = []string{".nextflow.log", ".nextflow_*.log"}

func ResolveRunDir(input string) (domain.RunDir, error) {
	if input == "" {
		input = DefaultRunDir
	}

	path, err := filepath.Abs(input)
	if err != nil {
		return domain.RunDir{}, fmt.Errorf("resolve run dir %q: %w", input, err)
	}

	info, err := os.Stat(path)
	if err != nil {
		return domain.RunDir{}, fmt.Errorf("resolve run dir %q: %w", input, err)
	}
	if !info.IsDir() {
		return domain.RunDir{}, fmt.Errorf("resolve run dir %q: not a directory", input)
	}

	return domain.RunDir{Path: filepath.Clean(path)}, nil
}

func IndexPath(runDir domain.RunDir) string {
	return filepath.Join(runDir.Path, IndexDirName, IndexFileName)
}

func chooseArtifactSource(label string, sources []domain.SourceFingerprint) (*domain.SourceFingerprint, error) {
	source, err := ChooseNewestSource(sources)
	if err != nil {
		return nil, fmt.Errorf("choose %s artifact: %w", label, err)
	}
	return source, nil
}

func chooseOptionalArtifactSource(label string, sources []domain.SourceFingerprint) (*domain.SourceFingerprint, error) {
	if len(sources) == 0 {
		return nil, nil
	}
	return chooseArtifactSource(label, sources)
}

func DiscoverArtifacts(ctx context.Context, runDir domain.RunDir) (domain.ArtifactSet, error) {
	searchedPatterns := append([]string{}, TracePatterns...)
	searchedPatterns = append(searchedPatterns, LogPatterns...)

	artifacts := domain.ArtifactSet{
		RunDir:           runDir,
		SelectedAt:       time.Now().UTC(),
		SearchedPatterns: searchedPatterns,
	}

	if ctx == nil {
		return artifacts, fmt.Errorf("discover artifacts: nil context")
	}
	if err := ctx.Err(); err != nil {
		return artifacts, err
	}

	traceSources, err := FindCandidateSources(runDir, domain.SourceKindTrace, TracePatterns)
	if err != nil {
		return artifacts, fmt.Errorf("discover trace artifacts: %w", err)
	}
	if err := ctx.Err(); err != nil {
		return artifacts, err
	}

	logSources, err := FindCandidateSources(runDir, domain.SourceKindLog, LogPatterns)
	if err != nil {
		return artifacts, fmt.Errorf("discover log artifacts: %w", err)
	}
	if err := ctx.Err(); err != nil {
		return artifacts, err
	}

	log, err := chooseOptionalArtifactSource("log", logSources)
	if err != nil {
		return artifacts, err
	}

	if len(traceSources) > 0 {
		trace, err := chooseArtifactSource("trace", traceSources)
		if err != nil {
			return artifacts, err
		}
		artifacts.Mode = domain.IndexModeTraceBacked
		artifacts.Trace = trace
		artifacts.Log = log
		return artifacts, nil
	}

	if log != nil {
		artifacts.Mode = domain.IndexModeLogOnly
		artifacts.Log = log
		return artifacts, nil
	}

	artifacts.Mode = domain.IndexModeUnsupported
	artifacts.Diagnostics = UnsupportedArtifactDiagnostics(runDir)
	return artifacts, nil
}

func FindCandidateSources(runDir domain.RunDir, kind domain.SourceKind, patterns []string) ([]domain.SourceFingerprint, error) {
	sources := make([]domain.SourceFingerprint, 0)
	seenPaths := make(map[string]struct{})

	for _, pattern := range patterns {
		globPattern := filepath.Join(runDir.Path, pattern)
		matches, err := filepath.Glob(globPattern)
		if err != nil {
			return nil, fmt.Errorf("glob source pattern %q: %w", pattern, err)
		}

		for _, match := range matches {
			fingerprint, err := SourceFingerprintForPath(kind, match)
			if err != nil {
				info, statErr := os.Stat(match)
				if statErr == nil && !info.Mode().IsRegular() {
					continue
				}
				return nil, fmt.Errorf("fingerprint source %q: %w", match, err)
			}
			if _, seen := seenPaths[fingerprint.Path]; seen {
				continue
			}
			seenPaths[fingerprint.Path] = struct{}{}
			sources = append(sources, fingerprint)
		}
	}

	sort.Slice(sources, func(i, j int) bool {
		return sources[i].Path < sources[j].Path
	})
	return sources, nil
}

func SourceFingerprintForPath(kind domain.SourceKind, path string) (domain.SourceFingerprint, error) {
	if path == "" {
		return domain.SourceFingerprint{}, fmt.Errorf("stat source: empty path")
	}

	resolvedPath, err := filepath.Abs(path)
	if err != nil {
		return domain.SourceFingerprint{}, fmt.Errorf("resolve source %q: %w", path, err)
	}
	resolvedPath = filepath.Clean(resolvedPath)

	info, err := os.Stat(resolvedPath)
	if err != nil {
		return domain.SourceFingerprint{}, fmt.Errorf("stat source %q: %w", path, err)
	}
	if !info.Mode().IsRegular() {
		return domain.SourceFingerprint{}, fmt.Errorf("stat source %q: not a regular file", path)
	}

	return domain.SourceFingerprint{
		Kind:    kind,
		Path:    resolvedPath,
		ModTime: info.ModTime(),
		Size:    info.Size(),
	}, nil
}

func ChooseNewestSource(sources []domain.SourceFingerprint) (*domain.SourceFingerprint, error) {
	if len(sources) == 0 {
		return nil, fmt.Errorf("choose newest source: no candidate sources")
	}

	selectedIndex := 0
	for i := 1; i < len(sources); i++ {
		candidate := sources[i]
		selected := sources[selectedIndex]
		if candidate.ModTime.After(selected.ModTime) || (candidate.ModTime.Equal(selected.ModTime) && candidate.Path < selected.Path) {
			selectedIndex = i
		}
	}

	selected := sources[selectedIndex]
	return &selected, nil
}

func UnsupportedArtifactDiagnostics(runDir domain.RunDir) []domain.Diagnostic {
	joinPatterns := func(patterns []string) string {
		if len(patterns) == 0 {
			return ""
		}

		joined := patterns[0]
		for _, pattern := range patterns[1:] {
			joined += ", " + pattern
		}
		return joined
	}

	return []domain.Diagnostic{
		{
			Severity: domain.DiagnosticError,
			Code:     "unsupported_artifacts",
			Message:  "No supported Nextflow trace or log artifacts found in " + runDir.Path,
			Detail: "Searched trace patterns: " + joinPatterns(TracePatterns) + "\n" +
				"Searched log patterns: " + joinPatterns(LogPatterns),
		},
		{
			Severity: domain.DiagnosticInfo,
			Code:     "nextflow_with_trace_recommended",
			Message:  "Run future Nextflow workflows with -with-trace",
			Detail:   "Use `nextflow run ... -with-trace` for future runs so gosh can build a complete trace-backed task index.",
		},
	}
}
