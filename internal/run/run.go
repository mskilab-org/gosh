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
	DefaultRunDir         = "."
	DefaultResultsDirName = "results"
	PipelineInfoDirName   = "pipeline_info"
	IndexDirName          = ".gosh"
	IndexFileName         = "index.sqlite"
)

var TracePatterns = []string{"trace*.txt", "trace*.csv", "trace*.tsv"}
var PipelineInfoTracePatterns = []string{"execution_trace*.txt", "execution_trace*.tsv", "execution_trace*.csv"}
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

func ResolveResultsDir(runDir domain.RunDir, input string) (domain.ResultsDir, error) {
	if runDir.Path == "" {
		return domain.ResultsDir{}, fmt.Errorf("resolve results dir: empty run dir")
	}

	path := input
	if path == "" {
		path = filepath.Join(runDir.Path, DefaultResultsDirName)
	} else if !filepath.IsAbs(path) {
		path = filepath.Join(runDir.Path, path)
	}

	return domain.ResultsDir{Path: filepath.Clean(path)}, nil
}

func BuildArtifactSearchLocations(runDir domain.RunDir, resultsDir domain.ResultsDir) []domain.ArtifactSearchLocation {
	copyPatterns := func(patterns []string) []string {
		copied := make([]string, len(patterns))
		copy(copied, patterns)
		return copied
	}

	return []domain.ArtifactSearchLocation{
		{
			Kind:        domain.SourceKindTrace,
			BaseDir:     runDir.Path,
			Patterns:    copyPatterns(TracePatterns),
			Description: "run directory trace files",
		},
		{
			Kind:        domain.SourceKindTrace,
			BaseDir:     filepath.Join(resultsDir.Path, PipelineInfoDirName),
			Patterns:    copyPatterns(PipelineInfoTracePatterns),
			Description: "pipeline_info execution trace files",
		},
		{
			Kind:        domain.SourceKindLog,
			BaseDir:     runDir.Path,
			Patterns:    copyPatterns(LogPatterns),
			Description: "run directory log files",
		},
	}
}

type sourceFingerprintCollector struct {
	sources   []domain.SourceFingerprint
	seenPaths map[string]struct{}
}

func newSourceFingerprintCollector() sourceFingerprintCollector {
	return sourceFingerprintCollector{
		sources:   make([]domain.SourceFingerprint, 0),
		seenPaths: make(map[string]struct{}),
	}
}

func (collector *sourceFingerprintCollector) add(source domain.SourceFingerprint) {
	if _, seen := collector.seenPaths[source.Path]; seen {
		return
	}
	collector.seenPaths[source.Path] = struct{}{}
	collector.sources = append(collector.sources, source)
}

func (collector *sourceFingerprintCollector) addAll(sources []domain.SourceFingerprint) {
	for _, source := range sources {
		collector.add(source)
	}
}

func (collector *sourceFingerprintCollector) sorted() []domain.SourceFingerprint {
	sort.Slice(collector.sources, func(i, j int) bool {
		return collector.sources[i].Path < collector.sources[j].Path
	})
	return collector.sources
}

func FindCandidateSourcesInLocations(locations []domain.ArtifactSearchLocation) ([]domain.SourceFingerprint, error) {
	collector := newSourceFingerprintCollector()

	for _, location := range locations {
		locationSources, err := FindCandidateSources(domain.RunDir{Path: location.BaseDir}, location.Kind, location.Patterns)
		if err != nil {
			description := location.Description
			if description == "" {
				description = location.BaseDir
			}
			return nil, fmt.Errorf("find candidate sources in %q: %w", description, err)
		}

		collector.addAll(locationSources)
	}

	return collector.sorted(), nil
}

func DiscoverArtifactsWithResultsDir(ctx context.Context, runDir domain.RunDir, resultsDir domain.ResultsDir) (domain.ArtifactSet, error) {
	searchLocations := BuildArtifactSearchLocations(runDir, resultsDir)
	searchedPatterns := make([]string, 0)
	for _, location := range searchLocations {
		searchedPatterns = append(searchedPatterns, location.Patterns...)
	}

	artifacts := domain.ArtifactSet{
		RunDir:           runDir,
		ResultsDir:       resultsDir,
		SelectedAt:       time.Now().UTC(),
		SearchedPatterns: searchedPatterns,
		SearchLocations:  searchLocations,
	}

	if ctx == nil {
		return artifacts, fmt.Errorf("discover artifacts: nil context")
	}
	if err := ctx.Err(); err != nil {
		return artifacts, err
	}

	sources, err := FindCandidateSourcesInLocations(searchLocations)
	if err != nil {
		return artifacts, fmt.Errorf("discover artifact sources: %w", err)
	}
	if err := ctx.Err(); err != nil {
		return artifacts, err
	}

	traceSources := make([]domain.SourceFingerprint, 0)
	logSources := make([]domain.SourceFingerprint, 0)
	for _, source := range sources {
		switch source.Kind {
		case domain.SourceKindTrace:
			traceSources = append(traceSources, source)
		case domain.SourceKindLog:
			logSources = append(logSources, source)
		}
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
	resultsDir, err := ResolveResultsDir(runDir, "")
	if err != nil {
		return domain.ArtifactSet{RunDir: runDir}, err
	}

	return DiscoverArtifactsWithResultsDir(ctx, runDir, resultsDir)
}

func FindCandidateSources(runDir domain.RunDir, kind domain.SourceKind, patterns []string) ([]domain.SourceFingerprint, error) {
	collector := newSourceFingerprintCollector()

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
			collector.add(fingerprint)
		}
	}

	return collector.sorted(), nil
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
	formatLocations := func(locations []domain.ArtifactSearchLocation) string {
		if len(locations) == 0 {
			return ""
		}

		formatted := ""
		for i, location := range locations {
			if i > 0 {
				formatted += "\n"
			}
			formatted += "- " + location.BaseDir
			if location.Description != "" {
				formatted += " (" + location.Description + ")"
			}
			formatted += ": " + joinPatterns(location.Patterns)
		}
		return formatted
	}

	resultsDir, err := ResolveResultsDir(runDir, "")
	if err != nil {
		resultsDir = domain.ResultsDir{Path: filepath.Join(runDir.Path, DefaultResultsDirName)}
	}
	locations := BuildArtifactSearchLocations(runDir, resultsDir)
	traceLocations := make([]domain.ArtifactSearchLocation, 0)
	logLocations := make([]domain.ArtifactSearchLocation, 0)
	for _, location := range locations {
		switch location.Kind {
		case domain.SourceKindTrace:
			traceLocations = append(traceLocations, location)
		case domain.SourceKindLog:
			logLocations = append(logLocations, location)
		}
	}

	return []domain.Diagnostic{
		{
			Severity: domain.DiagnosticError,
			Code:     "unsupported_artifacts",
			Message:  "No supported Nextflow trace or log artifacts found in " + runDir.Path,
			Detail: "Searched trace locations:\n" + formatLocations(traceLocations) + "\n" +
				"Searched log locations:\n" + formatLocations(logLocations),
		},
		domain.NextflowTraceRecommendationDiagnostic(),
	}
}
