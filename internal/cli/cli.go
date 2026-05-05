package cli

import (
	"context"
	"fmt"
	"io"
	"strings"

	"github.com/mskilab-org/gosh/internal/domain"
	idx "github.com/mskilab-org/gosh/internal/index"
	"github.com/mskilab-org/gosh/internal/inspect"
	"github.com/mskilab-org/gosh/internal/nflog"
	"github.com/mskilab-org/gosh/internal/render"
	"github.com/mskilab-org/gosh/internal/run"
	"github.com/mskilab-org/gosh/internal/tasks"
)

type CommandName string

const (
	CommandStatus  CommandName = "status"
	CommandTasks   CommandName = "tasks"
	CommandInspect CommandName = "inspect"
	CommandIndex   CommandName = "index"
)

type CommandSpec struct {
	Name    CommandName
	Usage   string
	Purpose string
}

type RootCommand struct {
	Version       string
	DefaultRunDir string
	Commands      []CommandSpec
}

type GlobalOptions struct {
	RunDir     string
	ResultsDir string
	Format     domain.OutputFormat
}

type StatusOptions struct {
	Global GlobalOptions
}

type TasksOptions struct {
	Global GlobalOptions
	Query  domain.TaskQuery
}

type InspectOptions struct {
	Global   GlobalOptions
	Selector string
}

type IndexOptions struct {
	Global  GlobalOptions
	Refresh bool
}

func NewRootCommand(version string) RootCommand {
	return RootCommand{
		Version:       version,
		DefaultRunDir: run.DefaultRunDir,
		Commands: []CommandSpec{
			{Name: CommandStatus, Usage: "gosh status [--run-dir DIR] [--results-dir DIR] [--json]", Purpose: "summarize the selected Nextflow run"},
			{Name: CommandTasks, Usage: "gosh tasks [--run-dir DIR] [--results-dir DIR] [filters] [--json]", Purpose: "list trace-backed task rows"},
			{Name: CommandInspect, Usage: "gosh inspect <task> [--run-dir DIR] [--results-dir DIR] [--json]", Purpose: "show a bounded task dossier"},
			{Name: CommandIndex, Usage: "gosh index [--run-dir DIR] [--results-dir DIR] [--refresh] [--json]", Purpose: "show or rebuild index diagnostics"},
		},
	}
}

func Main(args []string, stdout io.Writer, stderr io.Writer, version string) int {
	if err := Execute(context.Background(), args, stdout, stderr, version); err != nil {
		_, _ = fmt.Fprintln(stderr, err)
		return 1
	}
	return 0
}

func Execute(ctx context.Context, args []string, stdout io.Writer, stderr io.Writer, version string) error {
	_ = stderr
	if stdout == nil {
		stdout = io.Discard
	}

	root := NewRootCommand(version)

	writeVersion := func() error {
		_, err := fmt.Fprintf(stdout, "gosh %s\n", root.Version)
		return err
	}

	writeRootHelp := func() error {
		header := "gosh"
		if root.Version != "" {
			header += " " + root.Version
		}
		if _, err := fmt.Fprintf(stdout, "%s\n\n", header); err != nil {
			return err
		}
		if _, err := fmt.Fprintln(stdout, "Usage:"); err != nil {
			return err
		}
		if _, err := fmt.Fprintln(stdout, "  gosh [global options] <command> [command options]"); err != nil {
			return err
		}
		if _, err := fmt.Fprintln(stdout); err != nil {
			return err
		}
		if _, err := fmt.Fprintln(stdout, "Global options:"); err != nil {
			return err
		}
		globalLines := []string{
			fmt.Sprintf("  --run-dir, -d DIR   Nextflow run directory (default %q)", root.DefaultRunDir),
			"  --results-dir DIR   Nextflow results directory (default <run-dir>/results)",
			"  --json              Write JSON output",
			"  --format FORMAT     Output format: human or json",
			"  --help, -h          Show help",
			"  --version           Show version",
		}
		for _, line := range globalLines {
			if _, err := fmt.Fprintln(stdout, line); err != nil {
				return err
			}
		}
		if _, err := fmt.Fprintln(stdout); err != nil {
			return err
		}
		if _, err := fmt.Fprintln(stdout, "Commands:"); err != nil {
			return err
		}
		for _, command := range root.Commands {
			if _, err := fmt.Fprintf(stdout, "  %-46s %s\n", command.Usage, command.Purpose); err != nil {
				return err
			}
		}
		if _, err := fmt.Fprintln(stdout); err != nil {
			return err
		}
		_, err := fmt.Fprintln(stdout, "Use \"gosh <command> --help\" for command-specific options.")
		return err
	}

	findCommand := func(name CommandName) (CommandSpec, bool) {
		for _, command := range root.Commands {
			if command.Name == name {
				return command, true
			}
		}
		return CommandSpec{}, false
	}

	writeCommandHelp := func(name CommandName) error {
		command, ok := findCommand(name)
		if !ok {
			return fmt.Errorf("unknown command %q", string(name))
		}

		if _, err := fmt.Fprintln(stdout, "Usage:"); err != nil {
			return err
		}
		if _, err := fmt.Fprintf(stdout, "  %s\n", command.Usage); err != nil {
			return err
		}
		if _, err := fmt.Fprintln(stdout); err != nil {
			return err
		}
		if _, err := fmt.Fprintln(stdout, "Purpose:"); err != nil {
			return err
		}
		if _, err := fmt.Fprintf(stdout, "  %s\n", command.Purpose); err != nil {
			return err
		}
		if _, err := fmt.Fprintln(stdout); err != nil {
			return err
		}
		if _, err := fmt.Fprintln(stdout, "Options:"); err != nil {
			return err
		}

		optionLines := []string{
			fmt.Sprintf("  --run-dir, -d DIR   Nextflow run directory (default %q)", root.DefaultRunDir),
			"  --results-dir DIR   Nextflow results directory (default <run-dir>/results)",
			"  --json              Write JSON output",
			"  --format FORMAT     Output format: human or json",
			"  --help, -h          Show command help",
		}
		switch name {
		case CommandTasks:
			optionLines = append(optionLines,
				"  --process VALUE     Filter by process substring",
				"  --name VALUE        Filter by task name substring",
				"  --sample VALUE      Filter by sample/tag substring",
				"  --status STATUS     Filter by task status",
			)
		case CommandInspect:
			optionLines = append(optionLines, "  <task>              Required canonical ID, workdir path, or unambiguous task selector")
		case CommandIndex:
			optionLines = append(optionLines, "  --refresh           Refresh run-local index metadata")
		}
		for _, line := range optionLines {
			if _, err := fmt.Fprintln(stdout, line); err != nil {
				return err
			}
		}
		return nil
	}

	for i := 0; i < len(args); i++ {
		arg := args[i]
		switch {
		case arg == "--help" || arg == "-h":
			return writeRootHelp()
		case arg == "--version":
			return writeVersion()
		case arg == "--":
			i = len(args)
		case arg == "--run-dir" || arg == "-d" || arg == "--results-dir" || arg == "--format":
			i++
		case arg == "--json" || strings.HasPrefix(arg, "--run-dir=") || strings.HasPrefix(arg, "--results-dir=") || strings.HasPrefix(arg, "--format="):
			// Continue scanning global options before the command.
		case strings.HasPrefix(arg, "-"):
			i = len(args)
		default:
			i = len(args)
		}
	}

	global, remaining, err := ParseGlobalOptions(args)
	if err != nil {
		return err
	}
	if len(remaining) == 0 {
		return writeRootHelp()
	}

	commandName := CommandName(remaining[0])
	commandArgs := remaining[1:]
	wantsCommandHelp := false
	for _, arg := range commandArgs {
		if arg == "--" {
			break
		}
		if arg == "--help" || arg == "-h" {
			wantsCommandHelp = true
			break
		}
	}

	switch commandName {
	case CommandStatus:
		if wantsCommandHelp {
			return writeCommandHelp(commandName)
		}
		options, err := ParseStatusOptions(global, commandArgs)
		if err != nil {
			return err
		}
		return RunStatus(ctx, options, stdout)

	case CommandTasks:
		if wantsCommandHelp {
			return writeCommandHelp(commandName)
		}
		options, err := ParseTasksOptions(global, commandArgs)
		if err != nil {
			return err
		}
		return RunTasks(ctx, options, stdout)

	case CommandInspect:
		if wantsCommandHelp {
			return writeCommandHelp(commandName)
		}
		options, err := ParseInspectOptions(global, commandArgs)
		if err != nil {
			return err
		}
		return RunInspect(ctx, options, stdout)

	case CommandIndex:
		if wantsCommandHelp {
			return writeCommandHelp(commandName)
		}
		options, err := ParseIndexOptions(global, commandArgs)
		if err != nil {
			return err
		}
		return RunIndex(ctx, options, stdout)

	default:
		return fmt.Errorf("unknown command %q", remaining[0])
	}
}

func parseOptionValue(args []string, index int, option string, valueName string) (string, bool, int, error) {
	arg := args[index]
	if arg == option {
		if index+1 >= len(args) || args[index+1] == "" {
			return "", true, index, fmt.Errorf("%s requires %s", option, valueName)
		}
		return args[index+1], true, index + 1, nil
	}

	prefix := option + "="
	if strings.HasPrefix(arg, prefix) {
		value := arg[len(prefix):]
		if value == "" {
			return "", true, index, fmt.Errorf("%s requires %s", option, valueName)
		}
		return value, true, index, nil
	}

	return "", false, index, nil
}

func parseCLIOutputFormat(value string) (domain.OutputFormat, error) {
	switch value {
	case string(domain.OutputFormatHuman):
		return domain.OutputFormatHuman, nil
	case string(domain.OutputFormatJSON):
		return domain.OutputFormatJSON, nil
	default:
		return "", fmt.Errorf("unsupported output format %q", value)
	}
}

func parseSharedGlobalOption(args []string, index int, options *GlobalOptions) (bool, int, error) {
	if options == nil {
		return false, index, fmt.Errorf("parse global option: nil options")
	}

	arg := args[index]
	if value, handled, next, err := parseOptionValue(args, index, "--run-dir", "DIR"); handled {
		if err != nil {
			return true, index, err
		}
		options.RunDir = value
		return true, next, nil
	}

	if arg == "-d" {
		if index+1 >= len(args) || args[index+1] == "" {
			return true, index, fmt.Errorf("-d requires DIR")
		}
		options.RunDir = args[index+1]
		return true, index + 1, nil
	}

	if value, handled, next, err := parseOptionValue(args, index, "--results-dir", "DIR"); handled {
		if err != nil {
			return true, index, err
		}
		options.ResultsDir = value
		return true, next, nil
	}

	if arg == "--json" {
		options.Format = domain.OutputFormatJSON
		return true, index, nil
	}

	if value, handled, next, err := parseOptionValue(args, index, "--format", "FORMAT"); handled {
		if err != nil {
			return true, index, err
		}
		format, err := parseCLIOutputFormat(value)
		if err != nil {
			return true, index, err
		}
		options.Format = format
		return true, next, nil
	}

	return false, index, nil
}

func normalizeOutputFormat(format domain.OutputFormat) (domain.OutputFormat, error) {
	if format == "" {
		return domain.OutputFormatHuman, nil
	}
	switch format {
	case domain.OutputFormatHuman, domain.OutputFormatJSON:
		return format, nil
	default:
		return "", fmt.Errorf("unsupported output format %q", format)
	}
}

func nextflowTraceRecommendation() domain.Diagnostic {
	return domain.NextflowTraceRecommendationDiagnostic()
}

func selectedLogPath(artifacts domain.ArtifactSet) string {
	if artifacts.Log != nil && artifacts.Log.Path != "" {
		return artifacts.Log.Path
	}
	return "selected Nextflow log"
}

func diagnosticsWithTraceRecommendation(artifacts domain.ArtifactSet, diagnostic domain.Diagnostic) []domain.Diagnostic {
	diagnostics := append([]domain.Diagnostic{}, artifacts.Diagnostics...)
	return append(diagnostics, diagnostic, nextflowTraceRecommendation())
}

func unsupportedArtifactDiagnosticsForCommand(runDir domain.RunDir, artifacts domain.ArtifactSet) []domain.Diagnostic {
	if len(artifacts.Diagnostics) > 0 {
		diagnostics := make([]domain.Diagnostic, len(artifacts.Diagnostics))
		copy(diagnostics, artifacts.Diagnostics)
		return diagnostics
	}
	if len(artifacts.SearchLocations) > 0 {
		return run.UnsupportedArtifactDiagnosticsFromSearchLocations(runDir, artifacts.SearchLocations)
	}
	return run.UnsupportedArtifactDiagnostics(runDir)
}

func diagnosticsForUnsupportedArtifacts(runDir domain.RunDir, artifacts domain.ArtifactSet, diagnostic domain.Diagnostic) []domain.Diagnostic {
	diagnostics := append([]domain.Diagnostic{}, artifacts.Diagnostics...)
	if len(diagnostics) == 0 {
		diagnostics = run.UnsupportedArtifactDiagnostics(runDir)
	}
	return append(diagnostics, diagnostic)
}

func renderUnsupportedCommand(writer io.Writer, diagnostics []domain.Diagnostic, format domain.OutputFormat, commandErr error) error {
	if err := render.RenderUnsupportedDiagnostics(writer, diagnostics, format); err != nil {
		return err
	}
	return commandErr
}

type commandContext struct {
	RunDir     domain.RunDir
	ResultsDir domain.ResultsDir
	Artifacts  domain.ArtifactSet
	Format     domain.OutputFormat
}

func loadCommandContext(ctx context.Context, global GlobalOptions) (commandContext, error) {
	runDir, err := run.ResolveRunDir(global.RunDir)
	if err != nil {
		return commandContext{}, err
	}

	artifacts, err := run.DiscoverArtifacts(ctx, runDir)
	if err != nil {
		return commandContext{}, err
	}

	format, err := normalizeOutputFormat(global.Format)
	if err != nil {
		return commandContext{}, err
	}

	return commandContext{RunDir: runDir, Artifacts: artifacts, Format: format}, nil
}

func loadCommandContextWithResultsDir(ctx context.Context, global GlobalOptions) (commandContext, error) {
	runDir, err := run.ResolveRunDir(global.RunDir)
	if err != nil {
		return commandContext{}, err
	}

	resultsDir, err := run.ResolveResultsDir(runDir, global.ResultsDir)
	if err != nil {
		return commandContext{}, err
	}

	artifacts, err := run.DiscoverArtifactsWithResultsDir(ctx, runDir, resultsDir)
	if err != nil {
		return commandContext{}, err
	}

	format, err := normalizeOutputFormat(global.Format)
	if err != nil {
		return commandContext{}, err
	}

	return commandContext{RunDir: runDir, ResultsDir: resultsDir, Artifacts: artifacts, Format: format}, nil
}

func ParseGlobalOptions(args []string) (GlobalOptions, []string, error) {
	options := GlobalOptions{RunDir: run.DefaultRunDir, Format: domain.OutputFormatHuman}

	for i := 0; i < len(args); i++ {
		arg := args[i]
		if arg == "--" {
			return options, args[i+1:], nil
		}

		handled, next, err := parseSharedGlobalOption(args, i, &options)
		if err != nil {
			return GlobalOptions{}, nil, err
		}
		if handled {
			i = next
			continue
		}

		if len(arg) > 0 && arg[0] == '-' {
			return GlobalOptions{}, nil, fmt.Errorf("unknown global option %q", arg)
		}
		return options, args[i:], nil
	}

	return options, args[len(args):], nil
}

func ParseStatusOptions(global GlobalOptions, args []string) (StatusOptions, error) {
	options := StatusOptions{Global: global}

	for i := 0; i < len(args); i++ {
		arg := args[i]
		if arg == "--" {
			if i+1 < len(args) {
				return StatusOptions{}, fmt.Errorf("unexpected status argument %q", args[i+1])
			}
			return options, nil
		}

		handled, next, err := parseSharedGlobalOption(args, i, &options.Global)
		if err != nil {
			return StatusOptions{}, err
		}
		if handled {
			i = next
			continue
		}

		if len(arg) > 0 && arg[0] == '-' {
			return StatusOptions{}, fmt.Errorf("unknown status option %q", arg)
		}
		return StatusOptions{}, fmt.Errorf("unexpected status argument %q", arg)
	}

	return options, nil
}

func ParseTasksOptions(global GlobalOptions, args []string) (TasksOptions, error) {
	options := TasksOptions{Global: global}

	for i := 0; i < len(args); i++ {
		arg := args[i]
		if arg == "--" {
			if i+1 < len(args) {
				return TasksOptions{}, fmt.Errorf("unexpected tasks argument %q", args[i+1])
			}
			return options, nil
		}

		handled, next, err := parseSharedGlobalOption(args, i, &options.Global)
		if err != nil {
			return TasksOptions{}, err
		}
		if handled {
			i = next
			continue
		}

		if value, handled, next, err := parseOptionValue(args, i, "--process", "VALUE"); handled {
			if err != nil {
				return TasksOptions{}, err
			}
			options.Query.ProcessSubstring = value
			i = next
			continue
		}
		if value, handled, next, err := parseOptionValue(args, i, "--name", "VALUE"); handled {
			if err != nil {
				return TasksOptions{}, err
			}
			options.Query.NameSubstring = value
			i = next
			continue
		}
		if value, handled, next, err := parseOptionValue(args, i, "--sample", "VALUE"); handled {
			if err != nil {
				return TasksOptions{}, err
			}
			options.Query.SampleSubstring = value
			i = next
			continue
		}
		if value, handled, next, err := parseOptionValue(args, i, "--status", "STATUS"); handled {
			if err != nil {
				return TasksOptions{}, err
			}
			options.Query.StatusRaw = value
			i = next
			continue
		}

		if len(arg) > 0 && arg[0] == '-' {
			return TasksOptions{}, fmt.Errorf("unknown tasks option %q", arg)
		}
		return TasksOptions{}, fmt.Errorf("unexpected tasks argument %q", arg)
	}

	return options, nil
}

func ParseInspectOptions(global GlobalOptions, args []string) (InspectOptions, error) {
	options := InspectOptions{Global: global}
	selectors := make([]string, 0, 1)

	for i := 0; i < len(args); i++ {
		arg := args[i]
		if arg == "--" {
			for _, selector := range args[i+1:] {
				if selector == "" {
					return InspectOptions{}, fmt.Errorf("inspect requires exactly one selector")
				}
				selectors = append(selectors, selector)
			}
			if len(selectors) != 1 {
				return InspectOptions{}, fmt.Errorf("inspect requires exactly one selector")
			}
			options.Selector = selectors[0]
			return options, nil
		}

		handled, next, err := parseSharedGlobalOption(args, i, &options.Global)
		if err != nil {
			return InspectOptions{}, err
		}
		if handled {
			i = next
			continue
		}

		if len(arg) > 0 && arg[0] == '-' {
			return InspectOptions{}, fmt.Errorf("unknown inspect option %q", arg)
		}
		if arg == "" {
			return InspectOptions{}, fmt.Errorf("inspect requires exactly one selector")
		}
		selectors = append(selectors, arg)
	}

	if len(selectors) != 1 {
		return InspectOptions{}, fmt.Errorf("inspect requires exactly one selector")
	}
	options.Selector = selectors[0]
	return options, nil
}

func ParseIndexOptions(global GlobalOptions, args []string) (IndexOptions, error) {
	options := IndexOptions{Global: global}

	for i := 0; i < len(args); i++ {
		arg := args[i]
		if arg == "--" {
			if i+1 < len(args) {
				return IndexOptions{}, fmt.Errorf("unexpected index argument %q", args[i+1])
			}
			return options, nil
		}

		handled, next, err := parseSharedGlobalOption(args, i, &options.Global)
		if err != nil {
			return IndexOptions{}, err
		}
		if handled {
			i = next
			continue
		}

		if arg == "--refresh" {
			options.Refresh = true
			continue
		}
		if len(arg) > 0 && arg[0] == '-' {
			return IndexOptions{}, fmt.Errorf("unknown index option %q", arg)
		}
		return IndexOptions{}, fmt.Errorf("unexpected index argument %q", arg)
	}

	return options, nil
}

func RunStatus(ctx context.Context, options StatusOptions, writer io.Writer) error {
	command, err := loadCommandContextWithResultsDir(ctx, options.Global)
	if err != nil {
		return err
	}
	runDir := command.RunDir
	artifacts := command.Artifacts
	format := command.Format

	renderStatus := func(summary domain.StatusSummary) error {
		view := domain.StatusView{Summary: summary, Format: format}
		switch format {
		case domain.OutputFormatJSON:
			return render.RenderStatusJSON(writer, view)
		case domain.OutputFormatHuman:
			return render.RenderStatusHuman(writer, view)
		default:
			return fmt.Errorf("unsupported output format %q", format)
		}
	}

	switch artifacts.Mode {
	case domain.IndexModeTraceBacked:
		store, metadata, err := idx.EnsureFreshIndex(ctx, runDir, artifacts)
		if err != nil {
			return err
		}
		defer store.Close()

		counts, err := idx.CountTasksByStatus(ctx, store)
		if err != nil {
			return err
		}

		failedTasks, err := idx.QueryTasks(ctx, store, domain.TaskQuery{Status: domain.TaskStatusFailed})
		if err != nil {
			return err
		}
		abortedTasks, err := idx.QueryTasks(ctx, store, domain.TaskQuery{Status: domain.TaskStatusAborted})
		if err != nil {
			return err
		}
		failedTasks = append(failedTasks, abortedTasks...)

		summary, err := tasks.BuildStatusSummary(runDir, metadata, counts, failedTasks)
		if err != nil {
			return err
		}
		return renderStatus(summary)

	case domain.IndexModeLogOnly:
		if artifacts.Log == nil {
			return fmt.Errorf("status log-only: missing log source")
		}

		evidence, err := nflog.ParseLogOnlyTaskEvidence(ctx, runDir, *artifacts.Log)
		if err != nil {
			return err
		}
		summary, err := nflog.BuildLogOnlyEvidenceStatus(runDir, artifacts, evidence)
		if err != nil {
			return err
		}
		return renderStatus(summary)

	case domain.IndexModeUnsupported:
		diagnostics := artifacts.Diagnostics
		if len(diagnostics) == 0 {
			diagnostics = run.UnsupportedArtifactDiagnostics(runDir)
		}
		if err := render.RenderUnsupportedDiagnostics(writer, diagnostics, format); err != nil {
			return err
		}
		return fmt.Errorf("status: no usable artifacts in %s", runDir.Path)

	default:
		return fmt.Errorf("status: unsupported artifact mode %q", artifacts.Mode)
	}
}

func RunTasks(ctx context.Context, options TasksOptions, writer io.Writer) error {
	command, err := loadCommandContextWithResultsDir(ctx, options.Global)
	if err != nil {
		return err
	}
	runDir := command.RunDir
	artifacts := command.Artifacts
	format := command.Format

	renderTasks := func(taskRows []domain.Task, metadata *domain.IndexMetadata) error {
		view := domain.TasksView{Tasks: taskRows, Query: options.Query, Metadata: metadata, Format: format}
		switch format {
		case domain.OutputFormatJSON:
			return render.RenderTasksJSON(writer, view)
		case domain.OutputFormatHuman:
			return render.RenderTasksHuman(writer, view)
		default:
			return fmt.Errorf("unsupported output format %q", format)
		}
	}

	switch artifacts.Mode {
	case domain.IndexModeTraceBacked:
		store, metadata, err := idx.EnsureFreshIndex(ctx, runDir, artifacts)
		if err != nil {
			return err
		}
		defer store.Close()

		taskRows, err := idx.QueryTasks(ctx, store, options.Query)
		if err != nil {
			return err
		}
		return renderTasks(taskRows, &metadata)

	case domain.IndexModeLogOnly:
		diagnostic := domain.Diagnostic{
			Severity: domain.DiagnosticError,
			Code:     "tasks_unavailable_log_only",
			Message:  "gosh tasks requires a trace-backed task index; complete task/resource/status data is unavailable in log-only mode",
			Detail:   "Selected log: " + selectedLogPath(artifacts) + "\nOnly deterministic log-only failure evidence may be available; complete task rows require a Nextflow trace file.",
		}
		if format == domain.OutputFormatHuman {
			diagnostic.Message = "task table unavailable without a trace file"
			diagnostic.Detail = strings.Join([]string{
				"Mode: " + string(artifacts.Mode),
				"Run dir: " + runDir.Path,
				"Selected log: " + selectedLogPath(artifacts),
				"Only deterministic log-only failure evidence may be available; complete task/resource/status data is unavailable in log-only mode.",
				"Complete task rows require a Nextflow trace file.",
			}, "\n")
		}
		diagnostics := diagnosticsWithTraceRecommendation(artifacts, diagnostic)
		return renderUnsupportedCommand(writer, diagnostics, format, fmt.Errorf("tasks: complete task/resource/status data is unavailable in log-only mode"))

	case domain.IndexModeUnsupported:
		diagnostics := diagnosticsForUnsupportedArtifacts(runDir, artifacts, domain.Diagnostic{
			Severity: domain.DiagnosticError,
			Code:     "tasks_unavailable_unsupported",
			Message:  "gosh tasks requires a trace-backed task index; complete task/resource/status data is unavailable because no usable artifacts were found",
			Detail:   "A complete task table requires a Nextflow trace file. No task rows were fabricated from unsupported or missing artifacts.",
		})
		return renderUnsupportedCommand(writer, diagnostics, format, fmt.Errorf("tasks: no usable artifacts in %s", runDir.Path))

	default:
		return fmt.Errorf("tasks: unsupported artifact mode %q", artifacts.Mode)
	}
}

func RunInspect(ctx context.Context, options InspectOptions, writer io.Writer) error {
	const (
		inspectSnippetMaxBytes int64 = 4096
		inspectSnippetMaxLines       = 80
	)

	command, err := loadCommandContextWithResultsDir(ctx, options.Global)
	if err != nil {
		return err
	}
	runDir := command.RunDir
	artifacts := command.Artifacts
	format := command.Format

	renderInspect := func(view domain.InspectView) error {
		view.Format = format
		switch format {
		case domain.OutputFormatJSON:
			return render.RenderInspectJSON(writer, view)
		case domain.OutputFormatHuman:
			return render.RenderInspectHuman(writer, view)
		default:
			return fmt.Errorf("unsupported output format %q", format)
		}
	}

	switch artifacts.Mode {
	case domain.IndexModeTraceBacked:
		store, _, err := idx.EnsureFreshIndex(ctx, runDir, artifacts)
		if err != nil {
			return err
		}
		defer store.Close()

		taskRows, err := idx.QueryTasks(ctx, store, domain.TaskQuery{})
		if err != nil {
			return err
		}

		resolution, err := tasks.ResolveSelector(options.Selector, taskRows)
		if err != nil {
			return err
		}

		view := domain.InspectView{Resolution: resolution, Format: format}
		if resolution.Kind != domain.SelectorResolutionExact {
			return renderInspect(view)
		}

		inventory := domain.CommandFileInventory{}
		diagnostics := []domain.Diagnostic{}
		if resolution.Task != nil && strings.TrimSpace(resolution.Task.Workdir) != "" {
			inventory, err = inspect.InventoryCommandFiles(ctx, resolution.Task.Workdir, inspect.SnippetOptions{
				MaxBytes: inspectSnippetMaxBytes,
				MaxLines: inspectSnippetMaxLines,
			})
			if err != nil {
				return fmt.Errorf("inspect command files for task %q: %w", resolution.Selector, err)
			}
		} else {
			diagnostics = append(diagnostics, domain.Diagnostic{
				Severity: domain.DiagnosticWarning,
				Code:     "inspect_workdir_unknown",
				Message:  "task workdir is unknown; command-file inventory skipped",
				Detail:   "The selected task did not include a resolvable workdir. No .command.* paths were guessed or read.",
			})
		}

		dossier, err := tasks.BuildTaskDossier(resolution, inventory)
		if err != nil {
			return err
		}
		dossier.Diagnostics = append(dossier.Diagnostics, diagnostics...)
		view.Dossier = &dossier
		return renderInspect(view)

	case domain.IndexModeLogOnly:
		if artifacts.Log == nil {
			return fmt.Errorf("inspect log-only: missing log source")
		}

		evidence, err := nflog.ParseLogOnlyTaskEvidence(ctx, runDir, *artifacts.Log)
		if err != nil {
			return err
		}

		resolution, err := tasks.ResolveLogOnlySelector(options.Selector, evidence)
		if err != nil {
			return err
		}

		view := domain.InspectView{
			EvidenceKind:      domain.InspectEvidenceLogOnly,
			LogOnlyResolution: &resolution,
			Format:            format,
		}
		if len(evidence) == 0 {
			view.Diagnostics = append(view.Diagnostics, domain.Diagnostic{
				Severity: domain.DiagnosticError,
				Code:     "log_only_no_parseable_evidence",
				Message:  "No parseable task evidence found in selected Nextflow log",
				Detail: strings.Join([]string{
					"Selected log: " + selectedLogPath(artifacts),
					"No parseable task, lifecycle, failure, or workdir evidence was found in the selected Nextflow log.",
					"The run may have failed before task evidence was emitted, or this log format is unsupported.",
					"Complete task rows and command-file workdirs require a Nextflow trace file.",
					"Hint: Run future Nextflow workflows with -with-trace to produce complete task/resource/status data.",
				}, "\n"),
			}, nextflowTraceRecommendation())
		}
		if resolution.Kind != domain.SelectorResolutionExact {
			return renderInspect(view)
		}

		inventory := domain.CommandFileInventory{}
		diagnostics := []domain.Diagnostic{
			{
				Severity: domain.DiagnosticWarning,
				Code:     "log_only_partial",
				Message:  "log-only inspect is partial",
				Detail: strings.Join([]string{
					"Selected log: " + selectedLogPath(artifacts),
					"Complete task rows require a trace file.",
					"Hint: Run future Nextflow workflows with -with-trace to produce complete task/resource/status data.",
				}, "\n"),
			},
		}
		if resolution.Evidence != nil && strings.TrimSpace(resolution.Evidence.Workdir) != "" {
			inventory, err = inspect.InventoryCommandFiles(ctx, resolution.Evidence.Workdir, inspect.SnippetOptions{
				MaxBytes: inspectSnippetMaxBytes,
				MaxLines: inspectSnippetMaxLines,
			})
			if err != nil {
				return fmt.Errorf("inspect command files for log-only selector %q: %w", resolution.Selector, err)
			}
		} else {
			diagnostics = append(diagnostics, domain.Diagnostic{
				Severity: domain.DiagnosticWarning,
				Code:     "inspect_workdir_unknown",
				Message:  "command-file inventory unavailable",
				Detail:   "The selected log-only evidence did not include a resolvable workdir.\nHint: Use the selected log error block or rerun with -with-trace.",
			})
		}

		dossier, err := tasks.BuildLogOnlyTaskDossier(resolution, inventory)
		if err != nil {
			return err
		}
		dossier.Diagnostics = append(dossier.Diagnostics, diagnostics...)
		view.LogOnlyDossier = &dossier
		return renderInspect(view)

	case domain.IndexModeUnsupported:
		diagnostics := diagnosticsForUnsupportedArtifacts(runDir, artifacts, domain.Diagnostic{
			Severity: domain.DiagnosticError,
			Code:     "inspect_unavailable_unsupported",
			Message:  "gosh inspect requires a trace-backed task index; complete task/resource/status data is unavailable because no usable artifacts were found",
			Detail:   "A task dossier requires indexed trace task rows and a known workdir. No task or command-file data were fabricated from unsupported or missing artifacts.",
		})
		return renderUnsupportedCommand(writer, diagnostics, format, fmt.Errorf("inspect: no usable artifacts in %s", runDir.Path))

	default:
		return fmt.Errorf("inspect: unsupported artifact mode %q", artifacts.Mode)
	}
}

func RunIndex(ctx context.Context, options IndexOptions, writer io.Writer) error {
	command, err := loadCommandContextWithResultsDir(ctx, options.Global)
	if err != nil {
		return err
	}
	runDir := command.RunDir
	artifacts := command.Artifacts
	format := command.Format

	renderIndex := func(diagnostics domain.IndexDiagnostics) error {
		view := domain.IndexView{Diagnostics: diagnostics, Format: format}
		switch format {
		case domain.OutputFormatJSON:
			return render.RenderIndexJSON(writer, view)
		case domain.OutputFormatHuman:
			return render.RenderIndexHuman(writer, view)
		default:
			return fmt.Errorf("unsupported output format %q", format)
		}
	}

	if options.Refresh {
		store, err := idx.OpenStore(ctx, runDir)
		if err != nil {
			return err
		}
		defer store.Close()

		var metadata domain.IndexMetadata
		if artifacts.Mode == domain.IndexModeTraceBacked {
			metadata, err = idx.RebuildTraceIndex(ctx, store, runDir, artifacts)
		} else {
			metadata, err = idx.RefreshMetadata(ctx, store, runDir, artifacts)
		}
		if err != nil {
			return err
		}

		return renderIndex(domain.IndexDiagnostics{
			RunDir:    runDir,
			Artifacts: artifacts,
			Metadata:  &metadata,
		})
	}

	diagnostics, err := idx.IndexDiagnostics(ctx, runDir, artifacts)
	if err != nil {
		return err
	}
	return renderIndex(diagnostics)
}
