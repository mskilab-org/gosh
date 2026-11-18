import csv
import datetime
import json
import os
import re
import sys
from shlex import quote
from os import makedirs, path, getenv
from shutil import rmtree
from sys import exit
from subprocess import run, CalledProcessError, TimeoutExpired
from tempfile import TemporaryDirectory
import click
from ..core.module_loader import get_environment_defaults
from ..utils.datasets import convert_to_datasets_path
from ..settings import ALL_TOOLS, ALL_TOOLS_STR, GENOME_MAP as genome_map
from ..core.presets import PRESET_NAMES, build_params, resolve_preset
from ..core.samplesheet import check_if_tumor_only


@click.group(name="run")
def run_cli():
    """Run pipeline commands"""
    pass


@run_cli.command(name="params")
@click.option("--preset", type=click.Choice(PRESET_NAMES), default="default")
@click.option("--samplesheet", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--outdir", required=True, type=click.Path())
@click.option("--reference", default="hg19", show_default=True)
@click.option("--output", required=True, type=click.Path(dir_okay=False))
def params_command(preset, samplesheet, outdir, reference, output):
    """Write a new parameter snapshot without prompts, credentials, or Slurm."""
    try:
        params = build_params(preset, samplesheet, outdir, reference)
        makedirs(path.dirname(path.abspath(output)), exist_ok=True)
        with open(output, "x") as stream:
            json.dump(params, stream, indent=4, allow_nan=False)
            stream.write("\n")
    except (OSError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(output)


def _ensure_oncokb_secret(nextflow, module_command, api_key):
    """Provision or check the named secret without exposing its value."""
    args = [nextflow, "secrets"]
    args += ["set", "ONCOKB_API_KEY", api_key] if api_key else ["list"]
    if module_command:
        # The key is only a positional argument, never part of shell source.
        args = ["/bin/bash", "-c", module_command + 'exec "$@"', "gosh-secrets", *args]
    error = (
        "Unable to provision Nextflow secret ONCOKB_API_KEY"
        if api_key
        else "Nextflow secret ONCOKB_API_KEY is required"
    )
    try:
        result = run(args, check=True, capture_output=True, text=True, timeout=60)
    except (OSError, CalledProcessError, TimeoutExpired):
        raise click.ClickException(error) from None
    if not api_key and "ONCOKB_API_KEY" not in result.stdout.splitlines():
        raise click.ClickException(error)


# ---------------------------------------------------------------------------
# gosh run trace
# ---------------------------------------------------------------------------

@run_cli.group(name="trace")
def trace_cli():
    """Parse and merge Nextflow execution traces."""
    pass


@trace_cli.command(name="merge")
@click.argument(
    "pipeline_info_dir",
    default="./results/pipeline_info",
    type=click.Path(file_okay=False),
)
@click.option(
    "-w",
    "--work-dir",
    default=None,
    type=click.Path(file_okay=False),
    help=(
        "Path to the Nextflow work/ directory. "
        "Defaults to ../../work relative to PIPELINE_INFO_DIR."
    ),
)
@click.option(
    "-o",
    "--output",
    default=None,
    type=click.Path(),
    help=(
        "Output TSV file path. "
        "Defaults to merged_trace_<YYYY-MM-DD_HH-MM-SS>.tsv inside PIPELINE_INFO_DIR. "
        "Pass '-' to write to stdout."
    ),
)
@click.option(
    "--no-work-dir",
    is_flag=True,
    default=False,
    help="Omit the work_dir column from the output.",
)
@click.option(
    "--no-trace-file",
    is_flag=True,
    default=False,
    help="Omit the trace_file provenance column from the output.",
)
def trace_merge(pipeline_info_dir, work_dir, output, no_work_dir, no_trace_file):
    """Merge all execution_trace_*.txt files in PIPELINE_INFO_DIR into one TSV.

    Each task row gains a resolved work_dir column (absolute path on disk)
    built from the short hash in the trace (e.g. 02/5852bd → work/02/<full>/),
    plus a trace_file column recording which source file the row came from.

    Rows are sorted by submit timestamp then task_id.  When the same task
    appears in multiple trace files (e.g. across resume runs) all copies are
    kept; use --hash to deduplicate externally if needed.

    \b
    Examples
    --------
    # Write a timestamped TSV into pipeline_info_dir (default)
    gosh run trace merge ./results/pipeline_info

    # Write to a specific file
    gosh run trace merge ./results/pipeline_info -o merged_trace.tsv

    # Print to stdout
    gosh run trace merge ./results/pipeline_info -o -
    """
    from .parse_trace import parse_traces

    try:
        records = parse_traces(pipeline_info_dir, work_dir=work_dir)
    except FileNotFoundError as exc:
        click.secho(f"Error: {exc}", fg="red", err=True)
        sys.exit(1)

    if not records:
        click.secho("No records found.", fg="yellow", err=True)
        sys.exit(0)

    # Determine column order: original trace columns first, then extras
    base_cols = [
        c for c in records[0].keys()
        if c not in ("work_dir", "trace_file")
    ]
    extra_cols = []
    if not no_work_dir:
        extra_cols.append("work_dir")
    if not no_trace_file:
        extra_cols.append("trace_file")
    fieldnames = base_cols + extra_cols

    # Strip columns not requested
    if no_work_dir:
        for r in records:
            r.pop("work_dir", None)
    if no_trace_file:
        for r in records:
            r.pop("trace_file", None)

    def _write(fh):
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(records)

    if output == "-":
        _write(sys.stdout)
    else:
        # Build a timestamped filename when no explicit path is given
        if output is None:
            ts = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            output = os.path.join(
                os.path.abspath(pipeline_info_dir),
                f"merged_trace_{ts}.tsv",
            )
        with open(output, "w", newline="", encoding="utf-8") as fh:
            _write(fh)
        n_resolved = sum(1 for r in records if r.get("work_dir"))
        n_total = len(records)
        click.secho(
            f"Wrote {n_total} records to {output} "
            f"({n_resolved} work dirs resolved).",
            fg="green",
            err=True,
        )




@run_cli.command(context_settings=dict(ignore_unknown_options=True))
@click.option("--pipeline-dir", help="Path to nf-gos pipeline repo")
@click.option("--samplesheet", help="Path to samplesheet CSV file")
@click.option(
    "-o", "--outdir", default="./results/", help="Path to pipeline outputs directory"
)
@click.option(
    "-r", "--reference", default="hg19", help="Genome reference version (hg19 or hg38)"
)
@click.option(
    "--params-file", default=None, help="Path to parameters JSON file (interactive default: ./params.json)"
)
@click.option("--non-interactive", is_flag=True, help="Require explicit existing params and samplesheet; never prompt.")
@click.option(
    "--profile",
    default="singularity",
    help="Pipeline profile(s) to use; comma-separated",
)
@click.option(
    "--resume/--no-resume", default=True, help="Resume previous run if possible"
)
@click.option(
    "--bg",
    "background",
    is_flag=True,
    help="Run Nextflow in the background and write output to a timestamped execution log.",
)
@click.option(
    "--gos-eye",
    is_flag=False,
    flag_value=8080,
    type=click.IntRange(1, 65535),
    metavar="[PORT]",
    help=(
        "Send Nextflow weblog events to gOS-eye at http://localhost:PORT/webhook. "
        "Omit PORT to use 8080; the monitor must already be running."
    ),
)
@click.option("--gos-eye-host", default="localhost", show_default=True, help="gOS-eye DNS hostname or IPv4 address.")
@click.option(
    "-p", "--processes", help="Comma-separated list of process names to rerun"
)
@click.option("-s", "--samples", help="Comma-separated list of sample IDs to rerun")
@click.option(
    "--skip-tools",
    help=f"Comma-separated list of tools to skip. Mutually exclusive with --include-tools. Available: {ALL_TOOLS_STR}",
)
@click.option(
    "--include-tools",
    help=f"Comma-separated list of tools to include. Mutually exclusive with --skip-tools. Available: {ALL_TOOLS_STR}",
)
@click.option(
    "--preset",
    default="default",
    type=click.Choice(PRESET_NAMES),
    help="Apply a preset to this run; default preserves existing parameter snapshots.",
)
@click.option(
    "--oncokb-api-key",
    help="Provision OncoKB credentials; otherwise reuse the named Nextflow secret.",
)
@click.option(
    "-c",
    "--config",
    "--configs",
    "config",
    help="Additional config files",
    multiple=True,
)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def pipeline(
    pipeline_dir,
    samplesheet,
    outdir,
    reference,
    params_file,
    non_interactive,
    profile,
    resume,
    background,
    gos_eye,
    gos_eye_host,
    processes,
    samples,
    skip_tools,
    include_tools,
    preset,
    oncokb_api_key,
    config,
    extra_args,
):
    # Validate mutually exclusive options
    preview = "-preview" in extra_args
    host_label = r"[A-Za-z0-9](?:[A-Za-z0-9-]{0,61}[A-Za-z0-9])?"
    if len(gos_eye_host) > 253 or not re.fullmatch(host_label + r"(?:\." + host_label + r")*", gos_eye_host):
        raise click.UsageError("--gos-eye-host requires a DNS hostname or IPv4 address")
    if gos_eye_host != "localhost" and gos_eye is None:
        raise click.UsageError("--gos-eye-host requires --gos-eye")
    native_options = {arg.split("=", 1)[0] for arg in extra_args}
    if gos_eye is not None and "-with-weblog" in native_options:
        raise click.UsageError("--gos-eye conflicts with native -with-weblog")
    if preview:
        if (background or processes is not None or samples is not None
                or preset != "default" or skip_tools is not None or include_tools is not None
                or gos_eye is not None or gos_eye_host != "localhost" or oncokb_api_key is not None):
            raise click.UsageError("Preview cannot use execution, preset, monitoring, or credential options")
        forbidden = {"-resume", "-bg", "-with-weblog", "-with-report", "-with-trace",
                     "-with-timeline", "-with-notification", "-N"}
        if native_options & forbidden:
            raise click.UsageError("Preview cannot use native execution/report/notification options")
        if not path.isfile(params_file or "./params.json"):
            raise click.UsageError("Preview requires an existing params file; it never runs the wizard")
    if skip_tools is not None and include_tools is not None:
        raise click.UsageError(
            "Options --skip-tools and --include-tools are mutually exclusive."
        )
    if non_interactive:
        if not params_file or not samplesheet:
            raise click.UsageError("--non-interactive requires explicit --params-file and --samplesheet")
        if processes is not None or samples is not None:
            raise click.UsageError("--non-interactive cannot use --processes or --samples rerun options")
        if not path.isfile(params_file):
            raise click.ClickException(f"Parameters file '{params_file}' not found")
        if not path.isfile(samplesheet):
            raise click.ClickException(f"Samplesheet '{samplesheet}' not found")
    params_file = params_file or "./params.json"

    from ..core.nextflow import NextflowRunner
    from ..core.params_wizard import create_params_file
    from ..core.nextflow_log import (
        get_entries_with_process_names,
        get_entries_with_sample_names,
    )
    from ..core.module_loader import load_required_modules

    is_samplesheet_provided = not samplesheet is None
    default_samplesheet = "./samplesheet.csv"
    if not is_samplesheet_provided:
        samplesheet = default_samplesheet

    # Retrieve environment defaults - already imported above
    env_defaults = get_environment_defaults()

    # Set pipeline_dir if not specified
    if not pipeline_dir:
        pipeline_dir = env_defaults.get("pipeline-dir", pipeline_dir)

    # Host defaults augment explicit profiles without duplicating them.
    profiles = profile.split(",") + env_defaults.get("profile", "").split(",")
    profile = ",".join(dict.fromkeys(p.strip() for p in profiles if p.strip()))


    # Check if params_file is provided and exists
    is_params_file_found = path.isfile(params_file)
    default_params = "./params.json"
    is_default_params_file_found = path.isfile(default_params)
    if not is_params_file_found:
        print(f"Parameters file '{params_file}' not found.")
        # Check if default params.json exists
        default_params = "./params.json"
        if not is_default_params_file_found:
            print("No params.json file found. Launching wizard to create one...")
            # Call the wizard to create params.json with the preset supplied in the run command
            create_params_file(
                preset=preset, samplesheet=samplesheet, outdir=outdir, genome=reference
            )
            params_file = default_params
        else:
            params_file = default_params

    # Determine the skip_tools value based on flags and presets
    skip_tools_value = None
    if include_tools is not None:
        # Calculate skip_tools based on include_tools
        included_set = {t.strip() for t in include_tools.split(",")} if include_tools else set()
        # Validate included tools
        invalid_tools = included_set - set(ALL_TOOLS)
        if invalid_tools:
            raise click.BadParameter(
                f"Invalid tools specified in --include-tools: {', '.join(invalid_tools)}. Available: {ALL_TOOLS_STR}"
            )
        skipped_set = set(ALL_TOOLS) - included_set
        skip_tools_value = ",".join(sorted(list(skipped_set)))
        print(f"Calculating skip_tools based on --include-tools: {skip_tools_value}")
    elif skip_tools is not None:
        # Validate skipped tools
        skipped_set = {t.strip() for t in skip_tools.split(",")} if skip_tools else set()
        invalid_tools = skipped_set - set(ALL_TOOLS)
        if invalid_tools:
            raise click.BadParameter(
                f"Invalid tools specified in --skip-tools: {', '.join(invalid_tools)}. Available: {ALL_TOOLS_STR}"
            )
        skip_tools_value = skip_tools  # Use provided skip_tools directly

    # Read and update params.json with CLI flag values
    try:
        with open(params_file, "r") as pf:
            params_data = json.load(pf)
    except (OSError, ValueError) as exc:
        raise click.ClickException(f"Cannot read parameters file: {params_file}") from exc
    if not isinstance(params_data, dict):
        raise click.ClickException("Parameters file must be a JSON object")

    # Overwrite keys from CLI flags

    if not is_samplesheet_provided:
        samplesheet = params_data.get("input", default_samplesheet)

    if not os.path.exists(samplesheet):
        raise ValueError(
            f'Path "{samplesheet}" does not exist! Please provide a valid samplesheet path'
        )
    if preview or non_interactive or preset == "tumor-only":
        try:
            tumor_only = check_if_tumor_only(samplesheet)
        except (OSError, ValueError) as exc:
            raise click.ClickException(str(exc)) from exc
        if (preview or non_interactive) and (
            type(params_data.get("tumor_only")) is not bool
            or params_data["tumor_only"] != tumor_only
        ):
            raise click.ClickException("params.tumor_only is inconsistent with the samplesheet")
        if preset == "tumor-only" and not tumor_only:
            raise click.ClickException("The tumor-only preset requires every sample to have status 1")
    # A prepared run deliberately omits --preset, keeping its immutable snapshot.
    if preset != "default":
        try:
            params_data.update(resolve_preset(preset))
        except (OSError, ValueError) as exc:
            raise click.ClickException(str(exc)) from exc
    params_data["input"] = samplesheet
    params_data["outdir"] = outdir

    # Convert the provided reference to the corresponding genome value
    genome_retrieved = params_data.get("genome", None)

    is_genome_still_not_set = genome_retrieved is None
    errmsg = f"reference/genome must be one of hg19, GRCh37, GATK.GRCh37, hg38, GRCh38, or GATK.GRCh38 in {params_file}"

    if is_genome_still_not_set:
        raise ValueError(errmsg)

    is_genome_invalid = genome_map.get(genome_retrieved.lower(), None) is None
    if is_genome_invalid:
        raise ValueError(errmsg)

    # Also update skip_tools if a skip_tools value was computed
    if skip_tools_value is not None:
        params_data["skip_tools"] = skip_tools_value

    if not preview:
        # Normal execution retains its params rewrite and reference directory.
        makedirs(reference, exist_ok=True)
        with open(params_file, "w") as pf:
            json.dump(params_data, pf, indent=4)


    # Print all parameters
    print("Running gOS with the following parameters:")
    print(f"Pipeline directory: {pipeline_dir}")
    print(f"Samplesheet: {samplesheet}")
    print(f"Parameters file: {params_file}")
    print(f"Profile: {profile}")
    print(f"Resume: {resume}")

    # Initialize a set to collect work directories to delete
    workdirs_to_delete = {}

    # Check if processes or samples are specified
    if processes or samples:
        # Parse the comma-separated lists into Python lists
        processes_list = [p.strip() for p in processes.split(",")] if processes else []
        samples_list = [s.strip() for s in samples.split(",")] if samples else []

        workdirs_processes = {}
        if processes_list:
            print("Retrieving work directories for specified processes...")
            entries = get_entries_with_process_names(processes_list)
            for entry in entries:
                workdir = entry["workdir"]
                name = entry["name"]
                if workdir not in workdirs_processes:
                    workdirs_processes[workdir] = set()
                workdirs_processes[workdir].add(name)

        workdirs_samples = {}
        if samples_list:
            print("Retrieving work directories for specified samples...")
            entries = get_entries_with_sample_names(samples_list)
            for entry in entries:
                workdir = entry["workdir"]
                name = entry["name"]
                if workdir not in workdirs_samples:
                    workdirs_samples[workdir] = set()
                workdirs_samples[workdir].add(name)

        if processes_list and samples_list:
            workdirs_set = set(workdirs_processes.keys()) & set(workdirs_samples.keys())
        elif processes_list:
            workdirs_set = set(workdirs_processes.keys())
        elif samples_list:
            workdirs_set = set(workdirs_samples.keys())
        else:
            workdirs_set = set()

        for workdir in workdirs_set:
            names = set()
            if workdir in workdirs_processes:
                names.update(workdirs_processes[workdir])
            if workdir in workdirs_samples:
                names.update(workdirs_samples[workdir])
            workdirs_to_delete[workdir] = names

        if workdirs_to_delete:
            print("The following work directories will be deleted:")
            print("Names, Work Directory")
            for workdir, names in sorted(workdirs_to_delete.items()):
                names_str = ", ".join(sorted(names))
                print(f"{names_str}, {workdir}")

            confirm = (
                input("Do you want to delete these directories? ([y]es/[n]o): ")
                .strip()
                .lower()
            )
            if confirm == "yes" or confirm == "y":
                confirm_again = (
                    input(
                        "Are you sure? This will cause the pipeline to rerun from these steps. ([y]es/[n]o): "
                    )
                    .strip()
                    .lower()
                )
                if confirm_again == "yes" or confirm_again == "y":
                    for workdir in workdirs_to_delete:
                        if workdir and path.exists(workdir):
                            print(f"Deleting work directory: {workdir}")
                            rmtree(workdir)
                    print("Directories deleted.")
                else:
                    print("Run cancelled.")
                    exit(0)
            else:
                print("Run cancelled.")
                exit(0)
        else:
            print("No matching work directories found to delete.")
            exit(0)

    load_modules_command = load_required_modules(env_defaults)
    if load_modules_command:
        load_modules_command = "source /etc/profile.d/modules.sh && " + load_modules_command

    runner = NextflowRunner()
    if not preview:
        _ensure_oncokb_secret(
            runner.cmd, load_modules_command, oncokb_api_key or getenv("ONCOKB_API_KEY")
        )


    configs_lst = []
    for conf in config:
        is_conf_absolute_path = conf.startswith("/")
        is_conf_user_path = conf.startswith("~")
        is_conf_path_resolved = is_conf_absolute_path or is_conf_user_path
        conf = os.path.expanduser(conf) if is_conf_user_path else conf
        out_conf = (
            conf if is_conf_path_resolved else pipeline_dir + "/conf/settings/" + conf
        )
        if not os.path.exists(out_conf):
            print(f"{out_conf} does not exist!")
            next
        configs_lst.append(out_conf)

    if preview:
        configs_lst.append(path.abspath(path.join(path.dirname(__file__), "..", "core", "preview.config")))
    configs = " ".join(f"-c {quote(resolved_config)}" for resolved_config in configs_lst)
    run_timestamp = runner.get_timestamp()
    if preview:
        # Keep the snapshot private and alive until Nextflow has exited. Relative
        # input/output paths still resolve in the unchanged study working directory.
        with TemporaryDirectory(prefix=".gosh-preview-", dir=path.dirname(path.abspath(params_file))) as temp_dir:
            effective_params = path.join(temp_dir, "params.json")
            with open(effective_params, "w") as stream:
                json.dump(params_data, stream, allow_nan=False)
            command = (
                f"{load_modules_command} NXF_SYNTAX_PARSER=v1 {quote(runner.cmd)} {configs} "
                f"-log .nextflow-preview_{run_timestamp}.log run {quote(pipeline_dir)} "
                f"-params-file {quote(effective_params)} -profile {quote(profile)} "
                + " ".join(quote(arg) for arg in extra_args)
            )
            print(f"command to run: {command}")
            result = run(command, shell=True, executable="/bin/bash")
            if result.returncode:
                raise click.exceptions.Exit(result.returncode)
        return
    # nf-gos uses legacy DSL constructs; Nextflow >=26.04 defaults to the strict parser.
    nextflow_command = f"NXF_SYNTAX_PARSER=v1 {runner.cmd} "
    pid_file = None
    if background:
        pid_file = f".nextflow_{run_timestamp}.pid"
        nextflow_command = f"NXF_PID_FILE={pid_file} {nextflow_command}"

    command = (
        f"{load_modules_command} "
        f"{nextflow_command}"
        f"{configs} "
        f"-log .nextflow_{run_timestamp}.log "
        f"run {quote(pipeline_dir)} "
        f"-params-file {quote(params_file)} "
        f"-profile {quote(profile)} "
        f"-with-report report_{run_timestamp}.html "
        f"-with-trace"
    )

    if resume:
        command += " -resume"

    if gos_eye is not None:
        command += " -with-weblog " + quote(f"http://{gos_eye_host}:{gos_eye}/webhook")

    if not extra_args is None and not len(extra_args) == 0:
        command += f' {" ".join(quote(arg) for arg in extra_args)} '


    if background:
        execution_log = datetime.datetime.now().strftime(
            "execution_%Y-%m-%d.%H.%M.%S.log"
        )
        command += f" -bg > {execution_log} 2>&1"
    print(f"command to run: {command}")
    runner.run(command)

    if background:
        runner.finalize_after_background_run(pid_file, execution_log)
        return

    runner.run("gosh run post")
    runner.run("gosh protect")


@run_cli.command()
@click.option(
    "-p",
    "--pipeline-output-dir",
    default="./results/",
    required=False,
    type=click.Path(exists=True),
    help="Directory containing pipeline outputs",
)
@click.option(
    "--prefer-outputs/--no-prefer-outputs",
    default=True,
    show_default=True,
    help="Prefer matching results-directory paths; retain samplesheet paths when no output matches.",
)
def post(pipeline_output_dir, prefer_outputs):
    from ..core.module_loader import (
        load_required_modules,
    )  # Keep this specific import if needed elsewhere

    env_defaults = get_environment_defaults()
    load_modules_command = load_required_modules(env_defaults)
    preference = "--prefer-outputs" if prefer_outputs else "--no-prefer-outputs"
    output_dir = quote(pipeline_output_dir)
    command = (
        f"{load_modules_command} "
        f"{{ gosh run outputs --force -p {output_dir} {preference}; "
        f"gosh run samplesheet --force -p {output_dir} {preference}; "
        f"gosh debug log > debug.log; }}"
    )
    try:
        run(command, shell=True, check=True, executable="/bin/bash")
    except CalledProcessError as e:
        raise RuntimeError(f"gosh post execution failed: {e}")


from ..core.outputs import OUTPUT_KEYS


@run_cli.command()
@click.option(
    "-p",
    "--pipeline-output-dir",
    default="./results/",
    required=False,
    type=click.Path(exists=True),
    help="Directory containing pipeline outputs",
)
@click.option(
    "-s",
    "--samplesheet",
    default="./samplesheet.csv",
    required=False,
    type=click.Path(exists=True),
    help="Path to the samplesheet CSV file",
)
@click.option(
    "--old",
    is_flag=True,
    default=False,
    help="Whether to use the old outputs mapping (default: False)",
)
@click.option(
    "--prefer-outputs/--no-prefer-outputs",
    default=True,
    show_default=True,
    help="Prefer matching results-directory paths; retain samplesheet paths when no output matches.",
)
@click.option(
    "-o",
    "--output",
    default="outputs.csv",
    type=click.Path(),
    help="CSV file to  save outputs default: (./outputs.csv)",
)
@click.option(
    "-f",
    "--force",
    default=False,
    is_flag=True,
    help="Overwrite output file specified by -o/--output flag (default: False)",
)
@click.option(
    "-c",
    "--include-columns",
    help="Comma-separated list of columns to include. Available: {}".format(
        ",".join(OUTPUT_KEYS)
    ),
)
@click.option(
    "-C",
    "--exclude-columns",
    help="Comma-separated list of columns to exclude. Available: see --include-columns",
)
def outputs(
    pipeline_output_dir,
    samplesheet,
    old,
    prefer_outputs,
    output,
    include_columns,
    exclude_columns,
    force,
):
    """
    Generate an outputs.csv file suitable for skilifting.
    Reads pipeline outputs and samplesheet to create a CSV mapping patient data to output file paths.
    """
    from ..core.outputs import (
        Outputs,
        OUTPUT_KEYS,
    )  # Import OUTPUT_KEYS for help text formatting

    if include_columns and exclude_columns:
        raise click.UsageError(
            "Options -c/--include-columns and -C/--exclude-columns are mutually exclusive."
        )

    include_list = include_columns.split(",") if include_columns else None
    exclude_list = exclude_columns.split(",") if exclude_columns else None

    # Validate column names if provided
    available_columns = set(OUTPUT_KEYS)
    if include_list:
        invalid_cols = [col for col in include_list if col not in available_columns]
        if invalid_cols:
            raise click.BadParameter(
                f"Invalid columns specified in --include-columns: {', '.join(invalid_cols)}. Available: {', '.join(OUTPUT_KEYS)}"
            )
    if exclude_list:
        invalid_cols = [col for col in exclude_list if col not in available_columns]
        if invalid_cols:
            raise click.BadParameter(
                f"Invalid columns specified in --exclude-columns: {', '.join(invalid_cols)}. Available: {', '.join(OUTPUT_KEYS)}"
            )

    is_pipeline_output_valid = os.path.isdir(pipeline_output_dir)

    errmsg = []
    raise_error = False
    if not is_pipeline_output_valid:
        raise_error = True
        errmsg.append(
            "-p/--pipeline-output-dir flag must be set to an existing outputs directory, default ./results/"
        )

    is_samplesheet_valid = os.path.isfile(samplesheet)

    if not is_samplesheet_valid:
        raise_error = True
        errmsg.append(
            "-s/--samplesheet flag must be set to an existing samplesheet, default ./samplesheet.csv"
        )

    if raise_error:
        raise Exception("; ".join(errmsg))

    is_output_present_already = isinstance(output, str) and os.path.isfile(output)
    do_overwrite = is_output_present_already and force
    do_error_out = is_output_present_already and not force
    if do_overwrite:
        print(
            f"Warning: --force flag is set, {output} file exists and will be overwritten!"
        )
    elif do_error_out:
        raise Exception(
            f'{output} file exists, but --force flag is not set! If you wish to overwrite the output file, explicitly set the flag "--force"'
        )

    outputs_obj = Outputs(pipeline_output_dir, samplesheet, old, prefer_outputs)
    outputs_obj.emit_output_csv(
        output, include_columns=include_list, exclude_columns=exclude_list
    )

    if output:
        click.echo(f"Outputs CSV generated at: {output}")


from ..core.outputs import SAMPLESHEET_FIELDNAMES


@run_cli.command()
@click.option(
    "-p",
    "--pipeline-output-dir",
    default="./results/",
    required=False,
    type=click.Path(exists=True),
    help="Directory containing pipeline outputs",
)
@click.option(
    "-s",
    "--samplesheet",
    default="./samplesheet.csv",
    required=False,
    type=click.Path(exists=True),
    help="Path to the samplesheet CSV file",
)
@click.option(
    "--old",
    is_flag=True,
    default=False,
    help="Whether to use the old outputs mapping (default: False)",
)
@click.option(
    "--prefer-outputs/--no-prefer-outputs",
    default=True,
    show_default=True,
    help="Prefer matching results-directory paths; retain samplesheet paths when no output matches.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    default="samplesheet_from_outputs.csv",
    help="CSV file to save samplesheet (default: stdout)",
)
@click.option(
    "-f",
    "--force",
    type=bool,
    default=False,
    is_flag=True,
    help="Overwrite output file specified by -o/--output (default: False)",
)
@click.option(
    "-c",
    "--include-columns",
    help="Comma-separated list of columns to include. Available: {}".format(
        ",".join(SAMPLESHEET_FIELDNAMES)
    ),
)
@click.option(
    "-C",
    "--exclude-columns",
    help="Comma-separated list of columns to exclude. Available: see --include-columns",
)
def samplesheet(
    pipeline_output_dir,
    samplesheet,
    old,
    prefer_outputs,
    output,
    include_columns,
    exclude_columns,
    force,
):
    """
    Generate a samplesheet.csv file suitable for a pipeline run.
    Reads pipeline outputs and the original samplesheet to create a new samplesheet
    with paths to generated files filled in.
    """
    from ..core.outputs import Outputs

    if include_columns and exclude_columns:
        raise click.UsageError(
            "Options -c/--include-columns and -C/--exclude-columns are mutually exclusive."
        )

    include_list = include_columns.split(",") if include_columns else None
    exclude_list = exclude_columns.split(",") if exclude_columns else None

    # Validate column names if provided
    available_columns = set(SAMPLESHEET_FIELDNAMES)
    if include_list:
        invalid_cols = [col for col in include_list if col not in available_columns]
        if invalid_cols:
            raise click.BadParameter(
                f"Invalid columns specified in --include-columns: {', '.join(invalid_cols)}. Available: {', '.join(SAMPLESHEET_FIELDNAMES)}"
            )
    if exclude_list:
        invalid_cols = [col for col in exclude_list if col not in available_columns]
        if invalid_cols:
            raise click.BadParameter(
                f"Invalid columns specified in --exclude-columns: {', '.join(invalid_cols)}. Available: {', '.join(SAMPLESHEET_FIELDNAMES)}"
            )

    is_pipeline_output_valid = os.path.isdir(pipeline_output_dir)

    errmsg = []
    raise_error = False
    if not is_pipeline_output_valid:
        raise_error = True
        errmsg.append(
            "-p/--pipeline-output-dir flag must be set to an existing outputs directory, default ./results/"
        )

    is_samplesheet_valid = os.path.isfile(samplesheet)

    if not is_samplesheet_valid:
        raise_error = True
        errmsg.append(
            "-s/--samplesheet flag must be set to an existing samplesheet, default ./samplesheet.csv"
        )

    if raise_error:
        raise Exception("; ".join(errmsg))

    is_output_present_already = isinstance(output, str) and os.path.isfile(output)
    do_overwrite = is_output_present_already and force
    output_path = output
    if do_overwrite:
        print(
            f"Warning: --force flag is set, {output} file exists and will be overwritten!"
        )
    elif is_output_present_already:
        timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        original_filename = output
        filename_base, file_extension = os.path.splitext(original_filename)
        new_filename = f"{filename_base}_{timestamp}{file_extension}"
        output_path = new_filename
        print(
            f"Warning: --force flag is not set, {output} file exists, so output will be written to {output_path}!"
        )

    outputs_obj = Outputs(pipeline_output_dir, samplesheet, old, prefer_outputs)
    outputs_obj.emit_samplesheet_csv(
        output_path, include_columns=include_list, exclude_columns=exclude_list
    )

    if output:
        click.echo(f"Samplesheet CSV generated at: {output_path}")


@run_cli.command()
@click.option(
    "-p",
    "--pipeline-output-dir",
    default="./results/",
    type=click.Path(exists=True),
    help="Directory containing pipeline outputs",
)
@click.option(
    "-s",
    "--samplesheet",
    default="./samplesheet.csv",
    type=click.Path(exists=True),
    help="Path to the samplesheet CSV file used to generate the outputs.csv",
)
@click.option(
    "--old",
    is_flag=True,
    default=False,
    help="Whether the outputs.csv was generated using the old process_name/patient_id nf-gos outputs mapping (default: False)",
)
@click.option(
    "--output-csv",
    default="./outputs.csv",
    type=click.Path(exists=True),
    help="Path to the outputs.csv file generated by 'gosh run outputs'.",
)
@click.option(
    "-t",
    "--cohort-type",
    type=click.Choice(["paired", "tumor_only", "heme"]),
    default="paired",
    help="Type of the cohort.",
)
@click.option(
    "-o",
    "--gos_dir",
    type=click.Path(),
    required=True,
    help="Path to where the skilifted outputs should be deposited.",
)
@click.option(
    "-l",
    "--skilift-repo",
    type=click.Path(),
    default="~/git/skilift",
    help="Path to the skilift repo (default: ~/git/skilift)",
)
@click.option("-c", "--cores", type=int, default=1, help="Number of cores to use.")
def skilift(
    cohort_type,
    gos_dir,
    skilift_repo,
    cores,
    old,
    output_csv=None,
    pipeline_output_dir=None,
    samplesheet=None,
):
    """Lift raw data into gOS compatible formats using Skilift."""
    from os import path
    import subprocess

    # Need these imports here for the help text formatting in the options above
    from ..core.outputs import OUTPUT_KEYS, SAMPLESHEET_FIELDNAMES

    # Expand the skilift_repo path
    skilift_repo = path.expanduser(skilift_repo)
    if not path.isdir(skilift_repo):
        if click.confirm(
            f"Skilift repository not found at {skilift_repo}. Would you like to clone it?",
            default=True,
        ):
            clone_url = "https://github.com/mskilab-org/skilift.git"
            click.echo(f"Cloning skilift from {clone_url} to {skilift_repo} ...")
            try:
                subprocess.run(
                    ["git", "clone", clone_url, skilift_repo],
                    check=True,
                    executable="/bin/bash",
                )
            except Exception as err:
                click.echo(f"Error cloning skilift repository: {err}")
                return
        else:
            click.echo("Skilift repository is required. Exiting.")
            return

    is_samplesheet_present = isinstance(samplesheet, str) and os.path.isfile(
        samplesheet
    )
    is_output_csv_present = isinstance(output_csv, str) and os.path.isfile(output_csv)
    is_pipeline_output_dir_present = isinstance(
        pipeline_output_dir, str
    ) and os.path.isdir(pipeline_output_dir)

    if (not pipeline_output_dir and not is_pipeline_output_dir_present) and (
        not output_csv and not is_output_csv_present
    ):
        click.echo(
            "Error: Either --pipeline-output-dir or --output-csv must be provided."
        )
        return

    is_output_csv_absent = not output_csv and not is_output_csv_present
    if is_output_csv_absent:
        if not is_samplesheet_present:
            click.echo(
                f"Error: The samplesheet '{samplesheet}' does not exist. Please provide a samplesheet with the -s option."
            )
            return
        if not path.isdir(pipeline_output_dir):
            click.echo(
                f"Error: The directory '{pipeline_output_dir}' does not exist. Please provide a valid pipeline output directory with the -p option."
            )
            return

        # Generate the outputs.csv file if not provided
        output_csv = "./outputs.csv"
        click.echo(f"Generating temporary outputs CSV at: {output_csv}")
        from ..core.outputs import Outputs

        try:
            outputs_obj = Outputs(pipeline_output_dir, samplesheet, old)
            outputs_obj.emit_output_csv(output_csv)
            click.echo(f"Outputs CSV generated successfully.")
        except Exception as e:
            click.echo(f"Error generating outputs CSV: {e}")
            return

    if not path.isfile(output_csv):
        click.echo(
            f"Error: The outputs CSV file '{output_csv}' does not exist or could not be generated."
        )
        return

    makedirs(path.expanduser(gos_dir), exist_ok=True)

    r_code = f"""
    devtools::load_all("{path.expanduser(skilift_repo)}")
    cohort <- Cohort$new("{output_csv}", cohort_type="{cohort_type}")
    saveRDS(cohort, "{gos_dir}/cohort.rds")
    cohort_mod <- lift_all(cohort, output_data_dir="{gos_dir}", cores={cores})
    message("Saving over cohort")
    saveRDS(cohort_mod, "{gos_dir}/cohort.rds")
    """

    try:
        run(["Rscript", "-e", r_code], check=True)
    except CalledProcessError as e:
        click.echo(f"Error executing R script: {e}")


@run_cli.command()
@click.option(
    "-d",
    "--datasets-json",
    type=click.Path(),
    help="Path to the datasets.json file. Defaults to value from environment config.",
)
@click.option(
    "-g",
    "--gos-data-dir",
    required=True,
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
    help="Path to the directory containing skilifted gOS data and datafiles.json.",
)
@click.option(
    "-t", "--title", required=True, type=str, help="Title for the new dataset cohort."
)
def datasets(datasets_json, gos_data_dir, title):
    """
    Add a new dataset entry to the datasets.json file for gOS webapp consumption.
    """
    env_defaults = get_environment_defaults()
    default_datasets_json = env_defaults.get("datasets_json")

    if not datasets_json:
        if not default_datasets_json:
            click.echo(
                "Error: --datasets-json path must be provided or configured in environment defaults.",
                err=True,
            )
            exit(1)
        datasets_json = default_datasets_json
        click.echo(
            f"Using default datasets.json path from environment: {datasets_json}"
        )
    else:
        # Ensure the directory for the datasets.json exists if a path is provided
        datasets_dir = os.path.dirname(os.path.abspath(datasets_json))
        os.makedirs(datasets_dir, exist_ok=True)

    # Validate gos_data_dir contains datafiles.json
    datafiles_path = os.path.join(gos_data_dir, "datafiles.json")
    if not os.path.isfile(datafiles_path):
        click.echo(
            f"Error: 'datafiles.json' not found in the specified gOS data directory: {gos_data_dir}",
            err=True,
        )
        exit(1)

    # Generate ID from title
    dataset_id = title.lower().replace(" ", "_")

    # Convert gos_data_dir to the required URL format
    try:
        data_path_url = convert_to_datasets_path(gos_data_dir)
    except ValueError as e:
        click.echo(f"Error converting gOS data directory path: {e}", err=True)
        click.echo(
            "Please ensure the path follows the expected structure (e.g., /path/to/{lab_name}/external/{rest_of_path})",
            err=True,
        )
        exit(1)

    datafiles_path_url = data_path_url + "datafiles.json"

    # Create the new dataset entry
    new_entry = {
        "id": dataset_id,
        "title": title,
        "dataPath": data_path_url,
        "datafilesPath": datafiles_path_url,
    }

    # Read existing datasets.json or initialize if it doesn't exist
    all_datasets = []
    if os.path.exists(datasets_json):
        try:
            with open(datasets_json, "r") as f:
                all_datasets = json.load(f)
                if not isinstance(all_datasets, list):
                    click.echo(
                        f"Error: Expected {datasets_json} to contain a JSON list.",
                        err=True,
                    )
                    exit(1)
        except json.JSONDecodeError:
            click.echo(
                f"Error: Could not decode JSON from {datasets_json}. Initializing with new entry.",
                err=True,
            )
            all_datasets = []  # Reset if file is corrupt
        except Exception as e:
            click.echo(f"Error reading {datasets_json}: {e}", err=True)
            exit(1)

    # Check if an entry with the same ID or title already exists
    existing_ids = {entry.get("id") for entry in all_datasets if "id" in entry}
    existing_titles = {entry.get("title") for entry in all_datasets if "title" in entry}

    if new_entry["id"] in existing_ids:
        click.echo(
            f"Error: Dataset with ID '{new_entry['id']}' already exists in {datasets_json}.",
            err=True,
        )
        exit(1)
    if new_entry["title"] in existing_titles:
        click.echo(
            f"Warning: Dataset with title '{new_entry['title']}' already exists in {datasets_json}. Proceeding with unique ID '{new_entry['id']}'.",
            err=True,
        )

    # Append the new entry
    all_datasets.append(new_entry)

    # Write the updated list back to datasets.json
    try:
        with open(datasets_json, "w") as f:
            json.dump(all_datasets, f, indent=4)
        click.echo(
            f"Successfully added dataset '{title}' (ID: {dataset_id}) to {datasets_json}"
        )
        click.echo(f"  dataPath: {data_path_url}")
        click.echo(f"  datafilesPath: {datafiles_path_url}")
    except Exception as e:
        click.echo(f"Error writing updated data to {datasets_json}: {e}", err=True)
        exit(1)
