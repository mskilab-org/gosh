import json
import os
from os import makedirs, path, getenv
from shutil import rmtree
from sys import exit
from subprocess import run, CalledProcessError
import click
from ..core.module_loader import get_environment_defaults
from ..utils.datasets import convert_to_datasets_path

@click.group(name='run')
def run_cli():
    """Run pipeline commands"""
    pass

# Define the full list of tools available in the pipeline
ALL_TOOLS = [
    "aligner", "collect_wgs_metrics", "collect_multiple_metrics",
    "estimate_library_complexity", "msisensorpro", "gridss", "amber",
    "fragcounter", "dryclean", "cbs", "sage", "purple", "jabba",
    "non_integer_balance", "lp_phased_balance", "events", "fusions", "snpeff",
    "snv_multiplicity", "oncokb", "signatures", "hrdetect", "onenesstwoness"
]
ALL_TOOLS_STR = ",".join(ALL_TOOLS)

@run_cli.command()
@click.option('--pipeline-dir',
              help='Path to nf-gos pipeline repo')
@click.option('--samplesheet',
              default='./samplesheet.csv',
              help='Path to samplesheet CSV file')
@click.option('-o', '--outdir',
              default='./results/',
              help='Path to pipeline outputs directory')
@click.option('-r', '--reference',
              default='hg19',
              help='Genome reference version (hg19 or hg38)')
@click.option('--params-file',
              default='./params.json',
              help='Path to parameters JSON file')
@click.option('--profile',
              default='singularity',
              help='Pipeline profile(s) to use; comma-separated')
@click.option('--resume/--no-resume',
              default=True,
              help='Resume previous run if possible')
@click.option('-p', '--processes',
              help='Comma-separated list of process names to rerun')
@click.option('-s', '--samples',
              help='Comma-separated list of sample IDs to rerun')
@click.option('--skip-tools',
              help=f'Comma-separated list of tools to skip. Mutually exclusive with --include-tools. Available: {ALL_TOOLS_STR}')
@click.option('--include-tools',
              help=f'Comma-separated list of tools to include. Mutually exclusive with --skip-tools. Available: {ALL_TOOLS_STR}')
@click.option('--preset',
              default='default',
              type=click.Choice(['default', 'jabba', 'hrd', 'heme']),
              help='Preset option: "default" (all tools), "jabba", "hrd", or "heme"')
@click.option('--oncokb-api-key',
              help='OncoKB API key for accessing OncoKB annotations. Required if using OncoKB')
def pipeline(
    pipeline_dir,
    samplesheet,
    outdir,
    reference,
    params_file,
    profile,
    resume,
    processes,
    samples,
    skip_tools,
    include_tools,
    preset,
    oncokb_api_key
):
    # Validate mutually exclusive options
    if skip_tools and include_tools:
        raise click.UsageError("Options --skip-tools and --include-tools are mutually exclusive.")

    from ..core.nextflow import NextflowRunner
    from ..core.params_wizard import create_params_file
    from ..core.nextflow_log import get_entries_with_process_names, get_entries_with_sample_names
    from ..core.nextflow import NextflowRunner
    from ..core.params_wizard import create_params_file
    from ..core.nextflow_log import get_entries_with_process_names, get_entries_with_sample_names
    from ..core.module_loader import load_required_modules # Keep this specific import if needed elsewhere

    # Retrieve environment defaults - already imported above
    env_defaults = get_environment_defaults()

    # Set pipeline_dir if not specified
    if not pipeline_dir:
        pipeline_dir = env_defaults.get('pipeline-dir', pipeline_dir)

    # Add the environment-specific profile to the default profile
    if 'profile' in env_defaults:
        profile = f"{profile},{env_defaults['profile']}"

    # Create hg19 directory (required for fragcounter)
    makedirs(reference, exist_ok=True)

    # Check if params_file is provided and exists
    if not path.isfile(params_file):
        print(f"Parameters file '{params_file}' not found.")
        # Check if default params.json exists
        default_params = './params.json'
        if not path.isfile(default_params):
            print("No params.json file found. Launching wizard to create one...")
            # Call the wizard to create params.json with the preset supplied in the run command
            create_params_file(
                preset=preset,
                samplesheet=samplesheet,
                outdir=outdir,
                genome=reference
            )
            params_file = default_params
        else:
            params_file = default_params

    # Determine the skip_tools value based on flags and presets
    skip_tools_value = None
    if include_tools:
        # Calculate skip_tools based on include_tools
        included_set = set(t.strip() for t in include_tools.split(','))
        # Validate included tools
        invalid_tools = included_set - set(ALL_TOOLS)
        if invalid_tools:
             raise click.BadParameter(f"Invalid tools specified in --include-tools: {', '.join(invalid_tools)}. Available: {ALL_TOOLS_STR}")
        skipped_set = set(ALL_TOOLS) - included_set
        skip_tools_value = ",".join(sorted(list(skipped_set)))
        print(f"Calculating skip_tools based on --include-tools: {skip_tools_value}")
    elif skip_tools:
        # Validate skipped tools
        skipped_set = set(t.strip() for t in skip_tools.split(','))
        invalid_tools = skipped_set - set(ALL_TOOLS)
        if invalid_tools:
             raise click.BadParameter(f"Invalid tools specified in --skip-tools: {', '.join(invalid_tools)}. Available: {ALL_TOOLS_STR}")
        skip_tools_value = skip_tools # Use provided skip_tools directly
    elif preset != 'default':
        # Calculate skip_tools based on preset if no specific flags are given
        if preset == 'jabba':
            skip_tools_value = "sage,snpeff,snv_multiplicity,signatures,hrdetect"
        elif preset == 'hrd':
            skip_tools_value = "non_integer_balance,lp_phased_balance,events,fusions"
        elif preset == 'heme':
            skip_tools_value = "msisensorpro,hrdetect,onenesstwoness"
        print(f"Using preset '{preset}', setting skip_tools to: {skip_tools_value}")
    # else: skip_tools_value remains None (or empty string equivalent in params)

    # Read and update params.json with CLI flag values
    with open(params_file, "r") as pf:
        params_data = json.load(pf)

    # Overwrite keys from CLI flags
    params_data["input"] = samplesheet
    params_data["outdir"] = outdir

    # Convert the provided reference to the corresponding genome value
    genome_map = {"hg19": "GATK.GRCh37", "hg38": "GATK.GRCh38"}
    params_data["genome"] = genome_map.get(reference, "GATK.GRCh37")

    # Also update skip_tools if a skip_tools value was computed
    if skip_tools_value is not None:
        params_data["skip_tools"] = skip_tools_value

    with open(params_file, "w") as pf:
        json.dump(params_data, pf, indent=4)

    oncokb_api_key = oncokb_api_key or getenv('ONCOKB_API_KEY', None)
    if skip_tools_value is not None and 'oncokb' not in skip_tools_value and not oncokb_api_key:
        print("OncoKB API key is required for accessing OncoKB annotations. Please provide it using the --oncokb-api-key flag or set ONCOKB_API_KEY in your environment.")
        exit(1)
    elif not oncokb_api_key:
        print("OncoKB API key is required for accessing OncoKB annotations. Please provide it using the --oncokb-api-key flag or set ONCOKB_API_KEY in your environment.")
        exit(1)

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
        processes_list = [p.strip() for p in processes.split(',')] if processes else []
        samples_list = [s.strip() for s in samples.split(',')] if samples else []

        workdirs_processes = {}
        if processes_list:
            print("Retrieving work directories for specified processes...")
            entries = get_entries_with_process_names(processes_list)
            for entry in entries:
                workdir = entry['workdir']
                name = entry['name']
                if workdir not in workdirs_processes:
                    workdirs_processes[workdir] = set()
                workdirs_processes[workdir].add(name)

        workdirs_samples = {}
        if samples_list:
            print("Retrieving work directories for specified samples...")
            entries = get_entries_with_sample_names(samples_list)
            for entry in entries:
                workdir = entry['workdir']
                name = entry['name']
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
                names_str = ', '.join(sorted(names))
                print(f"{names_str}, {workdir}")

            confirm = input("Do you want to delete these directories? ([y]es/[n]o): ").strip().lower()
            if confirm == 'yes' or confirm == 'y':
                confirm_again = input("Are you sure? This will cause the pipeline to rerun from these steps. ([y]es/[n]o): ").strip().lower()
                if confirm_again == 'yes' or confirm_again == 'y':
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

    runner = NextflowRunner()

    command = (
        f"{load_modules_command} "
        f"{runner.cmd} secrets set ONCOKB_API_KEY {oncokb_api_key} && "
        f"{runner.cmd} -log .nextflow_{runner.get_timestamp()}.log "
        f"run {pipeline_dir} "
        f"-params-file {params_file} "
        f"-profile {profile} "
        f"-with-report report_{runner.get_timestamp()}.html "
        f"-with-trace"
    )

    if resume:
        command += " -resume"

    runner.run(command)


from ..core.outputs import OUTPUT_KEYS
@run_cli.command()
@click.option('-p', '--pipeline-output-dir', required=True, type=click.Path(exists=True), help="Directory containing pipeline outputs")
@click.option('-s', '--samplesheet', default='./samplesheet.csv', required=True, type=click.Path(exists=True), help="Path to the samplesheet CSV file")
@click.option('--old', is_flag=True, default=False, help="Whether to use the old outputs mapping (default: False)")
@click.option('--prefer-outputs', is_flag=True, default=False, help="If output and samplesheet column are both found, use the output column (default: False)")
@click.option('-o', '--output', type=click.Path(), help="CSV file to save outputs (default: stdout)")
@click.option('-c', '--include-columns', help='Comma-separated list of columns to include. Available: {}'.format(",".join(OUTPUT_KEYS)))
@click.option('-C', '--exclude-columns', help='Comma-separated list of columns to exclude. Available: see --include-columns')
def outputs(pipeline_output_dir, samplesheet, old, prefer_outputs, output, include_columns, exclude_columns):
    """
    Generate an outputs.csv file suitable for skilifting.
    Reads pipeline outputs and samplesheet to create a CSV mapping patient data to output file paths.
    """
    from ..core.outputs import Outputs, OUTPUT_KEYS # Import OUTPUT_KEYS for help text formatting

    if include_columns and exclude_columns:
        raise click.UsageError("Options -c/--include-columns and -C/--exclude-columns are mutually exclusive.")

    include_list = include_columns.split(',') if include_columns else None
    exclude_list = exclude_columns.split(',') if exclude_columns else None

    # Validate column names if provided
    available_columns = set(OUTPUT_KEYS)
    if include_list:
        invalid_cols = [col for col in include_list if col not in available_columns]
        if invalid_cols:
            raise click.BadParameter(f"Invalid columns specified in --include-columns: {', '.join(invalid_cols)}. Available: {', '.join(OUTPUT_KEYS)}")
    if exclude_list:
        invalid_cols = [col for col in exclude_list if col not in available_columns]
        if invalid_cols:
            raise click.BadParameter(f"Invalid columns specified in --exclude-columns: {', '.join(invalid_cols)}. Available: {', '.join(OUTPUT_KEYS)}")


    outputs_obj = Outputs(pipeline_output_dir, samplesheet, old, prefer_outputs)
    outputs_obj.emit_output_csv(output, include_columns=include_list, exclude_columns=exclude_list)
    if output:
        click.echo(f"Outputs CSV generated at: {output}")


from ..core.outputs import SAMPLESHEET_FIELDNAMES
@run_cli.command()
@click.option('-p', '--pipeline-output-dir', required=True, type=click.Path(exists=True), help="Directory containing pipeline outputs")
@click.option('-s', '--samplesheet', default='./samplesheet.csv', required=True, type=click.Path(exists=True), help="Path to the samplesheet CSV file")
@click.option('--old', is_flag=True, default=False, help="Whether to use the old outputs mapping (default: False)")
@click.option('-o', '--output', type=click.Path(), help="CSV file to save samplesheet (default: stdout)")
@click.option('-c', '--include-columns', help='Comma-separated list of columns to include. Available: {}'.format(",".join(SAMPLESHEET_FIELDNAMES)))
@click.option('-C', '--exclude-columns', help='Comma-separated list of columns to exclude. Available: see --include-columns')
def samplesheet(pipeline_output_dir, samplesheet, old, output, include_columns, exclude_columns):
    """
    Generate a samplesheet.csv file suitable for a pipeline run.
    Reads pipeline outputs and the original samplesheet to create a new samplesheet
    with paths to generated files filled in.
    """
    from ..core.outputs import Outputs

    if include_columns and exclude_columns:
        raise click.UsageError("Options -c/--include-columns and -C/--exclude-columns are mutually exclusive.")

    include_list = include_columns.split(',') if include_columns else None
    exclude_list = exclude_columns.split(',') if exclude_columns else None

    # Validate column names if provided
    available_columns = set(SAMPLESHEET_FIELDNAMES)
    if include_list:
        invalid_cols = [col for col in include_list if col not in available_columns]
        if invalid_cols:
            raise click.BadParameter(f"Invalid columns specified in --include-columns: {', '.join(invalid_cols)}. Available: {', '.join(SAMPLESHEET_FIELDNAMES)}")
    if exclude_list:
        invalid_cols = [col for col in exclude_list if col not in available_columns]
        if invalid_cols:
            raise click.BadParameter(f"Invalid columns specified in --exclude-columns: {', '.join(invalid_cols)}. Available: {', '.join(SAMPLESHEET_FIELDNAMES)}")

    outputs_obj = Outputs(pipeline_output_dir, samplesheet, old)
    outputs_obj.emit_samplesheet_csv(output, include_columns=include_list, exclude_columns=exclude_list)
    if output:
        click.echo(f"Samplesheet CSV generated at: {output}")


@run_cli.command()
@click.option('-p', '--pipeline-output-dir', type=click.Path(exists=True), help="Directory containing pipeline outputs")
@click.option('-s', '--samplesheet', type=click.Path(exists=True), help="Path to the samplesheet CSV file used to generate the outputs.csv")
@click.option('--old', is_flag=True, default=False, help="Whether the outputs.csv was generated using the old process_name/patient_id nf-gos outputs mapping (default: False)")
@click.option('--output-csv', type=click.Path(exists=True), help="Path to the outputs.csv file generated by 'gosh run outputs'.")
@click.option('-t', '--cohort-type', type=click.Choice(['paired', 'tumor_only', 'heme']), default='paired', help='Type of the cohort.')
@click.option('-o', '--gos_dir', type=click.Path(), required=True, help='Path to where the skilifted outputs should be deposited.')
@click.option('-l', '--skilift-repo', type=click.Path(), default="~/git/skilift", help='Path to the skilift repo (default: ~/git/skilift)')
@click.option('-c', '--cores', type=int, default=1, help='Number of cores to use.')
def skilift(
    cohort_type,
    gos_dir,
    skilift_repo,
    cores,
    old,
    output_csv=None,
    pipeline_output_dir=None,
    samplesheet=None
):
    """Lift raw data into gOS compatible formats using Skilift."""
    from os import path
    import subprocess
    # Need these imports here for the help text formatting in the options above
    from ..core.outputs import OUTPUT_KEYS, SAMPLESHEET_FIELDNAMES

    # Expand the skilift_repo path
    skilift_repo = path.expanduser(skilift_repo)
    if not path.isdir(skilift_repo):
        if click.confirm(f"Skilift repository not found at {skilift_repo}. Would you like to clone it?", default=True):
            clone_url = "https://github.com/mskilab-org/skilift.git"
            click.echo(f"Cloning skilift from {clone_url} to {skilift_repo} ...")
            try:
                subprocess.run(["git", "clone", clone_url, skilift_repo], check=True)
            except Exception as err:
                click.echo(f"Error cloning skilift repository: {err}")
                return
        else:
            click.echo("Skilift repository is required. Exiting.")
            return

    if not pipeline_output_dir and not output_csv:
        click.echo("Error: Either --pipeline-output-dir or --output-csv must be provided.")
        return

    if not output_csv:
        if not path.isfile(samplesheet):
            click.echo(f"Error: The samplesheet '{samplesheet}' does not exist. Please provide a samplesheet with the -s option.")
            return
        if not path.isdir(pipeline_output_dir):
            click.echo(f"Error: The directory '{pipeline_output_dir}' does not exist. Please provide a valid pipeline output directory with the -p option.")
            return

        # Generate the outputs.csv file if not provided
        output_csv = './outputs.csv'
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
        click.echo(f"Error: The outputs CSV file '{output_csv}' does not exist or could not be generated.")
        return

    makedirs(path.expanduser(gos_dir), exist_ok=True)

    r_code = f'''
    devtools::load_all("{path.expanduser(skilift_repo)}")
    cohort <- Cohort$new("{output_csv}", cohort_type="{cohort_type}")
    saveRDS(cohort, "{gos_dir}/cohort.rds")
    cohort_mod <- lift_all(cohort, output_data_dir="{gos_dir}", cores={cores})
    message("Saving over cohort")
    saveRDS(cohort_mod, "{gos_dir}/cohort.rds")
    '''

    try:
        run(['Rscript', '-e', r_code], check=True)
    except CalledProcessError as e:
        click.echo(f"Error executing R script: {e}")


@run_cli.command()
@click.option('-d', '--datasets-json',
              type=click.Path(),
              help='Path to the datasets.json file. Defaults to value from environment config.')
@click.option('-g', '--gos-data-dir',
              required=True,
              type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True),
              help='Path to the directory containing skilifted gOS data and datafiles.json.')
@click.option('-t', '--title',
              required=True,
              type=str,
              help='Title for the new dataset cohort.')
def datasets(datasets_json, gos_data_dir, title):
    """
    Add a new dataset entry to the datasets.json file for gOS webapp consumption.
    """
    env_defaults = get_environment_defaults()
    default_datasets_json = env_defaults.get('datasets_json')

    if not datasets_json:
        if not default_datasets_json:
            click.echo("Error: --datasets-json path must be provided or configured in environment defaults.", err=True)
            exit(1)
        datasets_json = default_datasets_json
        click.echo(f"Using default datasets.json path from environment: {datasets_json}")
    else:
        # Ensure the directory for the datasets.json exists if a path is provided
        datasets_dir = os.path.dirname(os.path.abspath(datasets_json))
        os.makedirs(datasets_dir, exist_ok=True)


    # Validate gos_data_dir contains datafiles.json
    datafiles_path = os.path.join(gos_data_dir, 'datafiles.json')
    if not os.path.isfile(datafiles_path):
        click.echo(f"Error: 'datafiles.json' not found in the specified gOS data directory: {gos_data_dir}", err=True)
        exit(1)

    # Generate ID from title
    dataset_id = title.lower().replace(' ', '_')

    # Convert gos_data_dir to the required URL format
    try:
        data_path_url = convert_to_datasets_path(gos_data_dir)
    except ValueError as e:
        click.echo(f"Error converting gOS data directory path: {e}", err=True)
        click.echo("Please ensure the path follows the expected structure (e.g., /path/to/{lab_name}/external/{rest_of_path})", err=True)
        exit(1)

    datafiles_path_url = data_path_url + "datafiles.json"

    # Create the new dataset entry
    new_entry = {
        "id": dataset_id,
        "title": title,
        "dataPath": data_path_url,
        "datafilesPath": datafiles_path_url
    }

    # Read existing datasets.json or initialize if it doesn't exist
    all_datasets = []
    if os.path.exists(datasets_json):
        try:
            with open(datasets_json, 'r') as f:
                all_datasets = json.load(f)
                if not isinstance(all_datasets, list):
                    click.echo(f"Error: Expected {datasets_json} to contain a JSON list.", err=True)
                    exit(1)
        except json.JSONDecodeError:
            click.echo(f"Error: Could not decode JSON from {datasets_json}. Initializing with new entry.", err=True)
            all_datasets = [] # Reset if file is corrupt
        except Exception as e:
            click.echo(f"Error reading {datasets_json}: {e}", err=True)
            exit(1)

    # Check if an entry with the same ID or title already exists
    existing_ids = {entry.get('id') for entry in all_datasets if 'id' in entry}
    existing_titles = {entry.get('title') for entry in all_datasets if 'title' in entry}

    if new_entry['id'] in existing_ids:
        click.echo(f"Error: Dataset with ID '{new_entry['id']}' already exists in {datasets_json}.", err=True)
        exit(1)
    if new_entry['title'] in existing_titles:
        click.echo(f"Warning: Dataset with title '{new_entry['title']}' already exists in {datasets_json}. Proceeding with unique ID '{new_entry['id']}'.", err=True)


    # Append the new entry
    all_datasets.append(new_entry)

    # Write the updated list back to datasets.json
    try:
        with open(datasets_json, 'w') as f:
            json.dump(all_datasets, f, indent=4)
        click.echo(f"Successfully added dataset '{title}' (ID: {dataset_id}) to {datasets_json}")
        click.echo(f"  dataPath: {data_path_url}")
        click.echo(f"  datafilesPath: {datafiles_path_url}")
    except Exception as e:
        click.echo(f"Error writing updated data to {datasets_json}: {e}", err=True)
        exit(1)
