#!/usr/bin/env python3

from os.path import exists
from csv import DictWriter
from sys import stdout
import click

@click.group(name='debug')
def debug_cli():
    """Debug commands for analyzing Nextflow runs."""
    pass

@debug_cli.command()
@click.argument('log_file', type=click.Path(exists=True), required=False)
def eye(log_file):
    """Analyze Nextflow log file using AI assistance."""
    try:
        # If no log file specified, look for .nextflow.log in current directory
        if not log_file:
            default_log = '.nextflow.log'
            if exists(default_log):
                log_file = default_log
            else:
                click.secho("Error: No .nextflow.log file found in current directory.", fg='red')
                click.secho("Please specify a log file path or run from a directory containing .nextflow.log", fg='yellow')
                return

        # Read the log file
        with open(log_file, 'r') as f:
            log_content = f.read()

        from ..utils.ai_helper import extract_error_messages, get_error_analysis_and_solution
        # Extract error messages
        error_messages = extract_error_messages(log_content)

        if not error_messages:
            click.secho("No errors found in the log file.", fg='green')
            return

        # Get AI analysis and solution
        analysis = get_error_analysis_and_solution(error_messages)

        # Print the analysis with colors
        click.secho("\n=== AI Analysis of Nextflow Errors ===", fg='blue', bold=True)
        click.echo("\n" + analysis)

    except Exception as e:
        click.secho(f"Error analyzing log file: {str(e)}", fg='red')

@debug_cli.command()
@click.option('-s', '--sample_names', type=str, help='Comma-separated list of sample ID(s).')
@click.option('-p', '--process_names', type=str, help='Comma-separated list of process name(s) (in ALL CAPS).')
@click.option('-o', '--output', type=click.Path(), help='Output file to save the results as CSV.')
@click.option('--refresh', is_flag=True, help='Refresh log cache')
@click.option('--status', type=str, help='Comma-separated list of statuses to filter by (e.g., FAILED,COMPLETED, CACHED, ABORTED).')
@click.option('--work-dir-exists', is_flag=True, help='Filter results to only those entries whose workdir exists.')
@click.option('-c', '--columns', type=str, help='Comma-separated list of columns to filter by. Valid columns: name, process, workdir, status, cpus, pcpu, memory, pmem')
def log(sample_names, process_names, output, refresh, status, work_dir_exists, columns):
    """Retrieve Nextflow log entries based on samples or processes."""

    from ..core.nextflow_log import (
        get_all_entries,
        write_cache,
        read_cache
    )

    try:
        cache_file = ".gosh_debug_log_cache.csv"
        # Determine whether to (re)generate the cache
        if refresh or not exists(cache_file):
            click.secho("Generating log cache...", fg='blue')
            entries = get_all_entries()
            # Save the cache to file
            write_cache(entries, cache_file)
        else:
            click.secho("Using cached log entries...", fg='blue')
            entries = read_cache(cache_file)

        if sample_names:
            sample_list = [s.strip() for s in sample_names.split(',')]
            click.secho(f"Filtering entries for sample name(s): {', '.join(sample_list)}", fg='blue')
            entries_sample = [entry for entry in entries if any(sample in entry.get('name', '') for sample in sample_list)]
        else:
            entries_sample = entries

        if process_names:
            process_list = [p.strip() for p in process_names.split(',')]
            click.secho(f"Filtering entries for process name(s): {', '.join(process_list)}", fg='blue')
            entries_process = [entry for entry in entries if any(proc in entry.get('process', '') for proc in process_list)]
        else:
            entries_process = entries

        if sample_names and process_names:
            # Compute the intersection based on tuple conversion
            set_sample = set(tuple(entry.items()) for entry in entries_sample)
            set_process = set(tuple(entry.items()) for entry in entries_process)
            entries = [dict(item) for item in (set_sample & set_process)]
        elif sample_names:
            entries = entries_sample
        elif process_names:
            entries = entries_process

        if status:
            statuses = [s.strip().upper() for s in status.split(',')]
            click.secho(f"Filtering entries for statuses: {', '.join(statuses)}", fg='blue')
            entries = [entry for entry in entries if entry.get('status', '').upper() in statuses]

        if work_dir_exists:
            click.secho("Filtering entries for existing work directories.", fg='blue')
            entries = [entry for entry in entries if exists(entry.get('workdir', ''))]

        if not entries:
            click.secho("No log entries found.", fg='yellow')
            return

        # Remove duplicate entries if both sample and process names overlap
        unique_entries = [dict(t) for t in {tuple(d.items()) for d in entries}]

        # If columns flag is provided, filter the keys accordingly
        allowed_fields = ['name', 'process', 'workdir', 'status', 'cpus', 'pcpu', 'memory', 'pmem']
        if columns:
            columns_list = [col.strip() for col in columns.split(',')]
            # Validate that all columns are allowed
            invalid = [col for col in columns_list if col not in allowed_fields]
            if invalid:
                click.secho(f"Error: Invalid column(s) specified: {', '.join(invalid)}", fg='red')
                return
            fieldnames = columns_list
            # For each entry in unique_entries, filter only the allowed keys in the defined order
            unique_entries = [
                {col: entry.get(col, '') for col in columns_list}
                for entry in unique_entries
            ]
        else:
            fieldnames = list(unique_entries[0].keys())

        if output:
            with open(output, 'w', newline='') as csvfile:
                writer = DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(unique_entries)
            click.secho(f"Results saved to {output}", fg='green')
        else:
            writer = DictWriter(stdout, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(unique_entries)

    except Exception as e:
        click.secho(f"Error retrieving log entries: {str(e)}", fg='red')
