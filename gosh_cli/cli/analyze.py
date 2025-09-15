# gosh_cli/cli/analyze.py

import click
import json
from pathlib import Path
from ..core.pipeline_analyzer import ProcessDependencyAnalyzer

@click.group(name='analyze')
def analyze_cli():
    """Analyze pipeline structure and dependencies"""
    pass

@analyze_cli.command()
@click.option('--pipeline-dir', 
              type=click.Path(exists=True),
              help='Path to nf-gOS pipeline directory')
@click.option('-o', '--output', 
              default='pipeline_dependencies.json',
              help='Output file for dependency analysis')
@click.option('--format', 
              type=click.Choice(['json', 'yaml', 'both']),
              default='both',
              help='Output format')
@click.option('--show-tree/--no-tree',
              default=True,
              help='Display workflow hierarchy tree')
@click.option('--show-groups/--no-groups',
              default=False,
              help='Display execution groups')
def dependencies(pipeline_dir, output, format, show_tree, show_groups):
    """Analyze process dependencies and workflow structure"""
    from ..core.module_loader import get_environment_defaults
    
    # Use environment default if not specified
    if not pipeline_dir:
        env_defaults = get_environment_defaults()
        pipeline_dir = env_defaults.get('pipeline-dir')
        if not pipeline_dir:
            raise click.ClickException(
                "Pipeline directory not specified. Use --pipeline-dir or run from an environment with defaults."
            )
    
    click.echo(f"Analyzing pipeline at: {pipeline_dir}")
    
    # Run analysis
    analyzer = ProcessDependencyAnalyzer(pipeline_dir)
    analyzer.analyze_all()
    
    # Display summary
    if show_tree or show_groups:
        analyzer.print_summary(show_tree=show_tree, show_groups=show_groups)
    
    # Save results
    config = analyzer.generate_dependency_config()
    
    if format in ['json', 'both']:
        with open(output, 'w') as f:
            json.dump(config, f, indent=2)
        click.echo(f"Saved JSON to {output}")
    
    if format in ['yaml', 'both']:
        import yaml
        yaml_file = output.replace('.json', '.yaml')
        with open(yaml_file, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
        click.echo(f"Saved YAML to {yaml_file}")
    
    # Generate schema extension
    schema_ext = {}
    field_mappings = config.get('field_mappings', {})
    for field, info in field_mappings.items():
        schema_ext[field] = {
            'produced_by': [p['name'] for p in info.get('produced_by', [])],
            'required_by': [p['name'] for p in info.get('required_by', [])],
        }
    
    schema_file = 'schema_dependencies.json'
    with open(schema_file, 'w') as f:
        json.dump(schema_ext, f, indent=2)
    click.echo(f"Saved schema dependencies to {schema_file}")

@analyze_cli.command()
@click.option('--pipeline-dir',
              type=click.Path(exists=True),
              help='Path to nf-gOS pipeline directory')
@click.option('--field',
              required=True,
              help='Samplesheet field name to trace')
@click.option('--direction',
              type=click.Choice(['upstream', 'downstream', 'both']),
              default='both',
              help='Direction to trace dependencies')
def trace(pipeline_dir, field, direction):
    """Trace dependencies for a specific samplesheet field"""
    from ..core.module_loader import get_environment_defaults
    
    # Use environment default if not specified
    if not pipeline_dir:
        env_defaults = get_environment_defaults()
        pipeline_dir = env_defaults.get('pipeline-dir')
        if not pipeline_dir:
            raise click.ClickException(
                "Pipeline directory not specified. Use --pipeline-dir or run from an environment with defaults."
            )
    
    # Run analysis
    analyzer = ProcessDependencyAnalyzer(pipeline_dir)
    analyzer.analyze_all()
    config = analyzer.generate_dependency_config()
    
    field_mappings = config.get('field_mappings', {})
    
    if field not in field_mappings:
        click.echo(f"Field '{field}' not found in pipeline")
        click.echo("\nAvailable fields:")
        for f in sorted(field_mappings.keys()):
            click.echo(f"  - {f}")
        return
    
    info = field_mappings[field]
    
    click.echo(f"\n=== Dependency trace for '{field}' ===")
    
    if direction in ['upstream', 'both']:
        click.echo(f"\nProduced by:")
        producers = info.get('produced_by', [])
        if producers:
            for p in producers:
                click.echo(f"  - [{p['type']}] {p['name']}")
        else:
            click.echo("  - (primary input - not produced by pipeline)")
    
    if direction in ['downstream', 'both']:
        click.echo(f"\nRequired by:")
        consumers = info.get('required_by', [])
        if consumers:
            for c in consumers:
                click.echo(f"  - [{c['type']}] {c['name']}")
        else:
            click.echo("  - (not required by any process)")

@analyze_cli.command()
@click.option('--pipeline-dir',
              type=click.Path(exists=True),
              help='Path to nf-gOS pipeline directory')
@click.option('--process',
              required=True,
              help='Process name to inspect')
def process(pipeline_dir, process):
    """Show detailed information about a specific process"""
    from ..core.module_loader import get_environment_defaults
    
    # Use environment default if not specified
    if not pipeline_dir:
        env_defaults = get_environment_defaults()
        pipeline_dir = env_defaults.get('pipeline-dir')
        if not pipeline_dir:
            raise click.ClickException(
                "Pipeline directory not specified. Use --pipeline-dir or run from an environment with defaults."
            )
    
    # Run analysis
    analyzer = ProcessDependencyAnalyzer(pipeline_dir)
    analyzer.analyze_all()
    
    if process not in analyzer.processes:
        click.echo(f"Process '{process}' not found")
        click.echo("\nAvailable processes:")
        for p in sorted(analyzer.processes.keys()):
            click.echo(f"  - {p}")
        return
    
    proc_info = analyzer.processes[process]
    
    click.echo(f"\n=== Process: {process} ===")
    click.echo(f"File: {proc_info['file']}")
    
    if proc_info.get('directives'):
        click.echo("\nDirectives:")
        for key, val in proc_info['directives'].items():
            click.echo(f"  {key}: {val}")
    
    click.echo("\nInputs:")
    inputs = proc_info['inputs']
    if inputs['files']:
        click.echo("  Files:")
        for f in inputs['files']:
            click.echo(f"    - {f['name']} ({f['type']})")
    if inputs['values']:
        click.echo("  Values:")
        for v in inputs['values']:
            click.echo(f"    - {v}")
    if inputs['tuples']:
        click.echo("  Tuples:")
        for t in inputs['tuples']:
            components = ', '.join([f"{c['type']}({c['name']})" for c in t])
            click.echo(f"    - [{components}]")
    
    click.echo("\nOutputs:")
    outputs = proc_info['outputs']
    if outputs['files']:
        click.echo("  Files:")
        for f in outputs['files']:
            click.echo(f"    - {f['pattern']} ({f['type']})")
    if outputs['emits']:
        click.echo("  Emits:")
        for name, val in outputs['emits'].items():
            click.echo(f"    - {name} = {val}")
    
    if proc_info.get('conditions'):
        click.echo(f"\nWhen condition: {proc_info['conditions']}")
    
    deps = proc_info.get('dependencies', {})
    if deps.get('depends_on'):
        click.echo("\nDepends on:")
        for d in deps['depends_on']:
            click.echo(f"  - {d['name']} (via {d.get('field', 'unknown')})")

# Add to main.py:
# from .cli.analyze import analyze_cli
# cli.add_command(analyze_cli, name='analyze')