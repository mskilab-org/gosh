import click
from .cli.run import run_cli
from .cli.debug import debug_cli
from .cli.help import help_cli

@click.group()
def cli():
    """gOSh - gOS sHell"""
    pass

# Register command groups
cli.add_command(run_cli, name='run')
cli.add_command(debug_cli, name='debug')
cli.add_command(help_cli, name='help')
