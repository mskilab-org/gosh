import click
from .cli.run import run_cli
from .cli.debug import debug_cli
from .cli.help import help_cli
from gosh_cli import __version__

@click.group()
@click.version_option(version=__version__, prog_name="gosh")
def cli():
    """gOSh - gOS sHell"""
    pass

# Register command groups
cli.add_command(run_cli, name='run')
cli.add_command(debug_cli, name='debug')
cli.add_command(help_cli, name='help')
