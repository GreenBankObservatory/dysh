import argparse
import os
import sys
from importlib.resources import files

from jupyter_server.serverapp import ServerApp
from traitlets.config import Config

from dysh import __version__
from dysh.util.core import get_project_root


def parse_args():
    parser = argparse.ArgumentParser(description=("Dysh lab"))
    parser.add_argument("--version", help="Print version and exit", action="store_true")
    parser.add_argument("--no-browser", help="Don't open browser automatically", action="store_true")
    return parser.parse_known_args()


def main():
    args, remaining_args = parse_args()
    if args.version:
        print(f"dysh-lab: v{__version__}")
        sys.exit(0)

    template_notebooks = files("dysh") / "lab_templates"
    shell_dir = files("dysh.shell")

    example_notebooks = get_project_root().parent.parent / "notebooks"

    # Set environment variable so jupyter_app_launcher can find our config
    # The extension searches for jp_app_launcher*.yaml in JUPYTER_APP_LAUNCHER_PATH
    os.environ["JUPYTER_APP_LAUNCHER_PATH"] = str(shell_dir)

    c = Config()

    c.ServerApp.jpserver_extensions = {
        "jupyterlab": True,
        "jupyterlab_templates": True,
        "jupyter_app_launcher": True,
    }

    # Add both lab_templates and notebooks to template dirs
    # jupyterlab-templates will find subdirectories (dysh, examples) and show them as categories
    c.JupyterLabTemplates.template_dirs = [str(template_notebooks), str(example_notebooks)]
    c.JupyterLabTemplates.include_default = False

    # Default to JupyterLab interface
    c.ServerApp.default_url = "/lab"

    # Open browser by default unless --no-browser is specified
    c.ServerApp.open_browser = not args.no_browser

    ServerApp.launch_instance(config=c, argv=remaining_args)


if __name__ == "__main__":
    main()
