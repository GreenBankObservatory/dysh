import os
import sys
from importlib.resources import files

from jupyter_server.serverapp import ServerApp
from traitlets.config import Config


def main():
    templates_dir = files("dysh") / "lab_templates"
    shell_dir = files("dysh.shell")

    # Add notebooks directory so examples subdirectory shows as templates
    from pathlib import Path

    import dysh

    notebooks_dir = Path(dysh.__file__).parent.parent.parent / "notebooks"

    # Set environment variable so jupyter_app_launcher can find our config
    # The extension searches for jp_app_launcher*.yaml in JUPYTER_APP_LAUNCHER_PATH
    os.environ["JUPYTER_APP_LAUNCHER_PATH"] = str(shell_dir)

    c = Config()
    c.ServerApp.jpserver_extensions = {
        "jupyterlab": True,
        "jupyterlab_templates": True,
        "jupyter_app_launcher": True,
        # "jupyter_lsp": False,  # Disable LSP to avoid basedpyright errors
    }

    # Add both lab_templates and notebooks to template dirs
    # jupyterlab-templates will find subdirectories (dysh, examples) and show them as categories
    template_paths = [str(templates_dir)]
    if notebooks_dir.exists():
        template_paths.append(str(notebooks_dir))

    c.JupyterLabTemplates.template_dirs = template_paths
    c.JupyterLabTemplates.include_default = False

    # Default to JupyterLab interface
    c.ServerApp.default_url = "/lab"

    ServerApp.launch_instance(config=c, argv=sys.argv[1:])


if __name__ == "__main__":
    main()
