"""Workflow interface for PyPSA-USA that wraps Snakemake functionality."""

from pathlib import Path

import snakemake


def setup_workflow_directories(base_dir=None):
    """
    Set up the workflow directories for data, resources, results, logs, and benchmarks.

    Parameters
    ----------
    base_dir : str or Path, optional
        Base directory for the workflow. If None, uses the current working directory.

    Returns
    -------
    dict
        Dictionary containing the paths to the workflow directories.
    """
    if base_dir is None:
        base_dir = Path.cwd()
    else:
        base_dir = Path(base_dir)

    # Define directory structure
    dirs = {
        "data": base_dir / "data",
        "resources": base_dir / "resources",
        "results": base_dir / "results",
        "logs": base_dir / "logs",
        "benchmarks": base_dir / "benchmarks",
        "cutouts": base_dir / "cutouts",
    }

    # Create directories if they don't exist
    for dir_path in dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)

    return dirs


def run_workflow(
    config_file=None,
    until=None,
    cores=1,
    dryrun=False,
    base_dir=None,
    target="all",
):
    """
    Run the PyPSA-USA workflow using Snakemake.

    Parameters
    ----------
    config_file : str or Path, optional
        Path to the config file. If None, uses default config.
    target : str, optional
        Target rule or file to generate.
    cores : int, default=1
        Number of CPU cores to use.
    dryrun : bool, default=False
        If True, only show what would be done without executing.
    base_dir : str or Path, optional
        Base directory for the workflow. If None, uses the current working directory.

    Returns
    -------
    bool
        True if the workflow completed successfully.
    """
    # Set up directories
    setup_workflow_directories(base_dir)

    # Get the path to the Snakefile
    workflow_dir = Path(__file__).parent.parent
    snakefile = workflow_dir / "Snakefile"

    if not snakefile.exists():
        raise FileNotFoundError(f"Snakefile not found at {snakefile}")

    # Set up Snakemake arguments
    snakemake_args = {
        "snakefile": str(snakefile),
        "cores": cores,
        "dryrun": dryrun,
        "workdir": str(workflow_dir),
        "configfiles": [str(config_file)] if config_file else None,
        "targets": [target],
        "until": [until] if until else None,
    }
    # Run the workflow
    return snakemake.snakemake(**snakemake_args)


def load_network(config_file=None, cores=1, dryrun=False, base_dir=None):
    """
    Prepare the network using the workflow.

    Parameters
    ----------
    config_file : str or Path, optional
        Path to the config file. If None, uses default config.
    cores : int, default=1
        Number of CPU cores to use.
    dryrun : bool, default=False
        If True, only show what would be done without actually doing it.
    base_dir : str or Path, optional
        Base directory for the workflow. If None, uses the current working directory.

    Returns
    -------
    bool
        True if the network was prepared successfully.
    """
    return run_workflow(
        config_file=config_file,
        until="prepare_network",
        cores=cores,
        dryrun=dryrun,
        base_dir=base_dir,
    )


def solve_network(config_file=None, cores=1, dryrun=False, base_dir=None):
    """
    Solve the network using the workflow.

    Parameters
    ----------
    config_file : str or Path, optional
        Path to the config file. If None, uses default config.
    cores : int, default=1
        Number of CPU cores to use.
    dryrun : bool, default=False
        If True, only show what would be done without actually doing it.
    base_dir : str or Path, optional
        Base directory for the workflow. If None, uses the current working directory.

    Returns
    -------
    bool
        True if the network was solved successfully.
    """
    return run_workflow(
        config_file=config_file,
        until="solve_network",
        cores=cores,
        dryrun=dryrun,
        base_dir=base_dir,
    )
