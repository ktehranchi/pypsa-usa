"""Test script for PyPSA-USA workflow interface."""

from pathlib import Path

import pypsa_usa


def main():
    # Set up base directory for the workflow
    base_dir = Path.cwd()

    # Get the default config file
    config_file = Path("workflow/config/config.default.yaml")

    print("Testing PyPSA-USA workflow interface...")
    print("\n1. First, let's see what would happen (dryrun):")

    print("\nPreparing network (dryrun):")
    pypsa_usa.prepare_network(
        config_file=config_file,
        dryrun=True,
        cores=1,
        base_dir=base_dir,
    )

    print("\nSolving network (dryrun):")
    pypsa_usa.solve_network(
        config_file=config_file,
        dryrun=True,
        cores=1,
        base_dir=base_dir,
    )

    print("\nGenerating plots (dryrun):")
    pypsa_usa.plot_results(
        config_file=config_file,
        dryrun=True,
        cores=1,
        base_dir=base_dir,
    )

    response = input("\nWould you like to run the actual workflow? (y/n): ")

    if response.lower() == "y":
        print("\n2. Running the actual workflow:")

        print("\nPreparing network...")
        success = pypsa_usa.prepare_network(
            config_file=config_file,
            cores=1,
            base_dir=base_dir,
        )
        print(f"Network preparation {'succeeded' if success else 'failed'}")

        print("\nSolving network...")
        success = pypsa_usa.solve_network(
            config_file=config_file,
            cores=1,
            base_dir=base_dir,
        )
        print(f"Network solution {'succeeded' if success else 'failed'}")

        print("\nGenerating plots...")
        success = pypsa_usa.plot_results(
            config_file=config_file,
            cores=1,
            base_dir=base_dir,
        )
        print(f"Plot generation {'succeeded' if success else 'failed'}")

    print("\nTest complete!")


if __name__ == "__main__":
    main()
