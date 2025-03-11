"""Script used to compare outputs from multiple snakemake scenarios."""

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import yaml
from matplotlib import pyplot as plt


class ScenarioDataGetter:
    """Class for loading and processing scenario data."""

    def __init__(self, yaml_path: str | Path, force_regenerate: bool = False):
        """
        Initialize with configuration from YAML file.

        Parameters
        ----------
        yaml_path : str or Path
            Path to the YAML configuration file
        force_regenerate : bool, default False
            Whether to force regeneration of statistics files even if they exist
        """
        self.config = self._load_yaml_config(yaml_path)
        self.scenarios = self.config["scenarios"]
        self.alias_dict = self.config.get("alias_dict", None)
        self.new_order = self.config.get("new_order", None)
        self.reference_scenario = self.config.get("reference_scenario", None)
        self.force_regenerate = force_regenerate

        # Set output path
        output_folder = self.config.get("output_folder_name", "scenario_comparison")
        self.figures_path = Path.cwd() / f"results/{output_folder}"
        self.figures_path.mkdir(exist_ok=True)

        # Create cache directory
        self.cache_path = self.figures_path / "cached_statistics"
        self.cache_path.mkdir(exist_ok=True)

        # Load data
        self.raw_data = self._load_scenario_data()
        self.processed_data = self._process_data()

        # Load network
        self.network = pypsa.Network(self.config["network"]["path"])
        self.carriers = self._get_carriers()

    def _load_yaml_config(self, yaml_path: str | Path) -> dict:
        """Load configuration from YAML file."""
        with open(yaml_path) as file:
            return yaml.safe_load(file)

    def _load_scenario_data(self) -> dict[str, dict[str, pd.DataFrame]]:
        """
        Load CSV data for all scenarios, generating statistics first if needed.

        If cached statistics files exist, they are loaded instead of regenerating.
        """
        data = {}
        for scenario in self.scenarios:
            scenario_name = scenario["name"]
            path = Path(scenario["path"])

            # Check if cached statistics exist
            cached_file = self.cache_path / f"statistics_{scenario_name}.csv"

            # Generate statistics if forced or if cache doesn't exist
            if self.force_regenerate or not cached_file.exists():
                print(f"Generating statistics for {scenario_name}...")
                start_time = time.time()

                # Load the network for this scenario
                network_path = path
                if not network_path.exists():
                    print(f"Network file not found at {network_path}. Trying to find based on wildcard...")
                    # Try to find network file with a wildcard pattern
                    network_files = list(path.glob("networks/elec_*.nc"))
                    if network_files:
                        network_path = network_files[0]
                    else:
                        raise FileNotFoundError(f"Cannot find network file for {scenario_name}")

                print(f"Loading network from {network_path}")
                n = pypsa.Network(network_path)

                # Generate statistics
                stats = self._generate_statistics(n)
                # Save to cache
                stats.to_csv(cached_file)

                elapsed = time.time() - start_time
                print(f"Statistics for {scenario_name} generated in {elapsed:.1f} seconds and saved to cache")

                data[scenario_name] = {"statistics": stats}
            else:
                print(f"Loading cached statistics for {scenario_name}...")
                data[scenario_name] = {"statistics": pd.read_csv(cached_file, index_col=[0, 1], header=[0, 1])}

        return data

    def _generate_statistics(self, n: pypsa.Network) -> pd.DataFrame:
        """
        Generate statistics from a PyPSA network.

        Parameters
        ----------
        n : pypsa.Network
            Network to analyze

        Returns
        -------
        pd.DataFrame
            DataFrame containing statistics
        """
        return n.statistics()

    def _process_data(self) -> dict[str, dict[str, pd.DataFrame]]:
        """Process data to match the expected format."""
        stats = {}
        for scenario_name, files in self.raw_data.items():
            stats[scenario_name] = files  # Placeholder for specific transformations

        if self.alias_dict:
            stats_with_alias = {}
            for scenario_name, df in stats.items():
                alias_name = self.alias_dict.get(scenario_name, scenario_name)
                stats_with_alias[alias_name] = df
            stats = stats_with_alias

        if self.new_order:
            stats = {key: stats[key] for key in self.new_order if key in stats}

        return stats

    def _get_carriers(self) -> pd.DataFrame:
        """Get carriers information from the network."""
        carriers = self.network.carriers
        carriers["legend_name"] = carriers.nice_name
        carriers.loc["DC", "legend_name"] = "Transmission"
        carriers.loc["DC", "color"] = "#cf1dab"
        carriers.loc["battery", "legend_name"] = "Existing BESS"
        carriers = carriers.set_index("nice_name")
        carriers = carriers.sort_values(by="co2_emissions", ascending=False)
        carriers.to_csv(self.figures_path / "carriers.csv")
        return carriers

    def prepare_combined_dataframe(
        self,
        variable: str,
        include_link: bool = False,
        as_pct: bool = False,
        variable_units: str | None = None,
    ) -> pd.DataFrame:
        """
        Prepare a combined DataFrame for all scenarios and horizons.

        Parameters
        ----------
        variable : str
            The variable to extract (e.g., "Optimal Capacity")
        include_link : bool, default False
            Whether to include Link components
        as_pct : bool, default False
            Whether to convert values to percentages
        variable_units : str, optional
            Units for the variable (affects scaling)

        Returns
        -------
        pd.DataFrame
            Combined DataFrame with all scenario data
        """
        factor_units = {"GW": 1e3, "GWh": 1e3, "%": 1}.get(variable_units, 1e9)

        data = []
        for scenario, df in self.processed_data.items():
            df = df["statistics"].fillna(0)

            tech_filter = ["Generator", "StorageUnit", "Link"]

            df = df.loc[df.index.get_level_values(0).isin(tech_filter), variable]
            df.index = df.index.get_level_values(1)
            df = df.reindex(self.carriers.index).dropna()

            if "capacity" in variable.lower():
                df = df.loc[~(df.index == "Load shedding")]

            if as_pct:
                df = ((df / df.sum()) * 100).round(2)

            for horizon in df.columns:
                df_horizon = df[horizon] / factor_units
                df_horizon = df_horizon.reset_index()
                df_horizon["Scenario"] = scenario
                df_horizon["horizon"] = horizon
                data.append(df_horizon.rename(columns={horizon: "statistics"}))

        combined_df = pd.concat(data, ignore_index=True)
        combined_df["scenario_name"] = combined_df["Scenario"].apply(
            lambda x: x.split("_")[0],
        )
        combined_df["trans_expansion"] = combined_df["Scenario"].apply(
            lambda x: x.split("_")[1],
        )
        combined_df.to_csv(self.figures_path / f"{variable}_comparison.csv")
        return combined_df


class ScenarioPlotter:
    """Class for creating scenario comparison plots."""

    def __init__(self, data_getter: ScenarioDataGetter):
        """
        Initialize with a data getter instance.

        Parameters
        ----------
        data_getter : ScenarioDataGetter
            Instance with loaded scenario data
        """
        self.data_getter = data_getter
        self.figures_path = data_getter.figures_path
        self.carriers = data_getter.carriers
        self.colors = data_getter.carriers["color"]
        self.network = data_getter.network
        self.processed_data = data_getter.processed_data
        self.reference_scenario = data_getter.reference_scenario

    def plot_scenario_comparison(
        self,
        combined_df: pd.DataFrame,
        variable: str,
        variable_units: str,
        title: str,
        include_link: bool = False,
    ) -> None:
        """
        Plot scenario comparison as horizontal bars.

        Parameters
        ----------
        combined_df : pd.DataFrame
            Combined DataFrame with all scenario data
        variable : str
            Variable name for file naming
        variable_units : str
            Units for the variable (for axis label)
        title : str
            Plot title
        include_link : bool, default False
            Whether to include Link components
        """
        planning_horizons = combined_df["horizon"].unique()
        scenarios = combined_df["Scenario"].unique()
        if not include_link:
            combined_df = combined_df[~combined_df["nice_name"].isin(["Link", "Ac"])]

        fig, axes = plt.subplots(
            nrows=len(planning_horizons),
            ncols=1,
            figsize=(8, 1.2 * len(planning_horizons) + 0.2 * len(scenarios)),
            sharex=True,
        )
        axes = np.atleast_1d(axes)  # Ensure axes is iterable for single horizon

        for ax, horizon in zip(axes, planning_horizons):
            horizon_df = combined_df[combined_df["horizon"] == horizon]
            y_positions = np.arange(len(scenarios))

            for j, scenario in enumerate(scenarios):
                scenario_df = horizon_df[horizon_df["Scenario"] == scenario]
                bottoms = np.zeros(len(y_positions))
                for tech in scenario_df["nice_name"].unique():
                    values = scenario_df[scenario_df["nice_name"] == tech]["statistics"].values[0]
                    ax.barh(
                        y_positions[j],
                        values,
                        left=bottoms[j],
                        color=self.colors[tech],
                        label=tech if j == 0 else "",
                    )
                    bottoms[j] += values

            ax.text(
                1.01,
                0.5,
                horizon,
                transform=ax.transAxes,
                va="center",
                rotation="vertical",
            )
            ax.set_yticks(y_positions)
            ax.set_yticklabels(scenarios)
            ax.grid(True, axis="x", linestyle="--", alpha=0.5)

        plt.xlabel(f"{variable} [{variable_units}]")
        plt.subplots_adjust(hspace=0.5)

        # Add legend
        carriers_plotted = self.carriers.loc[self.carriers.index.intersection(combined_df["nice_name"].unique())]
        legend_handles = [plt.Rectangle((0, 0), 1, 1, color=self.colors[tech]) for tech in carriers_plotted.index]
        fig.legend(
            handles=legend_handles,
            labels=carriers_plotted.legend_name.tolist(),
            loc="lower center",
            bbox_to_anchor=(0.5, -0.4),
            ncol=4,
            title="Technologies",
        )

        fig.suptitle(title, fontsize=12, fontweight="bold")
        plt.tight_layout()
        plt.savefig(
            self.figures_path / f"{variable}_comparison.png",
            dpi=300,
            bbox_inches="tight",
        )

        if self.reference_scenario:
            self._plot_reference_comparison(
                combined_df,
                self.reference_scenario,
                variable,
                horizon,
            )

    def _plot_reference_comparison(
        self,
        combined_df: pd.DataFrame,
        reference_scenario: str,
        variable: str,
        horizon: str,
    ) -> None:
        """
        Plot comparison against a reference scenario.

        Parameters
        ----------
        combined_df : pd.DataFrame
            Combined DataFrame with all scenario data
        reference_scenario : str
            Name of the reference scenario
        variable : str
            Variable name for file naming
        horizon : str
            Planning horizon to use for comparison
        """
        combined_df = combined_df.set_index("Scenario")
        combined_df = combined_df.loc[combined_df.horizon == horizon]
        ref = combined_df.loc[reference_scenario].set_index("nice_name")
        combined_df = combined_df.reset_index().set_index("nice_name")
        for carrier in combined_df.index.unique():
            combined_df.loc[carrier, "statistics"] = (
                (combined_df.loc[carrier, "statistics"] - ref.loc[carrier, "statistics"]) / ref.statistics.sum() * 100
            )
        combined_df = combined_df.reset_index().set_index("Scenario")
        stacked_data = combined_df.reset_index().pivot(
            index="Scenario",
            columns="nice_name",
            values="statistics",
        )
        stacked_data.plot(
            kind="bar",
            stacked=True,
            figsize=(10, 7),
            color=[self.colors[tech] for tech in stacked_data.columns],
            legend=False,
        )
        plt.ylabel("∆ Capacity[%]")
        plt.savefig(
            self.figures_path / f"{variable}_pct_comparison.png",
            dpi=300,
            bbox_inches="tight",
        )
        combined_df.to_csv(self.figures_path / f"{variable}_pct_comparison.csv")

    def scenario_comparison(
        self,
        variable: str,
        variable_units: str,
        title: str,
        include_link: bool = False,
        as_pct: bool = False,
    ) -> pd.DataFrame:
        """
        Original scenario comparison function, maintained for compatibility.

        Parameters
        ----------
        variable : str
            The variable to extract (e.g., "Optimal Capacity")
        variable_units : str
            Units for the variable (for axis label)
        title : str
            Plot title
        include_link : bool, default False
            Whether to include Link components
        as_pct : bool, default False
            Whether to convert values to percentages

        Returns
        -------
        pd.DataFrame
            Combined DataFrame with all scenario data
        """
        combined_df = pd.DataFrame(
            columns=["Scenario", "horizon", "nice_name", "statistics"],
            index=[],
        )
        stats = self.processed_data
        planning_horizons = stats[next(iter(stats.keys()))]["statistics"][variable].columns

        if variable_units in ["GW", "GWh"]:
            factor_units = 1e3
        elif variable_units in ["%"]:
            factor_units = 1
        else:
            factor_units = 1e9

        fig, axes = plt.subplots(
            nrows=len(planning_horizons),
            ncols=1,
            figsize=(8, 1.5 * len(planning_horizons) + 0.2 * len(stats)),
            sharex=True,
        )
        if len(planning_horizons) == 1:
            axes = [axes]

        for ax, horizon in zip(axes, planning_horizons):
            y_positions = np.arange(len(stats))
            for j, (scenario, df) in enumerate(stats.items()):
                df = df["statistics"].fillna(0)
                bottoms = np.zeros(len(y_positions))
                if include_link:
                    df = df.loc[
                        df.index.get_level_values(0).isin(
                            ["Generator", "StorageUnit", "Link"],
                        ),
                        variable,
                    ]
                    df = df.loc[~(df.index.get_level_values(1) == "Ac")]
                else:
                    df = df.loc[
                        df.index.get_level_values(0).isin(["Generator", "StorageUnit"]),
                        variable,
                    ]
                df.index = df.index.get_level_values(1)
                df = df.reindex(self.carriers.index).dropna()
                if as_pct:
                    df = ((df / df.sum()) * 100).round(2)
                for i, tech in enumerate(df.index.unique()):
                    values = df.loc[tech, horizon] / factor_units
                    ax.barh(
                        y_positions[j],
                        values,
                        left=bottoms[j],
                        color=self.colors[tech],
                        label=tech if j == 0 else "",
                    )
                    bottoms[j] += values

                df_copy = df.copy()
                df_copy[["Scenario", "horizon"]] = scenario, horizon
                df_copy = df_copy.reset_index()
                df_copy = df_copy.rename(columns={horizon: "statistics"})
                combined_df = pd.concat(
                    [combined_df, df_copy[["Scenario", "nice_name", "statistics", "horizon"]]],
                )

            ax.text(
                1.01,
                0.5,
                f"{horizon}",
                transform=ax.transAxes,
                va="center",
                rotation="vertical",
            )
            ax.set_yticks(y_positions)
            ax.set_yticklabels(stats.keys())
            ax.grid(True, axis="x", linestyle="--", alpha=0.5)

        plt.xlabel(f"{variable} [{variable_units}]")
        plt.subplots_adjust(hspace=0)
        carriers_plotted = self.carriers.loc[self.carriers.index.intersection(df.index.unique())]
        legend_handles = [plt.Rectangle((0, 0), 1, 1, color=self.colors[tech]) for tech in carriers_plotted.index]
        fig.legend(
            handles=legend_handles,
            labels=carriers_plotted.legend_name.tolist(),
            loc="lower center",
            bbox_to_anchor=(0.5, -0.4),
            ncol=4,
            title="Technologies",
        )
        fig.suptitle(title, fontsize=12, fontweight="bold")
        plt.tight_layout()
        plt.savefig(
            self.figures_path / f"{variable}_comparison.png",
            dpi=300,
            bbox_inches="tight",
        )

        combined_df["scenario_name"] = combined_df["Scenario"].apply(
            lambda x: x.split("_")[0],
        )
        combined_df["trans_expansion"] = combined_df["Scenario"].apply(
            lambda x: x.split("_")[1],
        )
        combined_df.to_csv(self.figures_path / f"{variable}_comparison.csv")

        if self.reference_scenario:
            # last_horizon = planning_horizons[-1]
            combined_df = combined_df.set_index("Scenario")
            combined_df = combined_df.query("horizon == @last_horizon").drop(
                columns="horizon",
            )  # only plot last horizon
            ref = combined_df.loc[self.reference_scenario].set_index("nice_name")
            combined_df = combined_df.reset_index().set_index("nice_name")
            for scenario in combined_df.index.unique():
                combined_df.loc[scenario, "statistics"] = (
                    (combined_df.loc[scenario, "statistics"] - ref.loc[scenario, "statistics"])
                    / ref.statistics.sum()
                    * 100
                )
            combined_df = combined_df.reset_index().set_index("Scenario")
            stacked_data = combined_df.reset_index().pivot(
                index="Scenario",
                columns="nice_name",
                values="statistics",
            )
            stacked_data.plot(
                kind="bar",
                stacked=True,
                figsize=(10, 7),
                color=[self.colors[tech] for tech in stacked_data.columns],
                legend=False,
            )
            plt.ylabel("∆ Capacity[%]")
            plt.savefig(
                self.figures_path / f"{variable}_pct_comparison.png",
                dpi=300,
                bbox_inches="tight",
            )
            combined_df.to_csv(self.figures_path / f"{variable}_pct_comparison.csv")
        return combined_df

    def plot_cost_comparison(
        self,
        variable: str,
        variable_units: str,
        title: str,
    ) -> None:
        """
        Plot cost comparison between scenarios.

        Parameters
        ----------
        variable : str
            Variable name for file naming
        variable_units : str
            Units for the variable (for axis label)
        title : str
            Plot title
        """
        stats = self.processed_data
        n = self.network
        combined_df = pd.DataFrame(columns=["Scenario", "statistics"], index=[])

        for j, (scenario, stat) in enumerate(stats.items()):
            stat = stat["statistics"]
            combined_df = pd.concat(
                [
                    combined_df,
                    pd.DataFrame(
                        {
                            "Scenario": scenario,
                            "statistics": (
                                (stat["Capital Expenditure"].sum() + stat["Operational Expenditure"].sum())
                                * n.investment_period_weightings.objective.values
                            ).sum()
                            / 1e9,
                        },
                        index=[j],
                    ),
                ],
            )

        combined_df.plot(
            kind="bar",
            x="Scenario",
            y="statistics",
            title="Total System Costs",
            legend=False,
        )
        plt.ylabel("Annualized System Costs [B$]")
        plt.savefig(
            self.figures_path / f"{variable}_comparison.png",
            dpi=300,
            bbox_inches="tight",
        )

        if self.reference_scenario:
            combined_df = combined_df.set_index("Scenario")
            ref = combined_df.loc[self.reference_scenario]
            pct_df = (combined_df - ref) / combined_df * 100
            pct_df.plot(
                kind="bar",
                y="statistics",
                title="Total System Costs",
                legend=False,
            )
            plt.ylabel("∆ Annualized System Costs [%]")
            plt.savefig(
                self.figures_path / f"{variable}_pct_comparison.png",
                dpi=300,
                bbox_inches="tight",
            )
            pct_df.to_csv(self.figures_path / f"{variable}_pct_comparison.csv")
            combined_df.to_csv(self.figures_path / f"{variable}_comparison.csv")


def main():
    """Main execution function."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run scenario comparison script.")
    parser.add_argument(
        "yaml_name",
        type=str,
        help="Name of the YAML configuration file.",
    )
    parser.add_argument(
        "--force-regenerate",
        action="store_true",
        help="Force regeneration of statistics files even if they exist in cache",
    )
    args = parser.parse_args()

    yaml_path = Path.cwd() / args.yaml_name  # Path to the YAML file

    # Initialize data getter and plotter
    data_getter = ScenarioDataGetter(yaml_path, force_regenerate=args.force_regenerate)
    plotter = ScenarioPlotter(data_getter)

    # Generate plots for Optimal Capacity
    variable = "Optimal Capacity"
    variable_units = "GW"
    title = "Capacity Comparison"
    combined_df = data_getter.prepare_combined_dataframe(
        variable,
        as_pct=False,
        variable_units=variable_units,
    )
    plotter.plot_scenario_comparison(
        combined_df,
        variable,
        variable_units,
        title,
    )

    # Generate plots for Supply
    variable = "Supply"
    variable_units = "%"
    title = "Supply Comparison"
    plotter.scenario_comparison(
        variable,
        variable_units,
        title,
        as_pct=True,
    )

    # Generate plots for CAPEX
    variable = "Capital Expenditure"
    variable_units = "$B"
    title = "CAPEX Comparison"
    plotter.scenario_comparison(
        variable,
        variable_units,
        title,
        as_pct=False,
    )

    # Generate plots for OPEX
    variable = "Operational Expenditure"
    variable_units = "$B"
    title = "OPEX Comparison"
    plotter.scenario_comparison(
        variable,
        variable_units,
        title,
        as_pct=False,
    )

    # Generate plots for System Costs
    variable = "System Costs"
    variable_units = "$B"
    title = "Scenario Comparison"
    plotter.plot_cost_comparison(
        variable,
        variable_units,
        title,
    )


if __name__ == "__main__":
    main()
