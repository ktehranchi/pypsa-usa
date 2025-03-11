"""Script used to compare outputs from multiple snakemake scenarios."""

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import yaml
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgba

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


class ScenarioDataGetter:
    """Class for loading and processing scenario data."""

    def __init__(
        self,
        yaml_path: str | Path,
        force_regenerate: bool = False,
        skip_plots: bool = False,
    ):
        """
        Initialize with configuration from YAML file.

        Parameters
        ----------
        yaml_path : str or Path
            Path to the YAML configuration file
        force_regenerate : bool, default False
            Whether to force regeneration of statistics files even if they exist
        skip_plots : bool, default False
            Whether to skip generating plots and only process data
        """
        self.config = self._load_yaml_config(yaml_path)
        self.scenarios = self.config["scenarios"]
        self.alias_dict = self.config.get("alias_dict", None)
        self.new_order = self.config.get("new_order", None)
        self.reference_scenario = self.config.get("reference_scenario", None)
        self.force_regenerate = force_regenerate
        self.skip_plots = skip_plots
        # Set output path
        output_folder = self.config.get("output_folder_name", "scenario_comparison")
        self.figures_path = Path.cwd() / f"results/{output_folder}"
        self.figures_path.mkdir(exist_ok=True, parents=True)

        # Create cache directory
        self.cache_path = self.figures_path / "cached_statistics"
        self.cache_path.mkdir(exist_ok=True)

        self.plot_data_path = self.figures_path / "plot_data"
        self.plot_data_path.mkdir(exist_ok=True)

        # Create metadata file to store info about runs
        self.metadata_file = self.cache_path / "metadata.json"
        self.metadata = self._load_or_create_metadata()

        # Set color scheme from config or use default
        self.color_scheme = self.config.get("color_scheme", "default")

        # Load data
        try:
            self.raw_data = self._load_scenario_data()
            self.processed_data = self._process_data()
            # Load network
            self.network = pypsa.Network(self.config["network"]["path"])
            self.carriers = self._get_carriers()
        except Exception as e:
            logger.error(f"Error initializing data: {e}", exc_info=True)
            raise

    def _load_or_create_metadata(self) -> dict:
        """Load existing metadata or create new if doesn't exist."""
        if self.metadata_file.exists():
            try:
                with open(self.metadata_file) as f:
                    return json.load(f)
            except json.JSONDecodeError:
                logger.warning("Corrupted metadata file. Creating new.")
                return {"scenarios": {}, "last_run": None}
        else:
            return {"scenarios": {}, "last_run": None}

    def _update_metadata(self, scenario_name: str, network_path: str) -> None:
        """Update metadata with info about generated statistics."""
        self.metadata["scenarios"][scenario_name] = {
            "network_path": str(network_path),
            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "statistics_path": str(self.cache_path / f"statistics_{scenario_name}.csv"),
        }
        self.metadata["last_run"] = time.strftime("%Y-%m-%d %H:%M:%S")

        with open(self.metadata_file, "w") as f:
            json.dump(self.metadata, f, indent=2)

    def _load_yaml_config(self, yaml_path: str | Path) -> dict:
        """Load configuration from YAML file."""
        with open(yaml_path) as file:
            return yaml.safe_load(file)

    def _load_scenario_data(self) -> dict[str, dict[str, pd.DataFrame]]:
        """
        Load CSV data for all scenarios, generating statistics first if needed.

        If cached statistics exist, they are loaded instead of regenerating.
        """
        data = {}
        # Check for file existence before starting
        missing_networks = []
        for scenario in self.scenarios:
            network_path = Path(scenario["path"])
            if not network_path.exists():
                missing_networks.append((scenario["name"], str(network_path)))
        if missing_networks:
            error_msg = "The following network files could not be found:\n"
            for name, path in missing_networks:
                error_msg += f"  - {name}: {path}\n"
            logger.error(error_msg)
            if not self.config.get("ignore_missing_files", False):
                raise FileNotFoundError(error_msg)

        for scenario in self.scenarios:
            scenario_name = scenario["name"]
            network_path = Path(scenario["path"])

            # Skip if network doesn't exist and we're ignoring missing files
            if not network_path.exists() and self.config.get("ignore_missing_files", False):
                logger.warning(f"Skipping missing network for {scenario_name}")
                continue
            # Check if cached statistics exist
            cached_file = self.cache_path / f"statistics_{scenario_name}.csv"

            # Generate statistics if forced or if cache doesn't exist
            if self.force_regenerate or not cached_file.exists():
                logger.info(f"Generating statistics for {scenario_name}...")
                start_time = time.time()

                try:
                    logger.info(f"Loading network from {network_path}")
                    n = pypsa.Network(network_path)

                    # Generate statistics
                    stats = self._generate_statistics(n)

                    # Save to cache
                    stats.to_csv(cached_file)

                    # Update metadata
                    self._update_metadata(scenario_name, str(network_path))

                    elapsed = time.time() - start_time
                    logger.info(f"Statistics for {scenario_name} generated in {elapsed:.1f} seconds and saved to cache")

                    data[scenario_name] = {"statistics": stats}

                except Exception as e:
                    logger.error(f"Error processing {scenario_name}: {e}", exc_info=True)
                    if not self.config.get("continue_on_error", False):
                        raise
            else:
                logger.info(f"Loading cached statistics for {scenario_name}...")
                try:
                    data[scenario_name] = {
                        "statistics": pd.read_csv(cached_file, index_col=[0, 1], header=[0, 1]),
                    }
                except Exception as e:
                    logger.error(f"Error loading cached statistics for {scenario_name}: {e}")
                    if self.config.get("continue_on_error", False):
                        continue
                    raise
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
        # Basic approach: use the built-in statistics method if available
        try:
            return n.statistics()
        except AttributeError:
            logger.warning("Network object does not have statistics method. Using custom implementation.")
            return self._custom_statistics(n)

    def _custom_statistics(self, n: pypsa.Network) -> pd.DataFrame:
        """
        Custom implementation of statistics for older PyPSA versions.

        Parameters
        ----------
        n : pypsa.Network
            Network to analyze

        Returns
        -------
        pd.DataFrame
            DataFrame containing statistics
        """
        # Create a dictionary to store statistics DataFrames
        stats = {}

        # Get investment periods or just a single period
        if hasattr(n, "investment_periods") and len(n.investment_periods) > 0:
            horizons = n.investment_periods
        else:
            # If no investment periods, use the year from the first snapshot
            horizons = [pd.Timestamp(n.snapshots[0]).year]

        # Calculate optimal capacity
        capacity = pd.DataFrame(
            index=pd.MultiIndex.from_product(
                [
                    ["Generator", "StorageUnit", "Link"],
                    n.carriers.index,
                ],
            ),
            columns=horizons,
        )

        # Generators
        for c in n.generators.carrier.unique():
            for p in horizons:
                capacity.loc[("Generator", c), p] = n.generators[n.generators.carrier == c].p_nom_opt.sum()

        # Storage units
        if not n.storage_units.empty:
            for c in n.storage_units.carrier.unique():
                for p in horizons:
                    capacity.loc[("StorageUnit", c), p] = n.storage_units[n.storage_units.carrier == c].p_nom_opt.sum()

        # Links
        if not n.links.empty:
            for c in n.links.carrier.unique():
                for p in horizons:
                    capacity.loc[("Link", c), p] = n.links[n.links.carrier == c].p_nom_opt.sum()

        stats["Optimal Capacity"] = capacity.fillna(0)

        # Implement other statistics like Supply, CAPEX, OPEX as needed
        # This is a simplified placeholder

        # Create a MultiIndex DataFrame with all statistics
        return pd.concat(stats, axis=1)

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
        carriers = self.network.carriers.copy()
        # Use custom carrier colors from config if available
        if "carrier_colors" in self.config:
            for carrier, color in self.config["carrier_colors"].items():
                if carrier in carriers.index:
                    carriers.loc[carrier, "color"] = color

        # Set default names and colors for common carriers
        carriers["legend_name"] = carriers.get("nice_name", carriers.index)

        # Set specific properties for common carriers if they exist
        carrier_properties = {
            "DC": {"legend_name": "Transmission", "color": "#cf1dab"},
            "battery": {"legend_name": "Existing BESS"},
            "solar": {"legend_name": "Solar PV"},
            "offwind": {"legend_name": "Offshore Wind"},
            "onwind": {"legend_name": "Onshore Wind"},
            "hydrogen": {"legend_name": "Hydrogen Storage"},
        }
        for carrier, properties in carrier_properties.items():
            if carrier in carriers.index:
                for key, value in properties.items():
                    carriers.loc[carrier, key] = value

        # Use index-based access for all values to prevent SettingWithCopyWarning
        carriers = carriers.set_index("nice_name") if "nice_name" in carriers.columns else carriers

        # Sort by emissions for better visualization
        try:
            carriers = carriers.sort_values(by="co2_emissions", ascending=False)
        except KeyError:
            # No co2_emissions column, just leave as is
            pass

        # Save for reference
        carriers.to_csv(self.plot_data_path / "carriers.csv")

        return carriers

    def prepare_combined_dataframe(
        self,
        variable: str,
        include_link: bool = False,
        as_pct: bool = False,
        variable_units: str | None = None,
        filters: dict[str, list[str]] | None = None,
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
        filters : Dict[str, List[str]], optional
            Filters to apply, e.g., {"carrier": ["solar", "wind"]}

        Returns
        -------
        pd.DataFrame
            Combined DataFrame with all scenario data
        """
        factor_units = {"GW": 1e3, "GWh": 1e3, "%": 1}.get(variable_units, 1e9)
        filters = filters or {}

        data = []
        for scenario, df_dict in self.processed_data.items():
            try:
                df = df_dict["statistics"].fillna(0)

                tech_filter = ["Generator", "StorageUnit", "Link"]
                if not include_link:
                    tech_filter = tech_filter[:2]  # Exclude Link

                # Get data for the specified variable and filter by component type
                try:
                    df = df.loc[df.index.get_level_values(0).isin(tech_filter), variable]
                except KeyError:
                    logger.warning(f"Variable '{variable}' not found in statistics for {scenario}")
                    continue

                # Extract carrier from index and set as new index
                df.index = df.index.get_level_values(1)

                # Apply any additional filters
                if "carrier" in filters:
                    df = df.loc[df.index.isin(filters["carrier"])]

                # Align with carriers
                df = df.reindex(self.carriers.index).dropna()

                # Check if we need to exclude certain carriers (like "Load shedding")
                if "capacity" in variable.lower():
                    exclude_carriers = self.config.get("exclude_from_capacity", ["Load shedding"])
                    df = df.loc[~df.index.isin(exclude_carriers)]

                if as_pct:
                    # Avoid division by zero
                    total = df.sum()
                    if (total > 0).all():
                        df = ((df / total) * 100).round(2)
                    else:
                        logger.warning(
                            f"Zero sum detected in {scenario} for {variable}, skipping percentage conversion",
                        )

                # Process each planning horizon
                for horizon in df.columns:
                    df_horizon = df[horizon] / factor_units
                    df_horizon = df_horizon.reset_index()
                    df_horizon["Scenario"] = scenario
                    df_horizon["horizon"] = horizon
                    data.append(df_horizon.rename(columns={horizon: "statistics"}))

            except Exception as e:
                logger.error(f"Error processing {scenario} for {variable}: {e}", exc_info=True)
                if not self.config.get("continue_on_error", False):
                    raise

        if not data:
            logger.warning(f"No valid data found for variable '{variable}'")
            return pd.DataFrame()  # Return empty DataFrame instead of None

        # Combine all data
        combined_df = pd.concat(data, ignore_index=True)
        # Add scenario components
        try:
            combined_df["scenario_name"] = combined_df["Scenario"].apply(
                lambda x: x.split("_")[0],
            )
            combined_df["trans_expansion"] = combined_df["Scenario"].apply(
                lambda x: x.split("_")[1],
            )
        except (IndexError, KeyError) as e:
            logger.warning(f"Could not parse scenario components: {e}")

        # Save for reference
        combined_df.to_csv(self.plot_data_path / f"{variable}_comparison.csv")

        return combined_df

    def export_summary_table(self, variables: list[str]) -> None:
        """
        Export a summary table with key metrics for all scenarios.

        Parameters
        ----------
        variables : List[str]
            List of variables to include in summary
        """
        summary_data = []

        for scenario, df_dict in self.processed_data.items():
            row = {"Scenario": scenario}
            df = df_dict["statistics"]

            for variable in variables:
                if variable in df.columns:
                    # Get the sum of this variable across all components
                    try:
                        value = df[variable].sum().sum()
                        row[f"{variable} (Total)"] = value
                    except Exception:
                        row[f"{variable} (Total)"] = np.nan
                    # Get top contributors (carriers) for this variable
                    try:
                        # Group by carrier and sum
                        carriers = df[variable].groupby(level=1).sum().sum()
                        top_carriers = carriers.nlargest(3)
                        for i, (carrier, val) in enumerate(top_carriers.items(), 1):
                            row[f"{variable} Top {i}"] = carrier
                            row[f"{variable} Top {i} Value"] = val
                    except Exception:
                        pass

            summary_data.append(row)

        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(self.plot_data_path / "summary_table.csv", index=False)

        # Also export as formatted Excel
        try:
            summary_df.to_excel(self.plot_data_path / "summary_table.xlsx", index=False)
        except ImportError:
            logger.warning("openpyxl not installed, skipping Excel export")


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
        self.plot_data_path = data_getter.plot_data_path
        self.carriers = data_getter.carriers
        self.colors = data_getter.carriers["color"]
        self.network = data_getter.network
        self.processed_data = data_getter.processed_data
        self.reference_scenario = data_getter.reference_scenario
        self.config = data_getter.config

        # Set plot style
        self._set_plot_style()

    def _set_plot_style(self) -> None:
        """Set matplotlib plot style based on configuration."""
        plt.style.use(self.config.get("plot_style", "default"))

        # Custom figure size
        plt.rcParams["figure.figsize"] = self.config.get("figure_size", [10, 6])

        # Custom font sizes
        font_sizes = self.config.get("font_sizes", {})
        if font_sizes:
            for key, size in font_sizes.items():
                plt.rcParams[f"font.{key}"] = size

    def _get_color_with_alpha(self, color: str, alpha: float = 0.8) -> tuple[float, float, float, float]:
        """Get color with specified alpha."""
        rgba = to_rgba(color)
        return rgba[0], rgba[1], rgba[2], alpha

    def _add_plot_decorations(
        self,
        fig: plt.Figure,
        title: str,
        add_timestamp: bool = True,
        filename_suffix: str = "",
    ) -> None:
        """Add common decorations to plot."""
        # Add title
        if title:
            fig.suptitle(title, fontsize=12, fontweight="bold")
        # Add timestamp if requested
        if add_timestamp:
            timestamp = time.strftime("%Y-%m-%d %H:%M")
            fig.text(
                0.99,
                0.01,
                f"Generated: {timestamp}",
                ha="right",
                va="bottom",
                fontsize=8,
                alpha=0.7,
            )

        # Add watermark/logo if configured
        logo_path = self.config.get("logo_path")
        if logo_path:
            try:
                logo_img = plt.imread(logo_path)
                ax_logo = fig.add_axes([0.02, 0.02, 0.1, 0.1], frameon=False)
                ax_logo.imshow(logo_img)
                ax_logo.axis("off")
            except Exception as e:
                logger.warning(f"Could not add logo: {e}")

    def plot_scenario_comparison(
        self,
        combined_df: pd.DataFrame,
        variable: str,
        variable_units: str,
        title: str,
        include_link: bool = False,
        filename_suffix: str = "",
        plot_format: str = "png",
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
        filename_suffix : str, default ""
            Optional suffix for output filename
        plot_format : str, default "png"
            Format for output plot (png, pdf, svg)
        """
        if combined_df.empty:
            logger.warning(f"Empty dataframe for {variable}, skipping plot")
            return
        # Apply scenario filter if specified
        if "plot_scenarios" in self.config:
            filtered_scenarios = self.config["plot_scenarios"]
            combined_df = combined_df[combined_df["Scenario"].isin(filtered_scenarios)]

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
                if scenario_df.empty:
                    logger.warning(f"No data for scenario {scenario} in horizon {horizon}")
                    continue
                bottoms = np.zeros(len(y_positions))
                for tech in scenario_df["nice_name"].unique():
                    try:
                        values = scenario_df[scenario_df["nice_name"] == tech]["statistics"].values[0]
                        ax.barh(
                            y_positions[j],
                            values,
                            left=bottoms[j],
                            color=self.colors.get(tech, "#999999"),  # Default gray if color not found
                            label=tech if j == 0 else "",
                            alpha=0.8,  # Add some transparency
                            edgecolor="white" if self.config.get("bar_edgecolor", False) else None,
                            linewidth=0.5 if self.config.get("bar_edgecolor", False) else 0,
                        )
                        bottoms[j] += values
                    except (IndexError, KeyError) as e:
                        logger.warning(f"Error plotting {tech} for {scenario}: {e}")

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
        legend_handles = [
            plt.Rectangle((0, 0), 1, 1, color=self.colors.get(tech, "#999999")) for tech in carriers_plotted.index
        ]

        # Get legend name or use index if not available
        legend_labels = [
            carriers_plotted.get("legend_name", {}).get(tech, tech) if hasattr(carriers_plotted, "get") else tech
            for tech in carriers_plotted.index
        ]

        fig.legend(
            handles=legend_handles,
            labels=legend_labels,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.4),
            ncol=min(4, len(legend_handles)),  # Adjust columns based on number of items
            title="Technologies",
        )

        # Add decorations
        self._add_plot_decorations(fig, title)

        plt.tight_layout()
        output_file = self.figures_path / f"{variable}_comparison{filename_suffix}.{plot_format}"
        plt.savefig(
            output_file,
            dpi=300,
            bbox_inches="tight",
        )
        logger.info(f"Saved plot to {output_file}")

        # Close the figure to free memory
        plt.close(fig)

        if self.reference_scenario and self.reference_scenario in scenarios:
            self._plot_reference_comparison(
                combined_df,
                self.reference_scenario,
                variable,
                planning_horizons[-1],  # Use last horizon for reference comparison
                filename_suffix,
                plot_format,
            )

    def _plot_reference_comparison(
        self,
        combined_df: pd.DataFrame,
        reference_scenario: str,
        variable: str,
        horizon: str,
        filename_suffix: str = "",
        plot_format: str = "png",
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
        filename_suffix : str, default ""
            Optional suffix for output filename
        plot_format : str, default "png"
            Format for output plot (png, pdf, svg)
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
        output_file = self.figures_path / f"pct_{variable}_comparison{filename_suffix}.{plot_format}"
        plt.savefig(
            output_file,
            dpi=300,
            bbox_inches="tight",
        )
        logger.info(f"Saved plot to {output_file}")
        combined_df.to_csv(self.plot_data_path / f"pct_{variable}_comparison.csv")

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
        combined_df.to_csv(self.plot_data_path / f"{variable}_comparison.csv")

        if self.reference_scenario:
            last_horizon = planning_horizons[-1]  # noqa
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
                self.figures_path / f"pct_{variable}_comparison.png",
                dpi=300,
                bbox_inches="tight",
            )
            combined_df.to_csv(self.plot_data_path / f"pct_{variable}_comparison.csv")
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
                self.figures_path / f"pct_{variable}_comparison.png",
                dpi=300,
                bbox_inches="tight",
            )
            pct_df.to_csv(self.plot_data_path / f"pct_{variable}_comparison.csv")
            combined_df.to_csv(self.plot_data_path / f"{variable}_comparison.csv")


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
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip generating plots and only process data",
    )
    args = parser.parse_args()

    yaml_path = Path.cwd() / args.yaml_name  # Path to the YAML file

    # Initialize data getter and plotter
    data_getter = ScenarioDataGetter(yaml_path, force_regenerate=args.force_regenerate, skip_plots=args.skip_plots)
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
