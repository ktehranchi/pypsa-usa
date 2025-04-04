import logging  # noqa: D100

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging
from add_electricity import sanitize_carriers

logger = logging.getLogger(__name__)


def plot_stat_vs_clusters(stats_dict, stat_name, output_path=None, title=None, ylabel=None, carriers_to_plot=None):
    """
    Plot how a specific statistic changes with the number of clusters,
    with separate lines for each energy carrier.

    Parameters
    ----------
    stats_dict : dict
        Dictionary with cluster numbers as keys and lists of statistics as values
    stat_name : str
        Name of the statistic to plot (must be a column in the stats dataframes)
    output_path : str, optional
        Path to save the figure
    title : str, optional
        Plot title
    ylabel : str, optional
        Y-axis label
    carriers_to_plot : list, optional
        Specific carriers to include in the plot. If None, all carriers are plotted.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    # Set larger font sizes for all elements
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.titlesize": 18,
            "axes.labelsize": 16,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 12,
            "figure.titlesize": 20,
        },
    )

    clusters = sorted(stats_dict.keys())

    # Initialize a dictionary to store values for each carrier
    carrier_values = {}

    # First, discover all carriers across all networks
    all_carriers = set()
    for n_clusters in clusters:
        for network_stats in stats_dict[n_clusters]:
            stats_df = network_stats["stats"]
            if stat_name in stats_df.columns:
                # Get carriers from the first level of MultiIndex
                for carrier_tuple in stats_df.index:
                    all_carriers.add(carrier_tuple)

    # Filter carriers if specified
    if carriers_to_plot:
        carriers_to_include = [c for c in all_carriers if any(carrier in c for carrier in carriers_to_plot)]
    else:
        carriers_to_include = all_carriers

    # Initialize carrier_values with empty lists for each carrier
    for carrier in carriers_to_include:
        carrier_values[carrier] = [np.nan] * len(clusters)

    # Fill in values for each cluster
    for i, n_clusters in enumerate(clusters):
        cluster_carrier_values = {carrier: [] for carrier in carriers_to_include}

        for network_stats in stats_dict[n_clusters]:
            stats_df = network_stats["stats"]
            if stat_name in stats_df.columns:
                for carrier in carriers_to_include:
                    if carrier in stats_df.index:
                        value = stats_df.loc[carrier, stat_name]
                        if isinstance(value, pd.Series):
                            value = value.iloc[0]  # Take the first value if it's a Series
                        cluster_carrier_values[carrier].append(value)

        # Average values for this cluster
        for carrier in carriers_to_include:
            if cluster_carrier_values[carrier]:
                # breakpoint()
                carrier_values[carrier][i] = np.nanmean(cluster_carrier_values[carrier])

    # Create plot
    fig, ax = plt.subplots(figsize=(14, 10))  # Increased figure size

    # Colors for different carriers
    colors = plt.cm.tab20(np.linspace(0, 1, len(carriers_to_include)))

    # Plot each carrier
    for i, (carrier, values) in enumerate(carrier_values.items()):
        if any(not np.isnan(v) for v in values):  # Only plot if there's at least one non-NaN value
            label = f"{carrier[0]} - {carrier[1]}" if isinstance(carrier, tuple) else str(carrier)
            ax.plot(
                clusters,
                values,
                "o-",
                linewidth=2.5,
                markersize=8,  # Increased marker and line size
                color=colors[i],
                label=label,
            )

    # Labels and formatting
    ax.set_xlabel("Number of Clusters", fontsize=16, fontweight="bold")
    ax.set_ylabel(ylabel if ylabel else stat_name, fontsize=16, fontweight="bold")
    ax.set_title(
        title if title else f"{stat_name} vs Number of Clusters",
        fontsize=18,
        fontweight="bold",
        pad=20,
    )  # Add padding to title
    ax.grid(True, alpha=0.3, linestyle="--", linewidth=1)  # Improved grid

    # Set x-axis to only show integer values corresponding to cluster counts
    ax.set_xticks(clusters)
    ax.tick_params(axis="both", which="major", labelsize=14, width=2, length=6)  # Enhanced tick marks

    # Add legend with increased font size and improved placement
    if len(carriers_to_include) > 1:
        legend = ax.legend(
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            borderaxespad=0.0,
            fontsize=12,
            frameon=True,
        )
        legend.get_frame().set_edgecolor("black")
        legend.get_frame().set_linewidth(1)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Plot saved to {output_path}")

    # Reset matplotlib parameters to defaults for future plots
    plt.rcdefaults()

    return fig


def plot_cost_vs_clusters(
    stats_dict,
    output_path=None,
    title=None,
    ylabel=None,
):
    """
    Plot the relative cost of each network compared to the lowest cost network.

    Parameters
    ----------
    stats_dict : dict
        Dictionary with cluster numbers as keys and lists of statistics as values
    output_path : str, optional
        Path to save the figure
    title : str, optional
        Plot title
    ylabel : str, optional
        Y-axis label

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    # Set larger font sizes for all elements
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.titlesize": 18,
            "axes.labelsize": 16,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 12,
            "figure.titlesize": 20,
        },
    )

    # Get the system cost for the highest cluster size network
    min_cost = float("inf")
    for n_clusters, networks in stats_dict.items():
        for network_stats in networks:
            stats_df = network_stats["stats"]
            capex = stats_df["Capital Expenditure"].sum()
            opex = stats_df["Operational Expenditure"].sum()
            cost = capex + opex
            if cost.values[0] < min_cost:
                min_cost = cost.values[0]

    # Calculate relative costs for each network
    relative_costs = {}
    for n_clusters, networks in stats_dict.items():
        for network_stats in networks:
            stats_df = network_stats["stats"]
            capex = stats_df["Capital Expenditure"].sum().values[0]
            opex = stats_df["Operational Expenditure"].sum().values[0]
            cost = capex + opex
            relative_costs[n_clusters] = (cost - min_cost) / min_cost * 100

    # Create plot
    fig, ax = plt.subplots(figsize=(14, 10))  # Increased figure size

    # Sort the cluster numbers for proper line plotting
    sorted_clusters = sorted(relative_costs.keys())
    sorted_costs = [relative_costs[cluster] for cluster in sorted_clusters]

    # Plot relative costs as a line
    ax.plot(
        sorted_clusters,
        sorted_costs,
        "o-",
        linewidth=2.5,
        markersize=8,
        color="royalblue",
        label="Relative Cost",
    )

    # Labels and formatting
    ax.set_xlabel("Number of Clusters", fontsize=16, fontweight="bold")
    ax.set_ylabel(ylabel if ylabel else "Relative Cost [% of lowest cost]", fontsize=16, fontweight="bold")
    ax.set_title(title if title else "Relative Cost vs Number of Clusters", fontsize=18, fontweight="bold", pad=20)

    # Set x-axis to only show integer values corresponding to cluster counts
    ax.set_xticks(sorted_clusters)

    # Add grid lines and improve tick marks
    ax.grid(True, alpha=0.3, linestyle="--", linewidth=1)
    ax.tick_params(axis="both", which="major", labelsize=14, width=2, length=6)

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Plot saved to {output_path}")

    # Reset matplotlib parameters to defaults for future plots
    plt.rcdefaults()

    return fig


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_statistics",
            interconnect="texas",
            clusters=7,
            ll="v1.00",
            opts="REM-400SEG",
            sector="E",
        )
    configure_logging(snakemake)

    # extract shared plotting files
    networks = snakemake.input.networks
    # read in the networks into dict based on the number of clusters
    networks_dict = {}
    stats_dict = {}  # Dictionary to store statistics for each network

    for network in networks:
        n = pypsa.Network(network)
        # get the number of clusters from the filename
        clusters_str = network.split("_")[-8]
        clusters = int("".join(c for c in clusters_str if c.isdigit()))
        if clusters not in networks_dict:
            networks_dict[clusters] = []
            stats_dict[clusters] = []  # Initialize list for storing statistics
        networks_dict[clusters].append(network)

        # Extract and store statistics for this network
        network_stats = n.statistics().round(3)
        stats_dict[clusters].append({"network": network, "stats": network_stats})

    sanitize_carriers(n, snakemake.config)

    # carriers to plot
    carriers = (
        snakemake.params.electricity["conventional_carriers"]
        + snakemake.params.electricity["renewable_carriers"]
        + snakemake.params.electricity["extendable_carriers"]["Generator"]
        + snakemake.params.electricity["extendable_carriers"]["StorageUnit"]
        + snakemake.params.electricity["extendable_carriers"]["Store"]
        + snakemake.params.electricity["extendable_carriers"]["Link"]
    )
    carriers = list(set(carriers))  # remove any duplicates

    # Create Statistics Plots
    # Example usage of the new plotting function with carriers
    plot_stat_vs_clusters(
        stats_dict,
        "Optimal Capacity",
        output_path=snakemake.output["capacity_comparison.pdf"],
        title="Optimal Capacity vs Number of Clusters",
        ylabel="Optimal Capacity (MW)",
        # Optionally filter specific carriers
        # carriers_to_plot=["Onshore Wind", "Solar", "Battery Storage"]
    )

    plot_stat_vs_clusters(
        stats_dict,
        "Supply",
        output_path=snakemake.output["production_comparison.pdf"],
        title="Supply vs Number of Clusters",
        ylabel="Supply (MWh)",
    )

    plot_stat_vs_clusters(
        stats_dict,
        "Supply",
        output_path=snakemake.output["production_comparison.pdf"],
        title="Supply vs Number of Clusters",
        ylabel="Supply (MWh)",
    )

    plot_cost_vs_clusters(
        stats_dict,
        output_path=snakemake.output["objective_comparison.pdf"],
        title="Relative Cost vs Number of Clusters",
        ylabel="Difference in Cost [% of lowest cost]",
    )
