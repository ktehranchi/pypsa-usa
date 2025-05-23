"""Adds policy data to the network.

Data on policy constraints to add to the network:
- emission limits
- technology capacity targets
- portfolio standards
- regional reserve margins
"""

import logging
from functools import partial

import pandas as pd
import pypsa
from _helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)


def load_csv_data(config, policy_type, **kwargs):
    """
    Generic function to load CSV data for different policy types.

    Parameters
    ----------
    config : dict
        Configuration dictionary containing file paths
    policy_type : str
        Type of policy data to load (e.g., 'technology_capacity_targets')
    **kwargs : dict
        Additional arguments to pass to pd.read_csv

    Returns
    -------
    pd.DataFrame
        Loaded policy data
    """
    return pd.read_csv(config["electricity"][policy_type], **kwargs)


def process_reeds_data(filepath, carriers, value_col):
    """
    Process RPS or CES REEDS data.

    Parameters
    ----------
    filepath : str
        Path to the REEDS data file
    carriers : list
        List of carrier types
    value_col : str
        Name of the value column

    Returns
    -------
    pd.DataFrame
        Processed REEDS data
    """
    reeds = pd.read_csv(filepath)

    # Handle both wide and long formats
    if "rps_all" not in reeds.columns:
        reeds = reeds.melt(
            id_vars="st",
            var_name="planning_horizon",
            value_name=value_col,
        )

    # Standardize column names
    reeds = reeds.rename(
        columns={"st": "region", "t": "planning_horizon", "rps_all": "pct"},
    )
    reeds["carrier"] = [", ".join(carriers)] * len(reeds)

    # Extract and create new rows for `rps_solar` and `rps_wind`
    additional_rows = []
    for carrier_col, carrier_name in [
        ("rps_solar", "solar"),
        ("rps_wind", "onwind, offwind, offwind_floating"),
    ]:
        if carrier_col in reeds.columns:
            temp = reeds[["region", "planning_horizon", carrier_col]].copy()
            temp = temp.rename(columns={carrier_col: "pct"})
            temp["carrier"] = carrier_name
            additional_rows.append(temp)

    # Combine original data with additional rows
    if additional_rows:
        additional_rows = pd.concat(additional_rows, ignore_index=True)
        reeds = pd.concat([reeds, additional_rows], ignore_index=True)

    # Ensure the final dataframe has consistent columns
    reeds = reeds[["region", "planning_horizon", "carrier", "pct"]]
    reeds = reeds[reeds["pct"] > 0.0]  # Remove any rows with zero or negative percentages

    return reeds


def load_reserves_data(n, config, snakemake):
    """
    Load and process reserves data including regional PRM requirements and ReEDS data.

    Parameters
    ----------
    n : pypsa.Network
        Network object containing bus information
    config : dict
        Configuration dictionary
    snakemake : object
        Snakemake object containing inputs

    Returns
    -------
    dict
        Dictionary containing reserves data
    """
    reserves_data = {}

    # Load user-defined PRM requirements
    reserves_data["regional_prm"] = pd.read_csv(
        config["electricity"]["SAFE_regional_reservemargins"],
        index_col=[0],
    )

    # Process ReEDS PRM data
    reeds_prm = pd.read_csv(snakemake.input.safer_reeds, index_col=[0])

    # Map NERC regions to ReEDS zones
    nerc_memberships = n.buses.groupby("nerc_reg")["reeds_zone"].apply(lambda x: ", ".join(x)).to_dict()

    reeds_prm["region"] = reeds_prm.index.map(nerc_memberships)
    reeds_prm = reeds_prm.dropna(subset="region")
    reeds_prm = reeds_prm.drop(
        columns=["none", "ramp2025_20by50", "ramp2025_25by50", "ramp2025_30by50"],
    )
    reeds_prm = reeds_prm.rename(columns={"static": "prm", "t": "planning_horizon"})

    reserves_data["reeds_prm"] = reeds_prm

    # Load operational reserve configuration if available
    if "operational_reserve" in config["electricity"]:
        reserves_data["operational_reserve"] = config["electricity"]["operational_reserve"]

    return reserves_data


def load_policy_data(network, config, snakemake):
    """
    Load all policy-related data from configuration files.

    Parameters
    ----------
    config : dict
        Configuration dictionary
    snakemake : object
        Snakemake object containing inputs and parameters

    Returns
    -------
    dict
        Dictionary containing all policy data
    """
    # Define carriers for RPS and CES
    rps_carriers = [
        "onwind",
        "offwind",
        "offwind_floating",
        "solar",
        "hydro",
        "geothermal",
        "biomass",
        "EGS",
    ]
    ces_carriers = [*rps_carriers, "nuclear", "SMR"]

    # Create partial function for loading CSV data
    load_policy_csv = partial(load_csv_data, config)

    # Load all policy data
    policy_data = {
        "technology_capacity_targets": load_policy_csv("technology_capacity_targets"),
        "portfolio_standards": load_policy_csv("portfolio_standards"),
        "regional_co2_limits": load_policy_csv("regional_Co2_limits", index_col=[0]),
        "rps_reeds": process_reeds_data(snakemake.input.rps_reeds, rps_carriers, "pct"),
        "ces_reeds": process_reeds_data(snakemake.input.ces_reeds, ces_carriers, "pct"),
    }

    # Load reserves data
    policy_data["reserves"] = load_reserves_data(network, config, snakemake)

    return policy_data


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_policy_data",
            # simpl="",
            clusters="100",
            interconnect="western",
            ll="v1.0",
            opts="500SEG",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)
    params = snakemake.params

    n = pypsa.Network(snakemake.input.network)

    # Load all policy data
    policy_data = load_policy_data(n, snakemake.config, snakemake)

    # Add policy data to network
    n.meta["policy_data"] = policy_data

    # Export network
    n.export_to_netcdf(snakemake.output.network)
