"""
Solves optimal operation and capacity for a network with the option to
iteratively optimize while updating line reactances.
"""

import logging

import numpy as np
import pandas as pd
import pypsa
import xarray as xr
import yaml
from _helpers import (
    configure_logging,
    update_config_from_wildcards,
    update_config_with_sector_opts,
)
from pypsa.descriptors import get_switchable_as_dense as get_as_dense

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


def get_region_buses(n, region_list):
    return n.buses[
        (
            n.buses.country.isin(region_list)
            | n.buses.reeds_zone.isin(region_list)
            | n.buses.reeds_state.isin(region_list)
            | n.buses.interconnect.str.lower().isin(region_list)
            | n.buses.nerc_reg.isin(region_list)
            | (1 if "all" in region_list else 0)
        )
    ]


def filter_components(
    n: pypsa.Network,
    component_type: str,
    planning_horizon: str | int,
    carrier_list: list[str],
    region_buses: pd.Index,
    extendable: bool,
):
    """
    Filter components based on common criteria.

    Parameters
    ----------
    - n: pypsa.Network
        The PyPSA network object.
    - component_type: str
        The type of component (e.g., "Generator", "StorageUnit").
    - planning_horizon: str or int
        The planning horizon to filter active assets.
    - carrier_list: list
        List of carriers to filter.
    - region_buses: pd.Index
        Index of region buses to filter.
    - extendable: bool, optional
        If specified, filters by extendable or non-extendable assets.

    Returns
    -------
    - pd.DataFrame
        Filtered assets.
    """
    component = n.df(component_type)
    if planning_horizon != "all":
        ph = int(planning_horizon)
        iv = n.investment_periods
        active_components = n.get_active_assets(component.index.name, iv[iv >= ph][0])
    else:
        active_components = component.index

    filtered = component.loc[
        active_components
        & component.carrier.isin(carrier_list)
        & component.bus.isin(region_buses)
        & (component.p_nom_extendable == extendable)
    ]

    return filtered


def add_land_use_constraints(n):
    """
    Adds constraint for land-use based on information from the generators
    table.

    Constraint is defined by land-use per carrier and land_region. The
    definition of land_region enables sub-bus level land-use
    constraints.
    """
    model = n.model
    generators = n.generators.query(
        "p_nom_extendable & land_region != '' ",
    ).rename_axis(index="Generator-ext")

    if generators.empty:
        return
    p_nom = n.model["Generator-p_nom"].loc[generators.index]

    grouper = pd.concat([generators.carrier, generators.land_region], axis=1)
    lhs = p_nom.groupby(grouper).sum()

    maximum = generators.groupby(["carrier", "land_region"])["p_nom_max"].max()
    maximum = maximum[np.isfinite(maximum)]

    rhs = xr.DataArray(maximum).rename(dim_0="group")
    index = rhs.indexes["group"].intersection(lhs.indexes["group"])

    if not index.empty:
        logger.info("Adding land-use constraints")
        model.add_constraints(
            lhs.sel(group=index) <= rhs.loc[index],
            name="land_use_constraint",
        )


def prepare_network(
    n,
    solve_opts=None,
):
    if "clip_p_max_pu" in solve_opts:
        for df in (
            n.generators_t.p_max_pu,
            n.generators_t.p_min_pu,
            n.storage_units_t.inflow,
        ):
            df = df.where(df > solve_opts["clip_p_max_pu"], other=0.0)

    load_shedding = solve_opts.get("load_shedding")
    load_shedding_gens = n.generators.query("carrier == 'load'")
    if load_shedding and load_shedding_gens.empty:
        # intersect between macroeconomic and surveybased willingness to pay
        # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
        # TODO: retrieve color and nice name from config
        n.add("Carrier", "load", color="#dd2e23", nice_name="Load shedding")
        buses_i = n.buses.query("carrier == 'AC'").index
        if not np.isscalar(load_shedding):
            # TODO: do not scale via sign attribute (use Eur/MWh instead of Eur/kWh)
            load_shedding = 1e2  # Eur/kWh

        n.madd(
            "Generator",
            buses_i,
            " load",
            bus=buses_i,
            carrier="load",
            sign=1e-3,  # Adjust sign to measure p and p_nom in kW instead of MW
            marginal_cost=load_shedding,  # Eur/kWh
            p_nom=1e9,  # kW
        )

    if solve_opts.get("noisy_costs"):
        for t in n.iterate_components():
            if "marginal_cost" in t.df:
                t.df["marginal_cost"] += 1e-2 + 2e-3 * (np.random.random(len(t.df)) - 0.5)

        for t in n.iterate_components(["Line", "Link"]):
            t.df["capital_cost"] += (1e-1 + 2e-2 * (np.random.random(len(t.df)) - 0.5)) * t.df["length"]

    if solve_opts.get("nhours"):
        nhours = solve_opts["nhours"]
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760.0 / nhours

    return n


def add_technology_capacity_target_constraints(n, config):
    """
    Add Technology Capacity Target (TCT) constraint to the network.

    Add minimum or maximum levels of generator nominal capacity per carrier for individual regions.
    Each constraint can be designated for a specified planning horizon in multi-period models.
    Opts and path for technology_capacity_targets.csv must be defined in config.yaml.
    Default file is available at config/policy_constraints/technology_capacity_targets.csv.

    Parameters
    ----------
    n : pypsa.Network
    config : dict

    Example
    -------
    scenario:
        opts: [Co2L-TCT-24H]
    electricity:
        technology_capacity_target: config/policy_constraints/technology_capacity_target.csv
    """
    p_nom = n.model["Generator-p_nom"]
    tct_data = pd.read_csv(config["electricity"]["technology_capacity_targets"])
    if tct_data.empty:
        return

    for _, target in tct_data.iterrows():
        planning_horizon = target.planning_horizon
        region_list = [region_.strip() for region_ in target.region.split(",")]
        carrier_list = [carrier_.strip() for carrier_ in target.carrier.split(",")]
        region_buses = get_region_buses(n, region_list)

        lhs_gens_ext = filter_components(
            n=n,
            component_type="Generator",
            planning_horizon=planning_horizon,
            carrier_list=carrier_list,
            region_buses=region_buses.index,
            extendable=True,
        )
        lhs_gens_existing = filter_components(
            n=n,
            component_type="Generator",
            planning_horizon=planning_horizon,
            carrier_list=carrier_list,
            region_buses=region_buses.index,
            extendable=False,
        )

        lhs_storage_ext = filter_components(
            n=n,
            component_type="StorageUnit",
            planning_horizon=planning_horizon,
            carrier_list=carrier_list,
            region_buses=region_buses.index,
            extendable=True,
        )
        lhs_storage_existing = filter_components(
            n=n,
            component_type="StorageUnit",
            planning_horizon=planning_horizon,
            carrier_list=carrier_list,
            region_buses=region_buses.index,
            extendable=False,
        )

        if region_buses.empty or (lhs_gens_ext.empty and lhs_storage_ext.empty):
            continue

        if not lhs_gens_ext.empty:
            grouper_g = pd.concat(
                [lhs_gens_ext.bus.map(n.buses.country), lhs_gens_ext.carrier],
                axis=1,
            ).rename_axis(
                "Generator-ext",
            )
            lhs_g = p_nom.loc[lhs_gens_ext.index].groupby(grouper_g).sum().rename(bus="country")
        else:
            lhs_g = None

        if not lhs_storage_ext.empty:
            grouper_s = pd.concat(
                [lhs_storage_ext.bus.map(n.buses.country), lhs_storage_ext.carrier],
                axis=1,
            ).rename_axis(
                "StorageUnit-ext",
            )
            lhs_s = n.model["StorageUnit-p_nom"].loc[lhs_storage_ext.index].groupby(grouper_s).sum()
        else:
            lhs_s = None

        if lhs_g is None and lhs_s is None:
            continue
        elif lhs_g is None:
            lhs = lhs_s.sum()
        elif lhs_s is None:
            lhs = lhs_g.sum()
        else:
            lhs = (lhs_g + lhs_s).sum()

        lhs_existing = lhs_gens_existing.p_nom.sum() + lhs_storage_existing.p_nom.sum()

        if target["max"] == "existing":
            target["max"] = round(lhs_existing, 2) + 0.01
        else:
            target["max"] = float(target["max"])

        if target["min"] == "existing":
            target["min"] = round(lhs_existing, 2) - 0.01
        else:
            target["min"] = float(target["min"])

        if not np.isnan(target["min"]):
            rhs = target["min"] - round(lhs_existing, 2)

            n.model.add_constraints(
                lhs >= rhs,
                name=f"GlobalConstraint-{target.name}_{target.planning_horizon}_min",
            )

            logger.info(
                "Adding TCT Constraint:\n"
                f"Name: {target.name}\n"
                f"Planning Horizon: {target.planning_horizon}\n"
                f"Region: {target.region}\n"
                f"Carrier: {target.carrier}\n"
                f"Min Value: {target['min']}\n"
                f"Min Value Adj: {rhs}",
            )

        if not np.isnan(target["max"]):
            assert (
                target["max"] >= lhs_existing
            ), f"TCT constraint of {target['max']} MW for {target['carrier']} must be at least {lhs_existing}"

            rhs = target["max"] - round(lhs_existing, 2)

            n.model.add_constraints(
                lhs <= rhs,
                name=f"GlobalConstraint-{target.name}_{target.planning_horizon}_max",
            )

            logger.info(
                "Adding TCT Constraint:\n"
                f"Name: {target.name}\n"
                f"Planning Horizon: {target.planning_horizon}\n"
                f"Region: {target.region}\n"
                f"Carrier: {target.carrier}\n"
                f"Max Value: {target['max']}\n"
                f"Max Value Adj: {rhs}",
            )

        if not np.isnan(target["equals"]):
            rhs = target["equals"] - round(lhs_existing, 2)

            n.model.add_constraints(
                lhs <= rhs,
                name=f"GlobalConstraint-{target.name}_{target.planning_horizon}_equals",
            )

            logger.info(
                "Adding TCT Constraint:\n"
                f"Name: {target.name}\n"
                f"Planning Horizon: {target.planning_horizon}\n"
                f"Region: {target.region}\n"
                f"Carrier: {target.carrier}\n"
                f"Max Value: {target['max']}\n"
                f"Max Value Adj: {rhs}",
            )


def add_RPS_constraints(n, config):
    """
    Add Renewable Portfolio Standards (RPS) constraints to the network.

    This function enforces constraints on the percentage of electricity generation
    from renewable energy sources for specific regions and planning horizons.
    It reads the necessary data from configuration files and the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object.
    config : dict
        A dictionary containing configuration settings and file paths.

    Returns
    -------
    None
    """

    def process_reeds_data(filepath, carriers, value_col):
        """Helper function to process RPS or CES REEDS data."""
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

    # Read portfolio standards data
    portfolio_standards = pd.read_csv(config["electricity"]["portfolio_standards"])

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

    # Process RPS and CES REEDS data
    rps_reeds = process_reeds_data(
        snakemake.input.rps_reeds,
        rps_carriers,
        value_col="pct",
    )
    ces_reeds = process_reeds_data(
        snakemake.input.ces_reeds,
        ces_carriers,
        value_col="pct",
    )

    # Concatenate all portfolio standards
    portfolio_standards = pd.concat([portfolio_standards, rps_reeds, ces_reeds])
    portfolio_standards = portfolio_standards[
        (portfolio_standards.pct > 0.0)
        & (
            portfolio_standards.planning_horizon.isin(
                snakemake.params.planning_horizons,
            )
        )
        & (portfolio_standards.region.isin(n.buses.reeds_state.unique()))
    ]

    # Iterate through constraints and add RPS constraints to the model
    for _, constraint_row in portfolio_standards.iterrows():
        region_list = [region.strip() for region in constraint_row.region.split(",")]
        region_buses = get_region_buses(n, region_list)

        if region_buses.empty:
            continue

        carriers = [carrier.strip() for carrier in constraint_row.carrier.split(",")]

        # Filter region generators
        region_gens = n.generators[n.generators.bus.isin(region_buses.index)]
        region_gens_eligible = region_gens[region_gens.carrier.isin(carriers)]

        if not region_gens_eligible.empty:
            # Eligible generation
            p_eligible = n.model["Generator-p"].sel(
                period=constraint_row.planning_horizon,
                Generator=region_gens_eligible.index,
            )
            lhs = p_eligible.sum()

            # Region demand
            region_demand = (
                n.loads_t.p_set.loc[
                    constraint_row.planning_horizon,
                    n.loads.bus.isin(region_buses.index),
                ]
                .sum()
                .sum()
            )

            rhs = constraint_row.pct * region_demand

            # Add constraint
            n.model.add_constraints(
                lhs >= rhs,
                name=f"GlobalConstraint-{constraint_row.name}_{constraint_row.planning_horizon}_rps_limit",
            )
            logger.info(
                f"Added RPS {constraint_row.name} for {constraint_row.planning_horizon}.",
            )


def add_regional_co2limit(n, sns, config):
    """Adding regional regional CO2 Limits Specified in the config.yaml."""
    regional_co2_lims = pd.read_csv(
        config["electricity"]["regional_Co2_limits"],
        index_col=[0],
    )
    logger.info("Adding regional Co2 Limits.")
    regional_co2_lims = regional_co2_lims[regional_co2_lims.planning_horizon.isin(snakemake.params.planning_horizons)]
    weightings = n.snapshot_weightings.loc[n.snapshots]

    for idx, emmission_lim in regional_co2_lims.iterrows():
        region_list = [region.strip() for region in emmission_lim.regions.split(",")]
        region_buses = get_region_buses(n, region_list)

        emissions = n.carriers.co2_emissions.fillna(0)[lambda ds: ds != 0]
        region_gens = n.generators[n.generators.bus.isin(region_buses.index)]
        region_gens_em = region_gens.query("carrier in @emissions.index")

        if region_buses.empty or region_gens_em.empty:
            continue

        region_co2lim = emmission_lim.limit
        planning_horizon = emmission_lim.planning_horizon

        efficiency = get_as_dense(
            n,
            "Generator",
            "efficiency",
            inds=region_gens_em.index,
        )  # mw_elect/mw_th
        em_pu = region_gens_em.carrier.map(emissions) / efficiency  # tonnes_co2/mw_electrical
        em_pu = em_pu.multiply(weightings.generators, axis=0).loc[planning_horizon].fillna(0)

        # Emitting Gens
        p_em = n.model["Generator-p"].loc[:, region_gens_em.index].sel(period=planning_horizon)
        lhs = (p_em * em_pu).sum()
        rhs = region_co2lim

        n.model.add_constraints(
            lhs <= rhs,
            name=f"GlobalConstraint-{emmission_lim.name}_{planning_horizon}co2_limit",
        )

        logger.info(
            f"Adding regional Co2 Limit for {emmission_lim.name} in {planning_horizon}",
        )


def add_PRM_constraints(n, config):
    """
    Add Planning Reserve Margin (PRM) constraints for regional capacity adequacy.

    This function enforces that each region has sufficient firm capacity to meet
    peak demand plus a reserve margin. Only firm resources (not variable renewables
    or storage) contribute to meeting this requirement.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object
    config : dict
        Configuration dictionary containing PRM parameters
    """
    # Load regional PRM requirements
    regional_prm = _get_combined_prm_requirements(n, config)

    # Apply constraints for each region and planning horizon
    for _, prm in regional_prm.iterrows():
        # Skip if no valid planning horizon or region
        if prm.planning_horizon not in n.investment_periods:
            continue

        region_list = [region_.strip() for region_ in prm.region.split(",")]
        region_buses = get_region_buses(n, region_list)

        if region_buses.empty:
            continue

        # Calculate peak demand and required reserve margin
        regional_demand = _get_regional_demand(n, prm.planning_horizon, region_buses)
        peak_demand = regional_demand.max()
        planning_reserve = peak_demand * (1.0 + prm.prm)

        # Get capacity contribution from resources
        lhs_capacity = _calculate_capacity_accredidation(
            n,
            prm.planning_horizon,
            region_buses,
            peak_demand_hour=regional_demand.idxmax(),
        )

        # Add the constraint to the model
        n.model.add_constraints(
            lhs_capacity >= planning_reserve,
            name=f"GlobalConstraint-{prm.name}_{prm.planning_horizon}_PRM",
        )

        logger.info(
            f"Added PRM constraint for {prm.name} in {prm.planning_horizon}: "
            f"Peak demand: {peak_demand:.2f} MW, "
            f"Required capacity: {planning_reserve:.2f} MW",
        )


def _get_combined_prm_requirements(n, config):
    """
    Combine PRM requirements from different sources into a single dataframe.

    Parameters
    ----------
    n : pypsa.Network
    config : dict

    Returns
    -------
    pd.DataFrame
        Combined PRM requirements with columns: name, region, prm, planning_horizon
    """
    # Load user-defined PRM requirements
    regional_prm = pd.read_csv(
        config["electricity"]["SAFE_regional_reservemargins"],
        index_col=[0],
    )

    # Process ReEDS PRM data if available
    try:
        reeds_prm = pd.read_csv(snakemake.input.safer_reeds, index_col=[0])

        # Map NERC regions to ReEDS zones
        nerc_memberships = (
            n.buses.groupby("nerc_reg")["reeds_zone"]
            .apply(
                lambda x: ", ".join(x),
            )
            .to_dict()
        )

        reeds_prm["region"] = reeds_prm.index.map(nerc_memberships)
        reeds_prm = reeds_prm.dropna(subset="region")
        reeds_prm = reeds_prm.drop(
            columns=["none", "ramp2025_20by50", "ramp2025_25by50", "ramp2025_30by50"],
        )
        reeds_prm = reeds_prm.rename(columns={"static": "prm", "t": "planning_horizon"})

        # Combine both data sources
        regional_prm = pd.concat([regional_prm, reeds_prm])
    except (FileNotFoundError, AttributeError):
        logger.info("ReEDS PRM data not available, using only user-defined PRM values")

    # Filter for relevant planning horizons
    return regional_prm[regional_prm.planning_horizon.isin(n.investment_periods)]


def _get_regional_demand(n, planning_horizon, region_buses):
    """
    Calculate hourly demand for a specific region and planning horizon.

    Parameters
    ----------
    n : pypsa.Network
    planning_horizon : int or str
        Planning horizon year
    region_buses : pd.DataFrame
        DataFrame containing buses in the region

    Returns
    -------
    pd.Series
        Hourly demand series for the region
    """
    return n.loads_t.p_set.loc[
        planning_horizon,
        n.loads.bus.isin(region_buses.index),
    ].sum(axis=1)


#  n.loads_t.p_set.loc[planning_horizon, n.loads.bus.isin(region_buses.index)].sum(axis=1)
def _calculate_capacity_accredidation(n, planning_horizon, region_buses, peak_demand_hour):
    """
    Calculate firm capacity contribution from all resources in a region.

    This function accounts for:
    1. Extendable resources with appropriate capacity credit
    2. Non-extendable existing resources

    Parameters
    ----------
    n : pypsa.Network
    planning_horizon : int or str
    region_buses : pd.DataFrame
    peak_demand_hour : pd.Timestamp
        Hour of peak demand used for calculating capacity credits

    Returns
    -------
    float or xarray.DataArray
        Total firm capacity contribution
    """
    # Get active generators during this planning period
    active_gens = n.get_active_assets("Generator", planning_horizon)
    extendable_gens = n.generators.p_nom_extendable
    region_gens = n.generators.bus.isin(region_buses.index)

    # Extendable capacity with capacity credit
    region_active_ext_gens = region_gens & active_gens & extendable_gens
    region_active_ext_gens = n.generators[region_active_ext_gens]

    if not region_active_ext_gens.empty:
        ext_p_nom = n.model["Generator-p_nom"].loc[region_active_ext_gens.index]
        ext_p_max_pu = get_as_dense(
            n,
            "Generator",
            "p_max_pu",
            inds=region_active_ext_gens.index,
        ).loc[
            planning_horizon,
            peak_demand_hour,
        ]
        ext_contribution = ext_p_nom * ext_p_max_pu
    else:
        ext_contribution = 0

    # Non-extendable existing capacity
    region_active_nonext_gens = region_gens & active_gens & ~extendable_gens
    region_active_nonext_gens = n.generators[region_active_nonext_gens]

    if not region_active_nonext_gens.empty:
        non_ext_p_max_pu = get_as_dense(
            n,
            "Generator",
            "p_max_pu",
            inds=region_active_nonext_gens.index,
        ).loc[
            planning_horizon,
            peak_demand_hour,
        ]
        non_ext_p_nom = region_active_nonext_gens.p_nom
        non_ext_contribution = (non_ext_p_nom * non_ext_p_max_pu).sum()
    else:
        non_ext_contribution = 0

    return ext_contribution.sum() + non_ext_contribution


def add_operational_reserve_margin(n, sns, config):
    """
    Build reserve margin constraints based on the formulation given in
    https://genxproject.github.io/GenX/dev/core/#Reserves.

    Parameters
    ----------
        n : pypsa.Network
        sns: pd.DatetimeIndex
        config : dict

    Example:
    --------
    config.yaml requires to specify operational_reserve:
    operational_reserve: # like https://genxproject.github.io/GenX/dev/core/#Reserves
        activate: true
        epsilon_load: 0.02 # percentage of load at each snapshot
        epsilon_vres: 0.02 # percentage of VRES at each snapshot
        contingency: 400000 # MW
    """
    reserve_config = config["electricity"]["operational_reserve"]
    eps_load = reserve_config["epsilon_load"]
    eps_vres = reserve_config["epsilon_vres"]
    contingency = reserve_config["contingency"]

    # Reserve Variables
    n.model.add_variables(
        0,
        np.inf,
        coords=[sns, n.generators.index],
        name="Generator-r",
    )
    reserve = n.model["Generator-r"]
    summed_reserve = reserve.sum("Generator")

    # Share of extendable renewable capacities
    ext_i = n.generators.query("p_nom_extendable").index
    vres_i = n.generators_t.p_max_pu.columns
    if not ext_i.empty and not vres_i.empty:
        capacity_factor = n.generators_t.p_max_pu[vres_i.intersection(ext_i)]
        p_nom_vres = n.model["Generator-p_nom"].loc[vres_i.intersection(ext_i)].rename({"Generator-ext": "Generator"})
        lhs = summed_reserve + (p_nom_vres * (-eps_vres * capacity_factor)).sum(
            "Generator",
        )
    else:  # if no extendable VRES
        lhs = summed_reserve

    # Total demand per t
    demand = get_as_dense(n, "Load", "p_set").sum(axis=1)

    # VRES potential of non extendable generators
    capacity_factor = n.generators_t.p_max_pu[vres_i.difference(ext_i)]
    renewable_capacity = n.generators.p_nom[vres_i.difference(ext_i)]
    potential = (capacity_factor * renewable_capacity).sum(axis=1)

    # Right-hand-side
    rhs = eps_load * demand + eps_vres * potential + contingency

    n.model.add_constraints(lhs >= rhs, name="reserve_margin")

    # additional constraint that capacity is not exceeded
    gen_i = n.generators.index
    ext_i = n.generators.query("p_nom_extendable").index
    fix_i = n.generators.query("not p_nom_extendable").index

    dispatch = n.model["Generator-p"]
    reserve = n.model["Generator-r"]

    capacity_fixed = n.generators.p_nom[fix_i]

    p_max_pu = get_as_dense(n, "Generator", "p_max_pu")

    if not ext_i.empty:
        capacity_variable = n.model["Generator-p_nom"].rename(
            {"Generator-ext": "Generator"},
        )
        lhs = dispatch + reserve - capacity_variable * p_max_pu[ext_i]
    else:
        lhs = dispatch + reserve

    rhs = (p_max_pu[fix_i] * capacity_fixed).reindex(columns=gen_i, fill_value=0)

    n.model.add_constraints(lhs <= rhs, name="Generator-p-reserve-upper")


# def remove_kvl(n):
#     """
#     Removes Kirchhoff's voltage law (KVL) constraints.

#     Function implemented for Kamran's research, and not added to default configs.
#     """
#     n.model.constraints.remove("Kirchhoff-Voltage-Law")


def extra_functionality(n, snapshots):
    """
    Collects supplementary constraints which will be passed to
    ``pypsa.optimization.optimize``.

    If you want to enforce additional custom constraints, this is a good
    location to add them. The arguments ``opts`` and
    ``snakemake.config`` are expected to be attached to the network.
    """
    opts = n.opts
    config = n.config
    if "RPS" in opts and n.generators.p_nom_extendable.any():
        add_RPS_constraints(n, config)
    if "REM" in opts and n.generators.p_nom_extendable.any():
        add_regional_co2limit(n, snapshots, config)
    if "PRM" in opts and n.generators.p_nom_extendable.any():
        add_PRM_constraints(n, config)
    if "TCT" in opts and n.generators.p_nom_extendable.any():
        add_technology_capacity_target_constraints(n, config)
    reserve = config["electricity"].get("operational_reserve", {})
    if reserve.get("activate"):
        add_operational_reserve_margin(n, snapshots, config)
    # if config.get("solving", {}).get("options", {}).get("remove_kvl", False):
    #     remove_kvl(n)
    add_land_use_constraints(n)


def solve_network(n, config, solving, opts="", **kwargs):
    set_of_options = solving["solver"]["options"]
    cf_solving = solving["options"]

    if len(n.investment_periods) > 1:
        kwargs["multi_investment_periods"] = config["foresight"] == "perfect"

    kwargs["solver_options"] = solving["solver_options"][set_of_options] if set_of_options else {}
    kwargs["solver_name"] = solving["solver"]["name"]
    kwargs["extra_functionality"] = extra_functionality
    kwargs["transmission_losses"] = cf_solving.get("transmission_losses", False)
    kwargs["linearized_unit_commitment"] = cf_solving.get(
        "linearized_unit_commitment",
        False,
    )
    kwargs["assign_all_duals"] = cf_solving.get("assign_all_duals", False)

    rolling_horizon = cf_solving.pop("rolling_horizon", False)
    skip_iterations = cf_solving.get("skip_iterations", False)
    if not n.lines.s_nom_extendable.any():
        skip_iterations = True
        logger.info("No expandable lines found. Skipping iterative solving.")

    # add to network for extra_functionality
    n.config = config
    n.opts = opts

    if rolling_horizon:
        kwargs["horizon"] = cf_solving.get("horizon", 365)
        kwargs["overlap"] = cf_solving.get("overlap", 0)
        n.optimize.optimize_with_rolling_horizon(**kwargs)
        status, condition = "", ""
    elif skip_iterations:
        status, condition = n.optimize(**kwargs)
    else:
        kwargs["track_iterations"] = cf_solving.get("track_iterations", False)
        kwargs["min_iterations"] = cf_solving.get("min_iterations", 4)
        kwargs["max_iterations"] = cf_solving.get("max_iterations", 6)
        status, condition = n.optimize.optimize_transmission_expansion_iteratively(
            **kwargs,
        )

    if status != "ok" and not rolling_horizon:
        logger.warning(
            f"Solving status '{status}' with termination condition '{condition}'",
        )
    if "infeasible" in condition:
        # n.model.print_infeasibilities()
        raise RuntimeError("Solving status 'infeasible'")

    return n


def solve_network_iterative(n, config, solving, opts="", **kwargs):
    """
    Iteratively solves the network alternating between Generation Expansion and
    Transmission Expansion.
    """
    logger.info("Solving network iteratively")
    metrics = []

    # Remove load shedding and global constraints
    load_shedding_gens = n.generators.query("carrier == 'load'")
    if not load_shedding_gens.empty and not config["solving"]["options"]["load_shedding"]:
        n.mremove("Generator", load_shedding_gens.index)
    n.global_constraints.drop(n.global_constraints.index, inplace=True)

    gens_extensible_mask = n.generators.p_nom_extendable
    storage_extensible_mask = n.storage_units.p_nom_extendable
    links_extensible_mask = n.links.p_nom_extendable

    # Set minimum capacities
    n.lines.s_nom_min = n.lines.s_nom
    n.links.p_nom_min = n.links.p_nom

    def track_metrics(iter_num, expansion_type):
        """Helper function to collect metrics at each iteration."""
        # Get optimal capacity statistics
        capacity_stats = n.statistics.optimal_capacity()
        lv0_stats = capacity_stats.groupby(level=0).sum()
        capacity_stats = capacity_stats.droplevel(0)
        capacity_stats.fillna(0, inplace=True)

        opex = n.statistics.opex().groupby(level=0).sum() / 1e6
        capex = n.statistics.capex().groupby(level=0).sum() / 1e6
        # append 'capex' and 'opex' to the metrics index names strings
        capex.index = capex.index.map(lambda x: x + "_capex")
        opex.index = opex.index.map(lambda x: x + "_opex")
        lv0_stats.index = lv0_stats.index.map(lambda x: x + " _summary_caps")

        # Create a dictionary for the metrics
        return {
            "iteration": iter_num,
            "expansion_type": expansion_type,
            "objective": n.objective + n.objective_constant,
            # unpack the stats to dicts
            **capex[n.investment_periods[0]].to_dict(),
            **opex[n.investment_periods[0]].to_dict(),
            **lv0_stats[n.investment_periods[0]].to_dict(),
            **capacity_stats[n.investment_periods[0]].to_dict(),
        }

    new_metrics = track_metrics(0, "Base")
    metrics.append(new_metrics)
    logger.info(f"New Metrics: {new_metrics}")

    for iter_ in range(1, config["iterative_solving"]["max_iter"] + 1):
        logger.info(f"SOLUTION ITERATION {iter_}")
        logger.info("Transmission Expansion")
        # Transmission Expansion step
        n.generators.p_nom = n.generators.p_nom_opt
        n.storage_units.p_nom = n.storage_units.p_nom_opt

        # Fix generation, allow transmission expansion
        n.generators.p_nom_extendable = False
        n.storage_units.p_nom_extendable = False

        n.lines.s_nom_extendable = True
        n.links.p_nom_extendable = links_extensible_mask

        n.lines.s_nom_max = 1e8
        n.links.p_nom_max = 1e8

        n = solve_network(n, config, solving, opts, **kwargs)
        new_metrics = track_metrics(iter_, "Transmission")
        metrics.append(new_metrics)
        logger.info(f"New Metrics: {new_metrics}")

        # Generation Expansion step
        logger.info("Generation Expansion")
        # Fix transmission, allow generation expansion
        n.lines.s_nom = n.lines.s_nom_opt
        n.links.p_nom = n.links.p_nom_opt
        n.lines.s_nom_extendable = False
        n.links.p_nom_extendable = False

        n.generators.p_nom_extendable = gens_extensible_mask
        n.storage_units.p_nom_extendable = storage_extensible_mask

        n = solve_network(n, config, solving, opts, **kwargs)
        new_metrics = track_metrics(iter_, "Generation")
        metrics.append(new_metrics)
        logger.info(f"New Metrics: {new_metrics}")

        if config["iterative_solving"]["TEP_only"]:
            break

        # Define the stopping criteria
        if iter_ > 1:
            # Compare the last two generation expansions
            if np.allclose(
                metrics[-1]["total_cost"],
                metrics[-3]["total_cost"],
                rtol=config["iterative_solving"]["rtol"],
            ):
                logger.info(f"Convergence reached after {iter_} iterations")
                break

    return n, pd.DataFrame(metrics)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_network_iterative",
            interconnect="western",
            simpl="70",
            clusters="4m",
            ll="v1.0",
            opts="1h-TCT",
            sector="E-G",
            planning_horizons="2030",
        )
    configure_logging(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)
    if "sector_opts" in snakemake.wildcards.keys():
        update_config_with_sector_opts(
            snakemake.config,
            snakemake.wildcards.sector_opts,
        )

    opts = snakemake.wildcards.opts

    if "iterative" in snakemake.params.keys():
        snakemake.config["solving"]["options"]["load_shedding"] = True
        snakemake.params["solving"]["options"]["load_shedding"] = True

    if "sector_opts" in snakemake.wildcards.keys():
        opts += "-" + snakemake.wildcards.sector_opts
    opts = [o for o in opts.split("-") if o != ""]
    solve_opts = snakemake.params.solving["options"]

    # sector specific co2 options
    if snakemake.wildcards.sector != "E":
        # sector co2 limits applied via config file, not through Co2L
        opts = [x for x in opts if not x.startswith("Co2L")]
        opts.append("sector")

    np.random.seed(solve_opts.get("seed", 123))

    n = pypsa.Network(snakemake.input.network)

    n = prepare_network(
        n,
        solve_opts,
    )

    n, metrics = solve_network_iterative(
        n,
        config=snakemake.config,
        solving=snakemake.params.solving,
        opts=opts,
        log_fn=snakemake.log.solver,
    )

    metrics.to_csv(snakemake.output.iterative_metrics)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])

    with open(snakemake.output.config, "w") as file:
        yaml.dump(
            n.meta,
            file,
            default_flow_style=False,
            allow_unicode=True,
            sort_keys=False,
        )
