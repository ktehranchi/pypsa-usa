"""Generic module to add a new energy network

Creates new sector ontop of existing one. Note: Currently, can only be built ontop of electricity sector

Marginal costs are handeled as follows:
- Links are the VOM of just the technology
- Replacement generators contain time varrying fuel costs
"""

import pypsa
import pandas as pd
from typing import List, Union
import logging
logger = logging.getLogger(__name__)
from _helpers import configure_logging
from add_electricity import load_costs

def add_carrier(n: pypsa.Network, carrier: str, costs: pd.DataFrame = pd.DataFrame()):
    """Adds new carrier into the network with emission factor"""
    if not carrier in n.carriers.index:
        try:
            n.add("Carrier", carrier, co2_emissions=costs.at[carrier, "co2_emissions"])
        except KeyError:
            logger.debug(f"{carrier} does not have an emission factor ")
            n.add("Carrier", carrier, co2_emissions=0)

def add_buses(n: pypsa.Network, carrier: str, carrier_follow: str = "AC"):
    """Creates buses with new carriers at same geographic location
    
    carrier: str
        New carrier type 
    carrier_follow: str
        Carrier of buses to follow
    """
    buses = n.buses[n.buses.carrier == carrier_follow]
    if buses.empty:
        logger.debug(f"Can not create new buses for type {carrier}")
    else:
        n.madd(
            "Bus", 
            names=buses.index,
            suffix=f" {carrier}",
            x=buses.x,
            y=buses.y,
            interconnect=buses.interconnect,
            country=buses.country,
            zone_id=buses.zone_id,
            carrier=carrier
        )

def add_links(n: pypsa.Network, new_carrier: str, old_carriers: List[str], costs: pd.DataFrame = pd.DataFrame()):
    """Replace Generators with cross sector links
    
    n: pypsa.Network, 
    new_carrier: str,
        New carrier of the network (ie. Gas)
    old_carriers: List[str],
        Old carriers to replace (ie. ["CCGT", "OCGT"])
    costs: pd.DataFrame = pd.DataFrame(),
    """
    
    for old_carrier in old_carriers:
        plants = n.generators[n.generators.carrier==old_carrier]
        try:
            marginal_costs=costs.at[old_carrier, "VOM"]
        except KeyError:
            logger.warning(f"Can not locate marginal costs for {old_carrier}")
            marginal_costs=plants.marginal_cost
        n.madd(
            "Link",
            names=plants.index,
            bus0=plants.bus + f" {new_carrier}", # incorrect with anything except electrical network
            bus1=plants.bus,
            carrier=plants.carrier,
            p_nom_min=plants.p_nom_min,
            p_nom=plants.p_nom,
            p_nom_max=plants.p_nom_max,
            p_nom_extendable=plants.p_nom_extendable,
            ramp_limit_up=plants.ramp_limit_up,
            ramp_limit_down=plants.ramp_limit_down,
            efficiency=plants.efficiency,
            marginal_cost=marginal_costs,
            capital_cost=plants.capital_cost,
            lifetime=plants.lifetime,
        )
        n.mremove("Generator", plants.index)
        n.generators_t["marginal_cost"] = n.generators_t["marginal_cost"].drop(columns=[plants.index])

def add_generators(n: pypsa.Network, new_carrier: str, old_carriers: Union[str, List[str]], costs: pd.DataFrame = pd.DataFrame()):
    """Attaches generators to new buses"""
    for old_carrier in old_carriers:
        old_carrier_names = [old_carrier, f"{old_carrier} new"]
        plants = n.generators[n.generators.carrier.isin(old_carrier_names)].drop_duplicates(subset="bus")
        
        # check if plants have already been added by an old carrier 
        plants["new_name"] = plants.bus.map(lambda x: f"{x} {new_carrier}")
        existing_plants = n.generators.bus.drop_duplicates().to_list()
        plants = plants.loc[~plants.new_name.isin(existing_plants)]
        if plants.empty:
            continue
        
        try:
            marginal_costs=costs.at[old_carrier, "VOM"]
            fuel_costs=plants.marginal_cost - marginal_costs
        except KeyError:
            logger.warning(f"Can not locate marginal costs for {old_carrier}")
            marginal_costs=plants.marginal_cost
            fuel_costs=0
        
        # bus names are done differently cause its not working if you just do 
        # > bus = plants.bus + f" {new_carrier}"
        n.madd(
            "Generator", 
            names=plants.bus,
            suffix=f" {new_carrier}",
            bus=[f"{x} {new_carrier}" for x in plants.bus],
            carrier=new_carrier,
            p_nom_extendable=True,
            marginal_costs=fuel_costs
        )
        
        # copy over time dependent marginal costs 
        # this is not done in madd in case some plants do not have time-dependent costs
        name_mapper = {f"{bus} {old_carrier}":f"{bus} {new_carrier}" for bus in plants.bus}
        generators = n.generators_t.marginal_cost[[x for x in name_mapper]]
        if generators.empty:
            logger.info(f"No generators with time varrying marginal cost for {old_carrier}")
        else:
            generators = generators - marginal_costs # extract fuel costs 
            generators = generators.rename(columns=name_mapper)
            n.generators_t.marginal_cost = n.generators_t.marginal_cost.join(generators)

def remove_generators(n: pypsa.Network, old_carriers: Union[str, List[str]],):
    """Remove generators from old buses"""
    plants = n.generators[n.generators.carrier.isin(old_carriers)]
    n.mremove("Generator", plants.index)
    n.generators_t["marginal_cost"] = n.generators_t["marginal_cost"].drop(columns=[plants.index])

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_sectors",
            interconnect="western",
            # simpl="",
            clusters="30",
            ll="v1.25",
            opts="Co2L0.75",
            sectors="G",
        )
    configure_logging(snakemake)
    
    n = pypsa.Network(snakemake.input.network)

    params = snakemake.params
    
    Nyears = n.snapshot_weightings.objective.sum() / 8760.0
    costs = load_costs(
        snakemake.input.tech_costs,
        params.costs,
        params.electricity["max_hours"],
        Nyears,
    )

    sectors = snakemake.wildcards.sectors.split("-")
    
    if "G" in sectors:
        new_carrier = "gas"
        old_carriers = ["CCGT", "OCGT"]
        add_carrier(n, new_carrier, costs)
        add_buses(n, new_carrier)
        add_generators(n, new_carrier, old_carriers, costs)
        add_links(n, new_carrier, old_carriers, costs)
        remove_generators(n, old_carriers) # must come last 
        
    n.export_to_netcdf(snakemake.output.network)