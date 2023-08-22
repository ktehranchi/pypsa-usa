'''
Preprocesses Historical and Forecasted Load, Solar, and Wind Data

Written by Kamran Tehranchi, Stanford University.
'''
import pandas as pd, glob, os, logging, pypsa
from _helpers import progress_retrieve, configure_logging

def add_breakthrough_demand_from_file(n, fn_demand):

    """
    Zone power demand is disaggregated to buses proportional to Pd,
    where Pd is the real power demand (MW).
    """

    demand = pd.read_csv(fn_demand, index_col=0)
    demand.columns = demand.columns.astype(int)
    demand.index = n.snapshots

    intersection = set(demand.columns).intersection(n.buses.zone_id.unique())
    demand = demand[list(intersection)]

    demand_per_bus_pu = (n.buses.set_index("zone_id").Pd / n.buses.groupby("zone_id").sum().Pd)
    demand_per_bus = demand_per_bus_pu.multiply(demand)
    demand_per_bus.columns = n.buses.index

    n.madd( "Load", demand_per_bus.columns, bus=demand_per_bus.columns, p_set=demand_per_bus, carrier='AC')
    return n


def prepare_ads_load_data(ads_load_path, data_year):
    '''
    Modify ads load data to match the balancing authority names in the network. Need to test with all potential ADS years
    '''
    df_ads = pd.read_csv(ads_load_path)
    df_ads.columns = df_ads.columns.str.removeprefix('Load_')
    df_ads.columns = df_ads.columns.str.removesuffix('.dat')
    df_ads.columns = df_ads.columns.str.removesuffix(f'_{data_year}')
    df_ads.columns = df_ads.columns.str.removesuffix(f'_[18].dat: {data_year}')
    df_ads['CISO-PGAE'] = df_ads.pop('CIPV') + df_ads.pop('CIPB') + df_ads.pop('SPPC')#hotfix see github issue #15
    df_ads['BPAT'] = df_ads.pop('BPAT') + df_ads.pop('TPWR') + df_ads.pop('SCL')
    df_ads['IPCO'] = df_ads.pop('IPFE') + df_ads.pop('IPMV') + df_ads.pop('IPTV')
    df_ads['PACW'] = df_ads.pop('PAID') + df_ads.pop('PAUT') + df_ads.pop('PAWY')
    df_ads['Arizona'] = df_ads.pop('SRP') + df_ads.pop('AZPS') 
    df_ads.drop(columns=['Unnamed: 44', 'TH_Malin', 'TH_Mead', 'TH_PV'],inplace=True)
    ba_list_map = {'CISC': 'CISO-SCE', 'CISD': 'CISO-SDGE','VEA': 'CISO-VEA','WAUW': 'WAUW_SPP'}
    df_ads.rename(columns=ba_list_map,inplace=True)
    # df_ads['datetime'] = pd.Timestamp(f'{data_year}-01-01')+pd.to_timedelta(df_ads.index, unit='H')
    df_ads.set_index('Index',inplace=True)
    # if len(df_ads.index) > 8761: #remove leap year day
    #     df_ads= df_ads[~(df_ads.index.date == pd.to_datetime(f'{data_year}-04-29'))]

    # not_in_list = df_ads.loc[:,~df_ads.columns.isin(ba_list)]
    return df_ads

def add_ads_demand(n, demand):
    """
    Zone power demand is disaggregated to buses proportional to Pd,
    where Pd is the real power demand (MW).
    """
    demand.index = n.snapshots 
    # n.buses['ba_load_data'] = n.buses.balancing_area.replace({'CISO-PGAE': 'CISO', 'CISO-SCE': 'CISO', 'CISO-VEA': 'CISO', 'CISO-SDGE': 'CISO'})
    n.buses['ba_load_data'] = n.buses.balancing_area.replace({'': 'missing_ba'})

    intersection = set(demand.columns).intersection(n.buses.ba_load_data.unique())
    demand = demand[list(intersection)]

    demand_per_bus_pu = (n.buses.set_index("ba_load_data").Pd / n.buses.groupby("ba_load_data").sum().Pd)
    demand_per_bus = demand_per_bus_pu.multiply(demand)
    demand_per_bus.fillna(0,inplace=True)
    demand_per_bus.columns = n.buses.index

    n.madd( "Load", demand_per_bus.columns, bus=demand_per_bus.columns, p_set=demand_per_bus, carrier='AC') 
    return n


def add_eia_demand(n, demand):
    """
    Zone power demand is disaggregated to buses proportional to Pd,
    where Pd is the real power demand (MW).
    """
    demand.set_index('timestamp', inplace=True)
    demand.index = n.snapshots #maybe add check to make sure they match?

    demand['Arizona'] = demand.pop('SRP') + demand.pop('AZPS')
    n.buses['ba_load_data'] = n.buses.balancing_area.replace({'CISO-PGAE': 'CISO', 'CISO-SCE': 'CISO', 'CISO-VEA': 'CISO', 'CISO-SDGE': 'CISO'})
    n.buses['ba_load_data'] = n.buses.ba_load_data.replace({'': 'missing_ba'})

    intersection = set(demand.columns).intersection(n.buses.ba_load_data.unique())
    demand = demand[list(intersection)]

    demand_per_bus_pu = (n.buses.set_index("ba_load_data").Pd / n.buses.groupby("ba_load_data").sum().Pd)
    demand_per_bus = demand_per_bus_pu.multiply(demand)
    demand_per_bus.fillna(0,inplace=True)
    demand_per_bus.columns = n.buses.index

    n.madd( "Load", demand_per_bus.columns, bus=demand_per_bus.columns, p_set=demand_per_bus, carrier='AC') 

    return n
    

if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input['network'])

    if snakemake.config['network_configuration']  == "pypsa-usa":
        snapshot_config = snakemake.config['snapshots']
        load_year = pd.to_datetime(snapshot_config['start']).year
        logger.info(f'Building Load Data using EIA demand data year {load_year}')
        eia_demand = pd.read_csv(snakemake.input['eia'][load_year%2015])
        n.set_snapshots(pd.date_range(freq="h", 
                                      start=snapshot_config['start'],
                                      end=snapshot_config['end'],
                                      inclusive=snapshot_config['inclusive'])
                        )
        n = add_eia_demand(n, eia_demand)
    elif snakemake.config['network_configuration']  == "ads2032":
        load_year = 2032
        demand = prepare_ads_load_data(f'data/WECC_ADS/processed/load_{load_year}.csv',load_year)
        n.set_snapshots(pd.date_range(freq="h", 
                                      start=f"{load_year}-01-01",
                                      end=f"{load_year+1}-01-01",
                                      inclusive="left")
                        )
        n = add_ads_demand(n,demand)
        n.export_to_netcdf(snakemake.output.network)
    elif snakemake.config['network_configuration']  == "breakthrough":
        logger.info("Adding Breakthrough Energy Network Demand data from 2016")
        n.set_snapshots(
            pd.date_range(freq="h", start="2016-01-01", end="2017-01-01", inclusive="left")
        )
        n = add_breakthrough_demand_from_file(n, snakemake.input["demand_breakthrough_2016"])
    else:
        raise ValueError(f"Unknown network_configuration: {snakemake.config['network_configuration']}")

    n.export_to_netcdf(snakemake.output.network)