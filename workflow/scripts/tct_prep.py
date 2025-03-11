"""Prepared TCT for CH2 Project."""
import pypsa
from _helpers import configure_logging


def clean_statistics_csv(df, output_file):
    """Clean the statistics CSV file by renaming columns, dropping unnecessary columns, and renaming values in the columns."""
    df.columns = df.columns.get_level_values(0)
    df.columns = ["component", "region", "carrier", *list(df.columns[3:])]

    # Drop components with value 'Load' and 'Line'
    df = df[~df["component"].isin(["Load", "Line"])]

    # Assign the second column the name 'region' and take only the first 3 characters
    df["region"] = df["region"].str[:3]

    # Use a dictionary to rename the values in the third column and assign this new column the name 'carriers'
    rename_dict = {
        "Ccgt-95Ccs": "CCGT-95CCS",
        "Combined-Cycle Gas": "CCGT",
        "Hydrogen Combustion Turbine": "hydrogen_ct",
        "Nuclear": "nuclear",
        "Oil": "oil",
        "Onshore Wind": "onwind",
        "Open-Cycle Gas": "OCGT",
        "Solar": "solar",
        "Biomass": "biomass",
        "Coal": "coal",
        "Fixed Bottom Offshore Wind": "offwind",
        "Reservoir & Dam": "hydro",
        "Battery Storage": "battery_storage",
        "4Hr_Battery_Storage": "4hr_battery_storage",
    }
    df["carrier"] = df["carrier"].map(rename_dict).fillna(df["carrier"])

    # Keep only the 'carriers', 'region', and 'Optimal Capacity' columns
    df = df[["carrier", "region", "Optimal Capacity"]]

    # Create new columns 'name' and 'planning_horizon'
    df["name"] = "tamu"
    df["planning_horizon"] = snakemake.params.planning_horizons[0]

    df["min"] = ""
    df["max"] = ""
    df["equals"] = df["Optimal Capacity"]
    df = df.drop(columns=["Optimal Capacity"])
    df["equals"] = df["equals"].fillna(0)

    # Reorder columns
    df = df[["name", "region", "carrier", "planning_horizon", "min", "max", "equals"]]

    # Write the cleaned data to a new CSV file
    df.to_csv(output_file, index=False)


def generate_statistics(n):
    groupers = n.statistics.groupers
    return n.statistics(groupby=groupers.get_bus_and_carrier).round(3).reset_index()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_network",
            interconnect="western",
            simpl="70",
            clusters="4m",
            ll="v1.0",
            opts="1h-TCT",
            sector="E-G",
            planning_horizons="2030",
        )
    configure_logging(snakemake)
    network = pypsa.Network(snakemake.input.network)
    stats = generate_statistics(network)
    output_file = snakemake.output.mapping_csv
    clean_statistics_csv(stats, output_file)
