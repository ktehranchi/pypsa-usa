# Rules to Optimize/Solve Network


def pop_layout_input(wildcards):
    if wildcards["sector"] != "E":
        return RESOURCES + "{interconnect}/pop_layout_elec_s{simpl}_c{clusters}.csv"
    else:
        return []


rule generate_tct_requirements:
    params:
        planning_horizons=config["scenario"]["planning_horizons"],
    input:
        network=RESULTS
        + "{interconnect}/networks/elec_s{simpl}_c{clusters}_ec_l{ll}_{opts}_{sector}.nc",
    output:
        mapping_csv=RESULTS
        + "{interconnect}/figures/s{simpl}_cluster_{clusters}/l{ll}_{opts}_{sector}/statistics/tct_inputs.csv",
    log:
        "logs/generate_tct_requirements/{interconnect}/figures/s{simpl}_cluster_{clusters}/l{ll}_{opts}_{sector}.log",
    script:
        "../scripts/tct_prep.py"


rule solve_network_mapping:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        planning_horizons=config["scenario"]["planning_horizons"],
        co2_sequestration_potential=config["sector"].get(
            "co2_sequestration_potential", 200
        ),
        transmission_network=config_provider("model_topology", "transmission_network"),
        mapping=config_provider("mapping"),
    input:
        network=RESOURCES
        + "{interconnect}/elec_s{simpl}_ch{clusters_hires}_ec_l{ll}_{opts}.nc",
        flowgates="repo_data/ReEDS_Constraints/transmission/transmission_capacity_init_AC_ba_NARIS2024.csv",
        safer_reeds="config/policy_constraints/reeds/prm_annual.csv",
        rps_reeds="config/policy_constraints/reeds/rps_fraction.csv",
        ces_reeds="config/policy_constraints/reeds/ces_fraction.csv",
        pop_layout=pop_layout_input,
        mapping_csv=RESULTS
        + "{interconnect}/figures/s{simpl}_cluster_{clusters}/l{ll}_{opts}_{sector}/statistics/tct_inputs.csv",
    output:
        network=RESULTS
        + "{interconnect}/networks/elec_s{simpl}_cl{clusters}_ch{clusters_hires}_ec_l{ll}_{opts}_{sector}_mapped.nc",
        config=RESULTS
        + "{interconnect}/configs/config.elec_s{simpl}_cl{clusters}_ch{clusters_hires}_l{ll}_{opts}_{sector}.yaml",
    log:
        solver=normpath(
            LOGS
            + "solve_network/{interconnect}/elec_s{simpl}_cl_{clusters}_ch{clusters_hires}_ec_l{ll}_{opts}_{sector}_solver.log"
        ),
        python=LOGS
        + "solve_network/{interconnect}/elec_s{simpl}_cl_{clusters}_ch{clusters_hires}_ec_l{ll}_{opts}_{sector}_python.log",
    threads: solver_threads
    resources:
        mem_mb=lambda wildcards, input, attempt: (input.size // 100000) * attempt * 80,
        walltime=config["solving"].get("walltime", "12:00:00"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"


rule solve_network_iterative:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        planning_horizons=config["scenario"]["planning_horizons"],
        co2_sequestration_potential=config["sector"].get(
            "co2_sequestration_potential", 200
        ),
        transmission_network=config_provider("model_topology", "transmission_network"),
        iterative=True,
    input:
        network=(
            RESULTS
            + "{interconnect}/networks/elec_s{simpl}_c{clusters}_ec_l{ll}_{opts}_{sector}.nc"
            if config_provider("mapping", "skip_mapping")
            else RESULTS
            + "{interconnect}/networks/elec_s{simpl}_cl{clusters}_ch{clusters_hires}_ec_l{ll}_{opts}_{sector}_mapped.nc"
        ),
        flowgates="repo_data/ReEDS_Constraints/transmission/transmission_capacity_init_AC_ba_NARIS2024.csv",
        safer_reeds="config/policy_constraints/reeds/prm_annual.csv",
        rps_reeds="config/policy_constraints/reeds/rps_fraction.csv",
        ces_reeds="config/policy_constraints/reeds/ces_fraction.csv",
        pop_layout=pop_layout_input,
    output:
        network=RESULTS
        + "{interconnect}/networks/elec_s{simpl}_cl{clusters}_ch{clusters_hires}_ec_l{ll}_{opts}_{sector}_mapped_TEP.nc",
        config=RESULTS
        + "{interconnect}/configs/config.elec_s{simpl}_cl{clusters}_ch{clusters_hires}_l{ll}_{opts}_{sector}_mapped_TEP.yaml",
        iterative_metrics=RESULTS
        + "{interconnect}/elec_s{simpl}_cl{clusters}_ch{clusters_hires}_ec_l{ll}_{opts}_{sector}iterative_metrics.csv",
    log:
        solver=normpath(
            LOGS
            + "solve_network/{interconnect}/elec_s{simpl}_cl{clusters}_ch{clusters_hires}_ec_l{ll}_{opts}_{sector}_solver_iter.log"
        ),
        python=LOGS
        + "solve_network/{interconnect}/elec_s{simpl}_cl{clusters}_ch{clusters_hires}_ec_l{ll}_{opts}_{sector}_python_iter.log",
    benchmark:
        (
            BENCHMARKS
            + "solve_network/{interconnect}/elec_s{simpl}_cl{clusters}_ch{clusters_hires}_ec_l{ll}_{opts}_{sector}"
        )
    threads: solver_threads
    resources:
        mem_mb=lambda wildcards, input, attempt: (input.size // 100000) * attempt * 80,
        walltime=config["solving"].get("walltime", "12:00:00"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network_iterative.py"


##### inputs for hires file #####55


rule cluster_network_hires:
    params:
        cluster_network=config_provider("clustering", "cluster_network"),
        conventional_carriers=config_provider("electricity", "conventional_carriers"),
        renewable_carriers=config_provider("electricity", "renewable_carriers"),
        aggregation_strategies=config_provider("clustering", "aggregation_strategies"),
        custom_busmap=config_provider("enable", "custom_busmap", default=False),
        focus_weights=config_provider("focus_weights", default=False),
        length_factor=config_provider("lines", "length_factor"),
        costs=config_provider("costs"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        transmission_network=config_provider("model_topology", "transmission_network"),
        topological_boundaries=config_provider(
            "model_topology", "topological_boundaries"
        ),
        topology_aggregation=config_provider("model_topology", "aggregate"),
    input:
        network=RESOURCES + "{interconnect}/elec_s{simpl}.nc",
        regions_onshore=RESOURCES
        + "{interconnect}/Geospatial/regions_onshore_s{simpl}.geojson",
        regions_offshore=RESOURCES
        + "{interconnect}/Geospatial/regions_offshore_s{simpl}.geojson",
        custom_busmap=(
            DATA + "{interconnect}/custom_busmap_{clusters_hires}.csv"
            if config["enable"].get("custom_busmap", False)
            else []
        ),
        tech_costs=RESOURCES
        + f"costs/costs_{config['scenario']['planning_horizons'][0]}.csv",
        itl_reeds_zone="repo_data/ReEDS_Constraints/transmission/transmission_capacity_init_AC_ba_NARIS2024.csv",
        itl_county="repo_data/ReEDS_Constraints/transmission/transmission_capacity_init_AC_county_NARIS2024.csv",
        itl_trans_grp="repo_data/ReEDS_Constraints/transmission/transmission_capacity_init_AC_transgrp_NARIS2024.csv",
        itl_costs_reeds_zone="repo_data/ReEDS_Constraints/transmission/transmission_distance_cost_500kVdc_ba.csv",
        itl_costs_county="repo_data/ReEDS_Constraints/transmission/transmission_distance_cost_500kVac_county.csv",
    output:
        network=RESOURCES + "{interconnect}/elec_s{simpl}_ch{clusters_hires}.nc",
        regions_onshore=RESOURCES
        + "{interconnect}/Geospatial/regions_onshore_s{simpl}_ch{clusters_hires}.geojson",
        regions_offshore=RESOURCES
        + "{interconnect}/Geospatial/regions_offshore_s{simpl}_ch{clusters_hires}.geojson",
        busmap=RESOURCES + "{interconnect}/busmap_s{simpl}_ch{clusters_hires}.csv",
        linemap=RESOURCES + "{interconnect}/linemap_s{simpl}_ch{clusters_hires}.csv",
    log:
        "logs/cluster_network/{interconnect}/elec_s{simpl}_ch{clusters_hires}.log",
    benchmark:
        "benchmarks/cluster_network/{interconnect}/elec_s{simpl}_ch{clusters_hires}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: (input.size // 100000) * attempt * 2,
    script:
        "../scripts/cluster_network.py"


rule add_extra_components_hires:
    input:
        **{
            f"phs_shp_{hour}": "repo_data/"
            + f"psh/40-100-dam-height-{hour}hr-no-croplands-no-ephemeral-no-highways.gpkg"
            for phs_tech in config["electricity"]["extendable_carriers"]["StorageUnit"]
            if "PHS" in phs_tech
            for hour in phs_tech.split("hr_")
            if hour.isdigit()
        },
        network=RESOURCES + "{interconnect}/elec_s{simpl}_ch{clusters_hires}.nc",
        tech_costs=lambda wildcards: expand(
            RESOURCES + "costs/costs_{year}.csv",
            year=config["scenario"]["planning_horizons"],
        ),
        regions_onshore=RESOURCES
        + "{interconnect}/Geospatial/regions_onshore_s{simpl}_{clusters_hires}.geojson",
    params:
        retirement=config["electricity"].get("retirement", "technical"),
        demand_response=config["electricity"].get("demand_response", {}),
    output:
        RESOURCES + "{interconnect}/elec_s{simpl}_ch{clusters_hires}_ec.nc",
    log:
        "logs/add_extra_components/{interconnect}/elec_s{simpl}_ch{clusters_hires}_ec.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: (input.size // 100000) * attempt * 2,
    group:
        "prepare"
    script:
        "../scripts/add_extra_components.py"


rule prepare_network_hires:
    params:
        time_resolution=config_provider("clustering", "temporal", "resolution_elec"),
        adjustments=False,
        links=config_provider("links"),
        lines=config_provider("lines"),
        co2base=config_provider("electricity", "co2base"),
        co2limit=config_provider("electricity", "co2limit"),
        co2limit_enable=config_provider("electricity", "co2limit_enable", default=False),
        gaslimit=config_provider("electricity", "gaslimit"),
        gaslimit_enable=config_provider("electricity", "gaslimit_enable", default=False),
        transmission_network=config_provider("model_topology", "transmission_network"),
        costs=config_provider("costs"),
        autarky=config_provider("electricity", "autarky"),
    input:
        network=(
            config["custom_files"]["files_path"]
            + config["custom_files"]["network_name"]
            if config["custom_files"].get("activate", False)
            else RESOURCES + "{interconnect}/elec_s{simpl}_ch{clusters_hires}_ec.nc"
        ),
        tech_costs=(
            config["custom_files"]["files_path"] + "costs_2030.csv"
            if config["custom_files"].get("activate", False)
            else RESOURCES
            + f"costs/costs_{config['scenario']['planning_horizons'][0]}.csv"
        ),
    output:
        RESOURCES + "{interconnect}/elec_s{simpl}_ch{clusters_hires}_ec_l{ll}_{opts}.nc",
    log:
        solver="logs/prepare_network/{interconnect}/elec_s{simpl}_ch{clusters_hires}_ec_l{ll}_{opts}.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: (input.size // 100000) * attempt * 2,
    group:
        "prepare"
    log:
        "logs/prepare_network",
    script:
        "../scripts/prepare_network.py"
