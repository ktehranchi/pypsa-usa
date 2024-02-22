# Rules to Optimize/Solve Network


rule solve_network:
    params:
        solving=config["solving"],
        foresight=config["foresight"],
        planning_horizons=config["scenario"]["planning_horizons"],
        co2_sequestration_potential=config["sector"].get(
            "co2_sequestration_potential", 200
        ),
    input:
        network=RESOURCES
        + "{interconnect}/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}.nc",
        config=RESULTS + "config.yaml",
    output:
        network=RESULTS
        + "{interconnect}/networks/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}.nc",
    log:
        solver=normpath(
            LOGS
            + "solve_network/{interconnect}/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}_solver.log"
        ),
        python=LOGS
        + "solve_network/{interconnect}/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}_python.log",
    benchmark:
        (
            BENCHMARKS
            + "solve_network/{interconnect}/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}"
        )
    threads: 4
    resources:
        mem_mb=memory,
        walltime=config["solving"].get("walltime", "12:00:00"),
    shadow:
        "minimal"
    conda:
        "../envs/environment.yaml"
    envmodules:
        "python/3.9.0",
        "gurobi/10.0.1_py39"
    script:
        "../scripts/subworkflows/pypsa-eur/scripts/solve_network.py"


rule solve_network_operations:
    params:
        solving=config["solving"],
        foresight=config["foresight"],
        planning_horizons=config["scenario"]["planning_horizons"],
        co2_sequestration_potential=config["sector"].get(
            "co2_sequestration_potential", 200
        ),
    input:
        network=RESOURCES
        + "{interconnect}/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}.nc",
        config=RESULTS + "config.yaml",
    output:
        network=RESULTS
        + "{interconnect}/networks/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}_operations.nc",
    log:
        solver=normpath(
            LOGS
            + "solve_network/{interconnect}/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}_solver.log"
        ),
        python=LOGS
        + "solve_network/{interconnect}/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}_python.log",
    benchmark:
        (
            BENCHMARKS
            + "solve_network/{interconnect}/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}"
        )
    threads: 4
    resources:
        mem_mb=memory,
        walltime=config["solving"].get("walltime", "12:00:00"),
    shadow:
        "minimal"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_operations_network.py"
