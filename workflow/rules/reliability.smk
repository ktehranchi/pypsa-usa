if config["run"]["reliability"]:

    rule generate_scenarios:
        params:
            n_samples = config['reliability']['MC_draws'],
        input:
            network=RESOURCES
            + "{interconnect}/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}.nc",
            config=RESULTS + "config.yaml",
            copula_runs= "repo_data/reliability_hw/simulated_failure_times_theta_{theta}.csv"
        output:
            samples = RESOURCES + "{interconnect}/c{clusters}_ec_l{ll}_{opts}_{sector}/samples_s_{theta}.csv",
        threads: 1,
        resources:
            mem_mb = memory,
        log:
            "logs/generate_scenarios/{interconnect}/c{clusters}_ec_l{ll}_{opts}_{sector}_{theta}.log",
        script:
            "../scripts/reliability/generate_scenarios.py"

    rule run_operations_model:
        params:
            draw_number="{draw_number}",
        input: 
            network=RESOURCES
            + "{interconnect}/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}.nc",
            samples = RESOURCES + "{interconnect}/c{clusters}_ec_l{ll}_{opts}_{sector}/samples_s_{theta}.csv",
        output:
            operations_results_raw = RESOURCES + "{interconnect}/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}/runs/{theta}/run_{draw_number}.csv",
        threads: 8,
        resources:
            mem_mb = memory,
        log:
            "logs/run_operations_model/{interconnect}/c{clusters}_ec_l{ll}_{opts}_{sector}/{theta}/draw{draw_number}.log",
        shell:
            '''
            # ml python/3.9.0
            # ml gurobi/10.0.1_py39
            # python3.9 scripts/reliability/run_operations_model.py {input} {wildcards.draw_number} 24 {output}
            python3 scripts/reliability/run_operations_model.py {input} {wildcards.draw_number} 24 {output}
            '''


    iterations = config['reliability']['MC_draws']

    rule compute_reliability_stats:
        input:
            operations_results_raw = lambda wildcards: expand(RESOURCES + "{interconnect}/elec_s_{clusters}_ec_l{ll}_{opts}_{sector}/runs/{theta}/run_{draw_number}.csv",
                draw_number=range(0, iterations),
                interconnect=wildcards.interconnect,
                clusters=wildcards.clusters,
                ll=wildcards.ll,
                opts=wildcards.opts,
                sector=wildcards.sector,
                theta=wildcards.theta,
                ),
        output:
            operations_results =RESULTS + "{interconnect}/c{clusters}_ec_l{ll}_{opts}_{sector}/reliability_results_{theta}.csv",
        threads: 1,
        log:
            "logs/combine_operations_results/{interconnect}/c{clusters}_ec_l{ll}_{opts}_{sector}_{theta}.log",
        script:
            "../scripts/reliability/combine_runs.py"