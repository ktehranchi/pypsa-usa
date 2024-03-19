import numpy as np
import pypsa
import pandas as pd
import logging
import matplotlib.pyplot as plt
import time
import sys


def extract_metrics(n):
    expected_unserved_energy = n.generators_t.p.loc[:,n.generators.carrier=='load'].round(3)
    expected_unserved_energy = expected_unserved_energy.loc[:,~expected_unserved_energy.columns.str.contains('British')]
    expected_unserved_energy = expected_unserved_energy.sum(axis=1)
    return expected_unserved_energy
    

def solve_operations_model(n, solver_options, snapshots, load_shedding = True, solver_name='gurobi'):
    """Solve the operations model for the given network and options"""
    print("Solving operations model")
    set_all_extendable_to(n, False)
    if load_shedding: n.optimize.add_load_shedding(sign=1, marginal_cost=10000,suffix=' load')
    n.optimize(snapshots, solver_name=solver_name, solver_options=solver_options)
    eue = extract_metrics(n)
    print(f"Expected unserved energy: {eue.sum()}")
    return n, eue

def set_all_extendable_to(network, val):
    """Sets components to be extendable or not"""
    for component in network.iterate_components():
        if hasattr(component.df, "p_nom_extendable"):
            component.df.p_nom_extendable = val
        if hasattr(component.df, "s_nom_extendable"):
            component.df.s_nom_extendable = val

def update_network_data(n, sample):
    """Update the network data with the given sample
    """
    # Sample is a 1D array with the each entry storing the time to failure for each component
    # Assume each component experiences a drops to X% of generation capacity
    derate_pct = 0.0
    n.loads_t.p_set = n.loads_t.p_set * 1.6
    # Let value = # of hours to failure
    hours_to_failure = (sample * 2).astype(int)

    #randomly choose from the lines
    num_lines = 8
    # generators_ca = generators_ca.sample(num_gens, random_state=np_gen)
    #take largest p_nom
    lines = n.lines.nlargest(num_lines, 's_nom')
    hours_to_failure.index = lines.index

    # Derate the capacity of the generators
    p_max_pu = pd.DataFrame(index=n.snapshots, columns=lines.index)
    for gen in lines.index:
        rating = np.ones(len(n.snapshots))
        rating[hours_to_failure[gen]:] = derate_pct
        p_max_pu[gen] = rating

    n.lines_t.s_max_pu = n.lines_t.s_max_pu.join(p_max_pu)

    return n

if __name__ == "__main__":
    logger = logging.getLogger(__name__)

    # Load Data
    n = pypsa.Network(sys.argv[1])
    samples = pd.read_csv(sys.argv[2])
    sample_num = int(sys.argv[3])
    sample = samples.iloc[sample_num, :]
    n_hours = int(24*7)
    output_path = sys.argv[5]
    solver_options =   {
                        'threads': 8,
                        'method': 2, # barrier,
                        'crossover': 0,
                        'BarConvTol': 1.e-5,
                        'FeasibilityTol': 1.e-6,
                        'AggFill': 0,
                        'PreDual': 0,
                        'GURO_PAR_BARDENSETHRESH': 200}

    #Run Operations Model
    operations_results = pd.DataFrame(np.zeros(shape=(1,3)), columns=['iteration', 'production_cost','runtime'])

    start_time = time.perf_counter()
    n = update_network_data(n, sample)
    n, eue = solve_operations_model(n, solver_options, n.snapshots[:n_hours])


    end_time = time.perf_counter()
    print(f"Finished iteration {sample_num} in {end_time - start_time:0.4f} seconds")
    results = np.append(sample_num, eue.values)

    # pd.DataFrame(results).T.to_csv(output_path, index=False) 
    results.tofile(output_path, sep=",", format="%s")
    # operations_results.to_csv(
    #     output_path, 
    #     index=False)