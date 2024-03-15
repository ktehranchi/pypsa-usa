import pandas as pd
import os
import sys
import pypsa

if __name__ == "__main__":
    n = snakemake.input.network
    samples = pd.read_csv(snakemake.input.copula_runs)
    samples.to_csv(snakemake.output[0], index=False)