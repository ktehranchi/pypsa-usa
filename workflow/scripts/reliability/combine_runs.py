import os
import pandas as pd


# Specify the name of the output CSV file
output_file = snakemake.output[0]

# Initialize an empty list to store DataFrames
data_frames = []

# Loop through all CSV files in the directory
for filename in snakemake.input:
    data = pd.read_csv(filename, header=None)
    data_frames.append(data)

# Concatenate all DataFrames in the list into one
combined_data = pd.concat(data_frames, ignore_index=True)

combined_data.rename(columns={0: 'iteration'}, inplace=True)

# Sort the combined data by the iteration column
combined_data = combined_data.sort_values(by=['iteration'], ascending=True)
# combined_data.set_index('iteration',inplace=True)

# df_melted = pd.melt(combined_data, id_vars='iteration', value_vars=combined_data.columns)

# Write the combined data to a new CSV file
combined_data.to_csv(output_file, index=False)

print(f'Combined data saved to {output_file}')


summary_table = combined_data.describe()
summary_table = summary_table.loc[['mean', 'std']]

# Calculate relative error in percent
relative_error = (summary_table.loc['std'] / summary_table.loc['mean']) * 100

# Add relative error to the summary table
summary_table.loc['relative_error'] = relative_error

summary_table.to_csv(output_file.replace('reliability_results.csv', 'reliability_summary.csv'), index=True)

