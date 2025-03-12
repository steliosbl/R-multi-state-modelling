import pandas as pd
import glob
import os
from tqdm import tqdm

# Define the path to the directory containing the CSV files
data_path = "data/results/"

# Get a list of all CSV files in the directory
csv_files = glob.glob(os.path.join(data_path, "results_*.csv"))

# Initialize an empty list to store individual dataframes
dataframes = []

# Loop through the list of CSV files with a progress bar
for file in tqdm(csv_files, desc="Loading CSV files"):
    # Extract the number from the filename
    filename = os.path.basename(file)
    test = filename.split("_")[1]

    # Read the CSV file into a dataframe
    df = pd.read_csv(file)

    # Append the dataframe to the list
    df = df.rename(
        columns={
            "parameters.n_elements": "parameters.elements",
            "parameters.n_cpu": "parameters.cpus",
            "parameters.n_rep": "parameters.rep",
        }
    )

    # Remove "parameters." from all column names that contain it
    df.columns = df.columns.str.replace("parameters.", "")

    dataframes.append(df)

# Concatenate all dataframes into a single dataframe
combined_df = pd.concat(dataframes, ignore_index=True)

combined_df.to_csv("data/results/combined_results.csv", index=False)
