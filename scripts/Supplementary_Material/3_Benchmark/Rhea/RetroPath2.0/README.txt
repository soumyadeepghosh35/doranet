Rhea_random_1000.json
Contains 1,000 randomly selected and filtered reactions from the Rhea database (retrieved April 2025).

source_rhea1000.csv
Lists the reactants from Rhea_random_1000.json. This file is used as input for the RetroPath2.0 program.

Benchmark_retropath.py
Runs RetroPath2.0 on each reactant and aggregates all predictions into a single output file: combined_result.csv.

Check_reproduction_rate.py
Compares the predictions in combined_result.csv with the original Rhea subset and prints the reproduction rate.