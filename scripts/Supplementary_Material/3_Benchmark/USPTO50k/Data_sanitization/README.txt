USPTO-50k Dataset Preparation and Sanitization
This folder contains code to process and sanitize the USPTO-50k dataset for benchmarking purposes.

Steps:
1. Download the data
   Obtain the dataset from the Supporting Information (SI) of the following publication:
   https://pubs.acs.org/doi/abs/10.1021/acs.jcim.6b00564
   Place the file ci6b00564_si_002.zip into the working directory.

2. Generate the dataset
   Run Prepare_dataset.py to automatically extract and generate dataset_subset.csv by pulling reactions from the SI file.

3. Sanitize the dataset
   The following scripts perform sanitization and generate cleaned datasets for benchmarking:
      Sanitization_onlyCHO.py
      Sanitization_onlyN.py
      Sanitization_onlyS.py