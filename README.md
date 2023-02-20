To investigate the impact of speed and directional persistence of breast cancer cells on tumor growth, we use a 2D stochastic cell-based model inspired by the persistence random walk model (38). However, besides cell migration, our model also accounts for other cellular processes including, proliferation, intravasation and cell death. 

The model implementation is described in https://www.biorxiv.org/content/10.1101/2022.04.25.489358v1 Also the necessary modeling code and data to reproduce the modeling results/figures in the publication are available in the "Simulation Outcome" folder 

The paper includes two models (M1 and M2):
The first model is available in the "Speed-Persistence Model (M1)" folder whereas the second model is available in the "Golgi_Fragmentation_Model (M2)" folder. In each of M1 simulation, which is characterized by certain speed and directional persistence, we account for the following outcomes: the number and distribution of cells in the 2D domain, the number of intravasating cells, the proliferation outcome and the mean time spent in the nourished zones. In this way we were able to estimate the optimal speed-persistence combination that leads to bigger tumors and more intravasation. We also use the model to study the pathological relevance of Golgi fragmentation in breast cancer (M2). Examplels of the simulation output for M1 and M1 are available in the "Examples of Simulation Output" folder with ReadMe files explaining each of the outputs. 

To run the M1 simulations, both files ("Modeling_backend" & "Run_model") should be saved in the same working directory. Load the "Modeling_backend" in R and open the "Run_model" file and execute the code.

To run the M2 simulations, both files ("Golgi_Fragmentation_model__backend.R" & "Golgi_Fragmentation_model__run.R") should be saved in the same working directory. Load the "Modeling_backend" in R and open the "Run_model" file and execute the code.

Computer simulations were performed mainly in HPC solutions provided by the Oslo Centre for Epidemiology and Biostatistics, at the University of Oslo and also on the Abel computer cluster, owned by the University of Oslo, and operated by the Department for Research Computing at USIT, the University of Oslo IT-department
