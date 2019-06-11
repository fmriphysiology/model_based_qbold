model-based-qbold:

MATLAB scripts used to generate and analyse simulated quantitative BOLD data, and to anaylse the results from FABBER Bayesian analysis of in vivo data. See model-based qBOLD paper for more information. 

ASE_Simulation:

These scripts are used to generate and analyse simulated ASE qBOLD data from a single model voxel. Simulat_qASE.m creates a .mat file of simulated ASE data (an example is included as data_qASE_example.mat), which can be analysed in a 2D grid search using gridSearchBayesian.m. The implementation of the qBOLD model is in qASE_model.m, it uses the static dephasing regime in the tissue compartment (Yablonskiy & Haacke, 1994, and He & Yablonskiy, 2007), the motional narrowing model in the blood compartment (Berman & Pike, 2017). 

Fabber_Analysis:

These scripts analyse in vivo GASE data and are designed to be used with the data provided  in:

    Cherukara MT, Stone AJ, Chappell MA, Blockley NP. Data acquired to demonstrate model-based Bayesian inference of brain oxygenation using quantitative BOLD, Oxford University Research Archive 2018. doi: <Please see ORA entry for DOI> 

Which is structured in BIDS format. Linear_sqBOLD.m calculates R2’ and DBV using a log-linear single-compartment qBOLD model, while Fabber_Averages is used to analyse data produced by FABBER variational Bayesian inference using the qBOLD model (included in the above dataset). These require the following scripts to run (not included): read_avw.m, save_avw.m, sigstar.m 