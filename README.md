# CMB_Pancreatic_Cancer_Article

This Mathematical Model for Pancreatic Cancer was developed for a Project of the Master of Computational and Mathematical Biology of the Aix-Marseill University.
All equations and parameters are based on the article:

Y. Louzoun, C. Xue, G. B. Lesinski, and A. Friedman. A mathematical model for
pancreatic cancer growth and treatments. Journal of Theoretical Biology, 351:74–82,
2014. ISSN 0022-5193. doi: https://doi.org/10.1016/j.jtbi.2014.02.028. URL https:
//www.sciencedirect.com/science/article/pii/S0022519314001039.

1. Lambda_c_Simulation:

   This file plots a graph to simulate the behavior of Cancer cells (C) for different values of Lambda_c.
   
4. Simulation_Drugs

   This File simulates the effects of different drugs treatments on the Cancer cells (C) over time:
   The simulation is plotted for four scenarios:
   (1) Untreated patient;
   (2) Patient after TGF_beta silencing treatment;
   (3) Patient after immune activation treatment;
   (4) Patient after both treatments.

   It also plots the graph foe the bahavior of Cancer cells, Stale Pancreatic cells, Ratio of Macrophages and T-Cells over time for an untreated patient.
   
3. Simulation_Mean_50_patients

   This file simulates a population of 50 patients with different values of lambda_c. After the model runs for all those patients, the mean Cancer cells value (C) is calculated. Finally, the graph of Cancer Cells over time is plotted.
