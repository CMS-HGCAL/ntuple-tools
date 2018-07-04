**General settings**

verbosity_level:  0

### range of ntuples to test (will be appended to the inputPath string below)
min_Ntuple:  10
max_Ntuple:  11

### stop analyzing each ntuple after that many events: 
analyze_events_per_tuple:     1000

### input and output paths

input_path: ../../data/_SingleGammaPt100Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO/NTUP/_SingleGammaPt100Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO_NTUP_

output_path: ../clusteringResultsCXX



**HGCal Imaging Algo parameters**

### should detector thickness be taken into account?
depend_sensor:  1

### ??? in cartesian coordiantes in cm, per detector (*clarification needed*)
deltac_EE:    2.0
deltac_FH:    2.0
deltac_BH:    5.0

### Request at least minClusters+1 2D clusters  (*clarification needed*)
min_clusters: 3

### cut on energy (relative to the noise):
energy_min:  3.0

### kappa (energy density threshold):
kappa:  9.0

### test only within this layers range:
min_layer: 0
max_layer: 40
