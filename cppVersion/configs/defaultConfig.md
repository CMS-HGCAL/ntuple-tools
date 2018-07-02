**General settings**

verbosity_level:  0

### range of ntuples to test (will be appended to the inputPath string below)
min_Ntuple:  10
max_Ntuple:  11

### stop analyzing each ntuple after that many events: 
analyze_events_per_tuple:     100

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

### cut on energy, with sensor dependance (relative to the noise):
energy_min_with_sensor_dependance:  3.0

### cut on energy, without sensor dependance (absolute units):
energy_min_no_sensor_dependance:    0.060

### kappa, with sensor dependance:
kappa_with_sensor_dependance:  9.0

### kappa, without sensor dependance:
kappa_no_sensor_dependance:    10.0

### test only within this layers range:
min_layer: 0
max_layer: 40




