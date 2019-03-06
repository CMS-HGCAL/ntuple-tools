**General settings**

verbosity_level:  0

### range of ntuples to test (will be appended to the inputPath string below)
min_Ntuple:  10
max_Ntuple:  11

### stop analyzing each ntuple after that many events: 
analyze_events_per_tuple:     1000

### input and output paths

input_path: ../../data/_SingleGammaPt100Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO/NTUP/_SingleGammaPt100Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO_NTUP_

output_path: ../../clusteringResultsCXX

**HGCal Imaging Algo parameters**

### should detector thickness be taken into account?
depend_sensor:  0

### Energy density inclusion funciton. Defines a way to include or reject hit when calculating local energy density.
### Possible options:
### *step* - include 100% of hit energy if distance smaller than *critical_distance*  
### *gaus* - include fraction of hit's energy based on the distance scaled by a gaussian distribution with width *critical_distribution*
energy_density_function:  step

### Critical distance for energy density ρ calculation (in cartesian coordiantes in cm, separately for each detector)
### Hits that are further than d_c from given hit will not be included in the energy density calculation for this hit.
critical_distance_EE:  2.0
critical_distance_FH:  2.0
critical_distance_BH:  2.0

### Maximum distance to the closest hit in a cluster to assign this hit to that cluster (in cm)
assignment_distance_EE:    2.0
assignment_distance_FH:    2.0
assignment_distance_BH:    2.0

### Critical distance to higher ρ hit (in cartesian coordiantes in cm, separately for each detector)
### Hits that are further than δ_c from any hit with higher ρ will be considered as potential cluster seeds.
deltac_EE:    2.0
deltac_FH:    2.0
deltac_BH:    5.0

### Critical energy density is defined as ρ_c = max(ρ)/κ. Hits with ρ > ρ_c will be considered as potential cluster seeds.
kappa:    10.0

### cut on energy (absolute units):
energy_min:    0.060

### test only within this layers range:
min_layer: 0
max_layer: 40

### Include only events were all particles reached EE before converting
reachedEE_only: 1

### Turn on/off the logic of marking hits as halo and removing them from a cluster
do_halo: 0

**Cluster Matching parameters**

### Sim clusters further than matching_max_distance from the clostest rec cluster will not be matched 
matching_max_distance: 5.0


