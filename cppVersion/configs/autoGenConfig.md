**General settings**

verbosity_level:  0

### range of ntuples to test (will be appended to the inputPath string below)
min_Ntuple:  0
max_Ntuple:  0

### stop analyzing each ntuple after that many events: 
analyze_events_per_tuple:	10

### input and output paths

input_path: ../../../data/MultiParticleInConeGunProducer_test_jniedzie_20180704/NTUP/merged_NTUP_

output_path: ../../clusteringResultsCXX/optimizer/


**HGCal Imaging Algo parameters**

### should detector thickness be taken into account?
depend_sensor:	0

### Energy density inclusion funciton. Defines a way to include or reject hit when calculating local energy density.
### Possible options:
### *step* - include 100% of hit energy if distance smaller than *critical_distance*  
### *gaus* - include fraction of hit's energy based on the distance scaled by a gaussian distribution with width *critical_distance*
### *exp* - exponential for distances smaller than *critical_distance*, don't include above
  energy_density_function:  gaus


### Critical distance for energy density ρ calculation (in cartesian coordiantes in cm, separately for each detector)
### Hits that are further than d_c from given hit will not be included in the energy density calculation for this hit.
critial_distance_EE:	1.97583
critial_distance_FH:	1.97583
critial_distance_BH:	1.97583

### Critical distance to higher ρ hit (in cartesian coordiantes in cm, separately for each detector)
### Hits that are further than δ_c from any hit with higher ρ will be considered as potential cluster seeds.
deltac_EE:	5.97579
deltac_FH:	11.9758
deltac_BH:	21.9758

### Critical energy density is defined as ρ_c = max(ρ)/κ. Hits with ρ > ρ_c will be considered as potential cluster seeds.
kappa:	9.97578

### cut on energy (relative to the noise):
energy_min:	0.06

### Request at least minClusters+1 2D clusters  (*clarification needed*)
min_clusters:	2.75905

### test only within this layers range:
min_layer: 0
max_layer: 40

### Include only events were all particles reached EE before converting
reachedEE_only:	1


**Cluster Matching parameters**

### Sim clusters further than matching_max_distance from the clostest rec cluster will not be matched 
matching_max_distance:	4.97579

