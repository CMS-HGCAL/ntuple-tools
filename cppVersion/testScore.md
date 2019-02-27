**General settings**

verbosity_level:  0

### range of ntuples to test (will be appended to the inputPath string below)
min_Ntuple:  1
max_Ntuple:  1

### stop analyzing each ntuple after that many events: 
analyze_events_per_tuple:	50

### input and output paths

input_path:	../../data/eventQCD_wf24031p0_Pt_80_120_14TeV_2023D28_noPU_jniedzie_20190208/NTUP/eventQCD_x1320_1.0To35.0_NTUP_

output_path: ../clusteringResultsCXX/testScore/

score_output_path:	testScore.out

**HGCal Imaging Algo parameters**

### should detector thickness be taken into account?
depend_sensor:	1

### Energy density inclusion funciton. Defines a way to include or reject hit when calculating local energy density.
### Possible options:
### *step* - include 100% of hit energy if distance smaller than *critical_distance*  
### *gaus* - include fraction of hit's energy based on the distance scaled by a gaussian distribution with width *critical_distance*
### *exp* - exponential for distances smaller than *critical_distance*, don't include above

energy_density_function:	step


### Critical distance for energy density ρ calculation (in cartesian coordiantes in cm, separately for each detector)
### Hits that are further than d_c from given hit will not be included in the energy density calculation for this hit.
critical_distance_EE:	18.9
critical_distance_FH:	10.7
critical_distance_BH:	39.4

### Critical distance to higher ρ hit (in cartesian coordiantes in cm, separately for each detector)
### Hits that are further than δ_c from any hit with higher ρ will be considered as potential cluster seeds.
deltac_EE:	0.10
deltac_FH:	29.5
deltac_BH:	18.9

### Critical energy density is defined as ρ_c = κ * noiseLevel. Hits with ρ > ρ_c will be considered as potential cluster
kappa:	29.0

### cut on energy (relative to the noise):
energy_min:	3.0

### test only within this layers range:
min_layer:	1
max_layer:	53

### Include only events were all particles reached EE before converting
reachedEE_only: 0


**Cluster Matching parameters**

### Sim clusters further than matching_max_distance from the clostest rec cluster will not be matched 
matching_max_distance:	10.0

