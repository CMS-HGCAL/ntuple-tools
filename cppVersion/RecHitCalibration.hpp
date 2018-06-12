//
//  RecHitCalibration.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//
//  Simple copy of the python script RecHitCalibration.py:
//
//   class for obtaining details from the RecHit calibration used in CMSSW
//   reproduces number listed at https://twiki.cern.ch/twiki/pub/CMS/HGCALSimulationAndPerformance/rechit.txt"""

#ifndef RecHitCalibration_hpp
#define RecHitCalibration_hpp

class RecHitCalibration {
public:
  RecHitCalibration(){}
    
  double MeVperMIP(int layer,int thicknessIndex){
    if(layer > 40){// no thickness correction for BH
      return dEdX_weights[layer];
    }
    return dEdX_weights[layer]/thicknessCorrection[thicknessIndex];
  }
  double MIPperGeV(int layer,int thicknessIndex){
    return 1000./MeVperMIP(layer, thicknessIndex);
  }
  double sigmaNoiseMIP(int layer,int thicknessIndex){
    if(layer > 40){// for BH, sigmaNoiseMIP = noise_MIP
      return noise_MIP;
    }
    return fC_per_ele * nonAgedNoises[thicknessIndex]/fCPerMIP[thicknessIndex];
  }
  double sigmaNoiseMeV(int layer,int thicknessIndex){
    return sigmaNoiseMIP(layer, thicknessIndex) * MeVperMIP(layer, thicknessIndex);
  }
private:
  // https://github.com/cms-sw/cmssw/blob/CMSSW_9_3_X/RecoLocalCalo/HGCalRecProducers/python/HGCalRecHit_cfi.py#L5
  double dEdX_weights[53] = {
    0.0,   // there is no layer zero
    8.603,  // Mev
    8.0675,
    8.0675,
    8.0675,
    8.0675,
    8.0675,
    8.0675,
    8.0675,
    8.0675,
    8.9515,
    10.135,
    10.135,
    10.135,
    10.135,
    10.135,
    10.135,
    10.135,
    10.135,
    10.135,
    11.682,
    13.654,
    13.654,
    13.654,
    13.654,
    13.654,
    13.654,
    13.654,
    38.2005,
    55.0265,
    49.871,
    49.871,
    49.871,
    49.871,
    49.871,
    49.871,
    49.871,
    49.871,
    49.871,
    49.871,
    62.005,
    83.1675,
    92.196,
    92.196,
    92.196,
    92.196,
    92.196,
    92.196,
    92.196,
    92.196,
    92.196,
    92.196,
    46.098};
  
  // https://github.com/cms-sw/cmssw/blob/CMSSW_9_3_X/RecoLocalCalo/HGCalRecProducers/python/HGCalRecHit_cfi.py#L86
  double thicknessCorrection[3] = {1.132,1.092,1.084};  // 100, 200, 300 um
  
  // Base configurations for HGCal digitizers
  // https://github.com/cms-sw/cmssw/blob/CMSSW_9_3_X/SimCalorimetry/HGCalSimProducers/python/hgcalDigitizer_cfi.py#L5
  // self.eV_per_eh_pair = 3.62
  double fC_per_ele = 1.6020506e-4;
  double nonAgedNoises[3] = {2100.0, 2100.0, 1600.0};  // 100,200,300 um (in electrons)
  
  // https://github.com/cms-sw/cmssw/blob/CMSSW_9_3_X/RecoLocalCalo/HGCalRecProducers/python/HGCalUncalibRecHit_cfi.py#L25
  double fCPerMIP[3] = {1.25, 2.57, 3.88};  // 100um, 200um, 300um
  
  // https://github.com/cms-sw/cmssw/blob/CMSSW_9_3_X/SimCalorimetry/HGCalSimProducers/python/hgcalDigitizer_cfi.py#L127
  double noise_MIP = 1.0/7.0; //expectation based on latest SiPM performance
  
};

#endif


