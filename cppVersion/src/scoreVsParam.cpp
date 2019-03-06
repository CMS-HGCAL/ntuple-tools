//
//  scoreVsParam.cpp
//
//  Created by Jeremi Niedziela on 27/02/2019.
//

#include "GeneticHelpers.hpp"
#include "Helpers.hpp"
#include "ConfigurationManager.hpp"

#include <iostream>

#include <TROOT.h>
#include <TH1D.h>

using namespace std;

const double severityFactor = 10.0;
const string baseConfigPath = "testScore.md";
string configPath;

double GetScore()
{
	string command = "./createQualityPlots " + configPath;
	system(command.c_str());
	
	auto config = my::make_unique<ConfigurationManager>(configPath);
	ClusteringOutput output = ReadOutput(config->GetScoreOutputPath());
	
	double distance = fabs(output.containmentMean-1)
									+      output.containmentSigma
									+ fabs(output.resolutionMean)
									+      output.resolutionSigma
									+      output.separationMean
									+      output.separationSigma
									+ fabs(output.deltaNclustersMean)
									+      output.deltaNclustersSigma
									+      output.nRecoFailed
									+      output.nCantMatchRecSim
									+      output.nFakeRec;
	
	double score = severityFactor/distance;
	
	if(   output.resolutionMean     > 1000
		 || output.separationMean     > 1000
		 || output.containmentMean    > 1000
		 || output.deltaNclustersMean > 1000
		 || output.nRecoFailed        > 1000
		 || output.nCantMatchRecSim   > 1000
		 || output.nFakeRec           > 1000
		 )
	{ // this means that clustering failed completely
		score = 0;
	}
	if(score < 1E-5){ // just round down to zero if score it extremaly poor
		score = 0;
	}
	
	if(score==0){
		cout<<"This configuration failed completely"<<endl;
	}
	return score;
}

int main(int argc, char* argv[])
{
	if(argc != 3){
		cout<<"Usage: param_name step"<<endl;
		cout<<"\nwhere param name can be:"<<endl;
		cout<<"\tcritical_distance_EE\n\tcritical_distance_FH\n\tcritical_distance_BH"<<endl;
    cout<<"\tassignment_distance_EE\n\tassignment_distance_FH\n\tassignment_distance_BH"<<endl;
		cout<<"\tdeltac_EE\n\tdeltac_FH\n\tdeltac_BH\n\tkappa"<<endl;
		exit(0);
	}
	string testParam = argv[1];
	double testParamStep = stod(argv[2]);
	
	EParam paramNum;
	if(			testParam == "critical_distance_EE") 	  paramNum = kCriticalDistanceEE;
	else if(testParam == "critical_distance_FH") 	  paramNum = kCriticalDistanceFH;
	else if(testParam == "critical_distance_BH") 	  paramNum = kCriticalDistanceBH;
  else if(testParam == "assignment_distance_EE")  paramNum = kAssignmentDistanceEE;
  else if(testParam == "assignment_distance_FH")  paramNum = kAssignmentDistanceFH;
  else if(testParam == "assignment_distance_BH")  paramNum = kAssignmentDistanceBH;
	else if(testParam == "deltac_EE") 						  paramNum = kDeltacEE;
	else if(testParam == "deltac_FH") 						  paramNum = kDeltacFH;
	else if(testParam == "deltac_BH") 						  paramNum = kDeltacBH;
	else if(testParam == "kappa") 								  paramNum = kKappa;
	else{
		cout<<"Passed unknown param "<<testParam<<endl;
		exit(1);
	}
	
	cout<<"Testing "<<testParam<<" with title "<<paramTitle[paramNum]<<" and step "<<testParamStep<<endl;
	
	cout<<"Starting score vs. param"<<endl;
	gROOT->ProcessLine(".L loader.C+");
	
	auto scoreVsParamHist = my::make_unique<TH1D>(("score vs "+testParam).c_str(),
																								("score vs "+testParam).c_str(),
																								(paramMax[paramNum]-paramMin[paramNum])/testParamStep+1,
																								paramMin[paramNum], paramMax[paramNum]);
	
	// Copy default config with some random name
	srand((uint)time(NULL));
	int r = RandInt(100, 10000000);
	configPath = "testScore_"+to_string(r)+".md";
	system(("cp "+baseConfigPath+" "+configPath).c_str());
	
	// Put random score output path in the random config
	string scoreOutputPath = "testScore_"+to_string(r)+".out";
	UpdateParamValue(configPath, "score_output_path", scoreOutputPath);
	
	for(double param = paramMin[paramNum]; param < paramMax[paramNum]; param += testParamStep){
		UpdateParamValue(configPath, testParam, param);
		
		double score = GetScore();
		scoreVsParamHist->SetBinContent(scoreVsParamHist->GetXaxis()->FindFixBin(param), score);
		cout<<"Score:"<<score<<endl;
	}
	
	// Remove random files
	system(("rm "+configPath).c_str());
	system(("rm "+scoreOutputPath).c_str());
	
	scoreVsParamHist->SaveAs(("score_vs_"+testParam+".root").c_str());
	return 0;
}


