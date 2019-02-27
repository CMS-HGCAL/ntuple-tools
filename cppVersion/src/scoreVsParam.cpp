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


string testParam = "critical_distance_EE";
EParam paramNum = kCriticalDistanceEE;
double testParamStep = 0.5;

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
	
	for(double param = paramMin[paramNum]; param < paramMax[paramNum]; param += testParamStep){
		UpdateParamValue(configPath, testParam, param);
		
		double score = GetScore();
		scoreVsParamHist->SetBinContent(scoreVsParamHist->GetXaxis()->FindFixBin(param), score);
		cout<<"Score:"<<score<<endl;
	}
	
	// Save the initial value back in the config
	system(("rm "+configPath).c_str());
	
	scoreVsParamHist->SaveAs(("score_vs_"+testParam+".root").c_str());
	return 0;
}


