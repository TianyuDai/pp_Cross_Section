#include <fstream>
#include <math.h>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"
#include "fastjet/Selector.hh"
using namespace Pythia8; 
ofstream fout("kk"); 
int main()
{
	int nEvent=10000, power=-1; // power=-1: anti-kT
	double jetRadius=0.4, jetpTMin=60., jetyMax=2.8; 
	size_t nBin = 16; 

	std::vector < std::vector <double> > pTBin, jet_cs, err, sqSum; 
	std::vector < std::vector <int> > jet_ct; 
	static const double arr[] = {28., 40., 60., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 2000., 3000., 7000.}; 
	vector <double> pTHatBin(arr, arr+sizeof(arr)/sizeof(arr[0])); 
	static const double arr1[] = {80., 110., 160., 210., 260., 310., 400., 500., 600., 800.}; 
	vector <double> leading(arr1, arr1+sizeof(arr1)/sizeof(arr1[0])); 
	pTBin.push_back(leading); 
	static const double arr2[] = {60., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800.}; 
	vector <double> second(arr2, arr2+sizeof(arr2)/sizeof(arr2[0])); 
	pTBin.push_back(second);
	static const double arr3[] = {60., 80., 110., 160., 210., 260., 310., 400., 500.}; 
	vector <double> third(arr3, arr3+sizeof(arr3)/sizeof(arr3[0])); 
	pTBin.push_back(third); 
	static const double arr4[] = {60., 80., 110., 160., 210.}; 
	vector <double> forth(arr4, arr4+sizeof(arr4)/sizeof(arr4[0])); 
	pTBin.push_back(forth); 
	
	jet_cs.reserve(pTBin.size()); 
	err.reserve(pTBin.size()); 
	sqSum.reserve(pTBin.size()); 
	//sq.reserve(pTBin.size()); 
	jet_ct.reserve(pTBin.size()); 
	//jet_ctp.reserve(pTBin.size()); 
	for (unsigned int i=0; i<pTBin.size(); i++)
	{
		jet_cs[i] = std::vector<double>(pTBin[i].size()-1,0.); 
		err[i] = std::vector<double>(pTBin[i].size()-1,0.); 
		sqSum[i] = std::vector<double>(pTBin[i].size()-1,0.); 
	}

	Pythia pythia;
	Event& event = pythia.event; 
	Info& info = pythia.info; 
	Settings& settings = pythia.settings;

	// Settings of Pythia 8 wrapper program
	pythia.readString("Main:numberOfEvents = 100000"); 
	pythia.readString("Next:numberShowEvent = 0"); 

	// Random seed
	pythia.readString("Random:setSeed = on"); 
	pythia.readString("Random:seed = 0");

	// Beam parameter settings 
	pythia.readString("Beams:idA = 2212"); 
	pythia.readString("Beams:idB = 2212"); 
	pythia.readString("Beams:eCM = 7000."); 

	// Process setup
	pythia.readString("HardQCD:all = on"); 

	// Set cuts
	// Use this for hard leading-jets in a certain mHat window
	settings.parm("PhaseSpace:mHatMin", 0.); 
	settings.parm("PhaseSpace:mHatMax", 7000.); 

	// Makes particles with c*tau > 10 mm stable
	settings.parm("ParticleDecays:limitTau0 = on"); 
	settings.parm("ParticleDecays:tau0Max = 10.0"); 
 
	//pythia.init(); 
	fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, jetRadius, power);
 	std::vector <fastjet::PseudoJet> fjInputs; 
	fastjet::Selector select_rapidity = fastjet::SelectorAbsRapMax(jetyMax); 
	
	for (size_t iBin=0; iBin<nBin; iBin++)
	{
		settings.parm("PhaseSpace:pTHatMin", pTHatBin[iBin]); 
		settings.parm("PhaseSpace:pTHatMax", pTHatBin[iBin+1]); 
		pythia.init(); 
		for (unsigned int i=0; i<pTBin.size(); i++)
		{
			jet_ct[i] = std::vector<int>(pTBin[i].size()-1, 0); 
			//sq[i] = std::vector<int>(pTBin[i].size()-1, 0); 
		}
		for (int iEvent=0; iEvent<nEvent; ++iEvent)
		{
			if (!pythia.next()) continue; 
			fjInputs.resize(0); 
			//for (unsigned int i=0; i<pTBin.size(); i++)
				//jet_ctp[i] = std::vector<int>(pTBin[i].size()-1, 0); 
			for (int i=0; i<event.size(); ++i)
			{
				if (event[i].isFinal() && event[i].isVisible())
				{
					fastjet::PseudoJet particleTemp = event[i]; 
					fjInputs.push_back(particleTemp); 
				}
			}
			vector <fastjet::PseudoJet> inclusiveJets, sortedJets; 
    			fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    			inclusiveJets = clustSeq.inclusive_jets(jetpTMin); 
			std::vector <fastjet::PseudoJet> selected_jets = select_rapidity(inclusiveJets); 
    			sortedJets    = sorted_by_pt(selected_jets); 
			if (sortedJets.size() > 1)
				if (sortedJets[0].perp() >= pTBin[0][0])
				{
					for (unsigned int j=0; j<pTBin[0].size()-1; j++)
						if (sortedJets[0].perp()>pTBin[0][j] && sortedJets[0].perp()<pTBin[0][j+1])
						{
							jet_ct[0][j]++; 
							break; 
						}
					for (unsigned int i=1; i<pTBin.size(); i++)
						if (sortedJets.size() > i)
							for (unsigned int j=0; j<pTBin[i].size()-1; j++)
								if (sortedJets[i].perp()>pTBin[i][j] && sortedJets[i].perp()<pTBin[i][j+1])
								{
									jet_ct[i][j]++; 
									break; 
								}
				}
			//for (unsigned int i=0; i<pTBin.size(); i++)
				//for (unsigned int j=0; j<pTBin[i].size()-1; j++)
				//{
					//jet_ct[i][j] += jet_ctp[i][j]; 
					//sq[i][j] += pow(jet_ctp[i][j], 2); 
				//}
		}
		//double sigmapb = info.sigmaGen() * 1.0e9; 
		double sigmapb_weight = info.sigmaGen() * 1.0e9 / nEvent; 
		for (unsigned int i=0; i<pTBin.size(); i++)
			for (unsigned int j=0; j<pTBin[i].size()-1; j++)
			{
				jet_cs[i][j] += jet_ct[i][j]*sigmapb_weight; 
				sqSum[i][j] += jet_ct[i][j]*pow(sigmapb_weight, 2); 
			}
	}
	for (unsigned int i=0; i<pTBin.size(); i++)
	{
		fout<<i+1<<endl; 
		for (unsigned int j=0; j<pTBin[i].size()-1; j++)
		{
			err[i][j] = jet_cs[i][j]/sqrt(pow(jet_cs[i][j], 2)/sqSum[i][j]); 
			//err[i][j] = sqrt((sqSum[i][j]/nEvent-pow(jet_cs[i][j], 2))/(nEvent-1)); 
			fout<<(pTBin[i][j]+pTBin[i][j+1])/2<<" "<<jet_cs[i][j]/(pTBin[i][j+1]-pTBin[i][j])<<" "<<err[i][j]/(pTBin[i][j+1]-pTBin[i][j])<<endl; 
		}
		fout<<endl; 
	}
	return 0; 
}
