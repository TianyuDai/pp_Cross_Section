#include <fstream>
#include <math.h>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"
#include "fastjet/Selector.hh"
//using namespace Pythia8; 
std::ofstream fout("non"); 
int main()
{
	int nEvent=10000000, power=-1; // power=-1: anti-kT
	double jetRadius=0.4, jetrapMin=0.2, jetrapMax=0.8, delta_phi=2*M_PI, delta_rap=(jetrapMax-jetrapMin); 

	std::vector <double> jet_cs, err, sqSum, pTMean, pTMeanBin; 
	std::vector <int> jet_ct; 
	//std::vector <double> pTHatBin{1., 5., 20., 40., 60., 80., 90., 100.}; 
	std::vector <double> pTHatBin{0., 200.}; 
	std::vector <double> pTBin{5., 6.2, 7.6, 9.3, 11.4, 14.1, 17.3, 21.3, 26.2, 32.2, 39.6, 48.7}; 
	size_t nBin = pTHatBin.size() - 1; 
	
	jet_cs = std::vector<double>(pTBin.size()-1, 0.); 
	err = std::vector<double>(pTBin.size()-1, 0.); 
	sqSum = std::vector<double>(pTBin.size()-1, 0.); 
	pTMean = std::vector<double>(pTBin.size()-1, 0.); 

	Pythia8::Pythia pythia;
	Pythia8::Event& event = pythia.event; 
	Pythia8::Info& info = pythia.info; 
	Pythia8::Settings& settings = pythia.settings;

	// Settings of Pythia 8 wrapper program
	pythia.readString("Main:numberOfEvents = 100000"); 
	pythia.readString("Next:numberShowEvent = 0"); 

	// Random seed
	pythia.readString("Random:setSeed = on"); 
	pythia.readString("Random:seed = 0");

	// Beam parameter settings 
	pythia.readString("Beams:idA = 2212"); 
	pythia.readString("Beams:idB = 2212"); 
	pythia.readString("Beams:eCM = 200."); 

	// Process setup
	//pythia.readString("HardQCD:all = on"); 
	//pythia.readString("SoftQCD:all = on"); 
	pythia.readString("SoftQCD:nonDiffractive = on"); 
	//pythia.readString("SoftQCD:singleDiffractive = on"); 
	//pythia.readString("SoftQCD:doubleDiffractive = on"); 

	// Set cuts
	// Use this for hard leading-jets in a certain mHat window
	settings.parm("PhaseSpace:mHatMin", 0.); 
	settings.parm("PhaseSpace:mHatMaz", 200.); 

	// Makes particles with c*tau > 10 mm stable
	settings.parm("ParticleDecays:limitTau0 = on"); 
	settings.parm("ParticleDecays:tau0Max = 10.0"); 
 
	fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, jetRadius, power);
 	std::vector <fastjet::PseudoJet> fjInputs; 
	fastjet::Selector select_maxrapidity = fastjet::SelectorRapMax(jetrapMax); 
	fastjet::Selector select_minrapidity = fastjet::SelectorRapMin(jetrapMin); 
	fastjet::Selector select_rapidity = select_minrapidity && select_maxrapidity; 
	
	std::vector<int> pT(nBin, 0); 
	double jetpTMin = pTBin[0]; 
	for (size_t iBin=0; iBin<nBin; iBin++)
	{
		settings.parm("PhaseSpace:pTHatMin", pTHatBin[iBin]); 
		settings.parm("PhaseSpace:pTHatMax", pTHatBin[iBin+1]); 
		pythia.init(); 
		jet_ct = std::vector<int>(pTBin.size()-1, 0); 
		pTMeanBin = std::vector<double>(pTBin.size()-1, 0.); 
		for (int iEvent=0; iEvent<nEvent; ++iEvent)
		{
			if (!pythia.next()) continue; 
			fjInputs.resize(0); 
			for (int i=0; i<event.size(); ++i)
			{
				if (event[i].isFinal() && event[i].isVisible())
				{
					fastjet::PseudoJet particleTemp = event[i]; 
					fjInputs.push_back(particleTemp); 
				}
			}
			Pythia8::vector <fastjet::PseudoJet> inclusiveJets, sortedJets; 
    			fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    			inclusiveJets = clustSeq.inclusive_jets(jetpTMin); 
			std::vector <fastjet::PseudoJet> selected_jets = select_rapidity(inclusiveJets); 
    			sortedJets    = sorted_by_pt(selected_jets); 
			for (unsigned int i=0; i < sortedJets.size(); i++)
				for (unsigned int j=0; j<pTBin.size()-1; j++)
					if (sortedJets[i].perp()>pTBin[j] && sortedJets[i].perp()<pTBin[j+1])
					{
						jet_ct[j]++; 
						pTMeanBin[j] += sortedJets[i].perp(); 
						break; 
					}
		}
		double sigmapb_weight = info.sigmaGen() * 1.0e9 / nEvent; 
		for (unsigned int i=0; i<pTBin.size()-1; i++)
		{
			jet_cs[i] += jet_ct[i]*sigmapb_weight; 
			sqSum[i] += jet_ct[i]*pow(sigmapb_weight, 2); 
			if (jet_ct[i]) pTMean[i] += pTMeanBin[i]/jet_ct[i]; 
		}
	}
	for (unsigned int i=0; i<pTBin.size()-1; i++)
	{
		err[i] = jet_cs[i]/sqrt(pow(jet_cs[i], 2)/sqSum[i]); 
		fout<<(pTBin[i]+pTBin[i+1])/2<<" "<<jet_cs[i]/(pTBin[i+1]-pTBin[i])/delta_phi/delta_rap<<" "<<err[i]/(pTBin[i+1]-pTBin[i])/delta_phi/delta_rap<<" "<<std::endl; 
	}
	return 0; 
}
