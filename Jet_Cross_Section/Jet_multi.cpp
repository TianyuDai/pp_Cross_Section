#include <fstream>
#include <math.h>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"
#include "fastjet/Selector.hh"
using namespace Pythia8; 
ofstream fout("jet_multiplicity_err"); 
int main()
{
	int nEvent=50000, power=-1; // power=-1: anti-kT
	double jetRadius=0.4, jetpTMin=60., jetyMax=2.8; 

	std::vector <int> multiBin{2, 3, 4, 5, 6}; 
	std::vector <double> cs(multiBin.size(), 0.); 
	std::vector <double> pTHatBin{28., 50., 100., 200., 500., 1000., 1500., 2000., 2500., 3000., 7000.}; 
	std::vector <double> err(multiBin.size(), 0.); 

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
	settings.parm("PhaseSpace:mHatMaz", 7000.); 

	// Makes particles with c*tau > 10 mm stable
	settings.parm("ParticleDecays:limitTau0 = on"); 
	settings.parm("ParticleDecays:tau0Max = 10.0"); 
 
	//pythia.init(); 
	fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, jetRadius, power);
 	std::vector <fastjet::PseudoJet> fjInputs; 
	fastjet::Selector select_rapidity = fastjet::SelectorAbsRapMax(jetyMax); 

	for (unsigned int iBin=0; iBin<pTHatBin.size()-1; iBin++)
	{
		settings.parm("PhaseSpace:pTHatMin", pTHatBin[iBin]); 
		settings.parm("PhaseSpace:pTHatMax", pTHatBin[iBin+1]); 
		pythia.init(); 
		std::vector <int> count(multiBin.size(), 0); 
		//std::vector <double> ierror(multiBin.size(), 0.); 
		for (int iEvent=0; iEvent<nEvent; ++iEvent)
		{
			if (!pythia.next()) continue; 
			fjInputs.resize(0); 
			for (int i=0; i<event.size(); ++i)
				if (event[i].isFinal() && event[i].isVisible())
				{
					fastjet::PseudoJet particleTemp = event[i]; 
					fjInputs.push_back(particleTemp); 
				}
			vector <fastjet::PseudoJet> inclusiveJets, sortedJets; 
    			fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    			inclusiveJets = clustSeq.inclusive_jets(jetpTMin); 
			std::vector <fastjet::PseudoJet> selected_jets = select_rapidity(inclusiveJets); 
    			sortedJets    = sorted_by_pt(selected_jets); 
			if (sortedJets.size() >= 2)
				if (sortedJets[0].perp() >= 80.0)
				{
					count[0]++; 
					if (sortedJets.size() >= 3) count[1]++; 
					if (sortedJets.size() >= 4) count[2]++; 
					if (sortedJets.size() >= 5) count[3]++; 
					if (sortedJets.size() >= 6) count[4]++; 
				}
		}
		double sigmapb_weight = info.sigmaGen() * 1.0e9; 
		for (unsigned int i = 0; i < multiBin.size(); i++)
		{
			double ierror = sqrt(((double)count[i]/nEvent - pow(count[i]/nEvent, 2))/nEvent); 
			err[i] += pow(ierror * sigmapb_weight, 2); 
			cs[i] += count[i] * sigmapb_weight / nEvent; 
		}
	}
	for (unsigned int i=0; i < multiBin.size(); i++)
	{
		err[i] = sqrt(err[i]); 
		fout << cs[i] << " " << err[i] << std::endl; 
	}
	return 0; 
}
