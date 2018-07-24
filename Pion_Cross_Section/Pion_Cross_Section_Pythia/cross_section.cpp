#include <fstream>
#include <math.h>
#include "Pythia8/Pythia.h"
using namespace Pythia8; 
ofstream fout("cs8_10000"); 
int main()
{
	int nEvent = 10000; 
	int N = 9, nBin = 8; //N is the number of pT bins in the plot
	double delta_pT = 2.0, delta_y = 0.35, delta_phi = 2*M_PI, pTBin[nBin+1] = {1., 3., 5., 10., 20., 40., 60., 80.}, cs[N] = {0.}, pTAvg[N], err[N]={0.}, sqSum[N] = {0.}; 
	Pythia pythia; 
	Settings& settings = pythia.settings; 
	Info& info = pythia.info; 
	pythia.readString("HardQCD:all = on"); 
	pythia.readString("Beams:eCM = 200."); 
	for (size_t iBin=0; iBin<nBin; ++iBin)
	{
		settings.parm("PhaseSpace:pTHatMin", pTBin[iBin]); 
		settings.parm("PhaseSpace:pTHatMax", pTBin[iBin+1]); 
		pythia.init(); 
		double dist[N] = {0.}; 
		double pTSum[N] = {0.}; 
		//double sigmaNorm; 
		//sigmaNorm = info.sigmaGen() / nEvent / delta_pT / (2*delta_y) / delta_phi / 2; 
		for (int iEvent=0; iEvent<nEvent; ++iEvent)
		{
			if (!pythia.next()) {continue; iEvent--; }	
			double weight = info.weight(); 
			double dist_[N]; 
			for (size_t j=0; j<N; j++) 
				dist_[j]= dist[j]; 
			for (int i=0; i<pythia.event.size(); ++i)
				if (pythia.event[i].isFinal() && (abs(pythia.event[i].id())==211) && (abs(pythia.event[i].y())<delta_y) && (pythia.event[i].pT()>4.0) && (pythia.event[i].pT()<40.0)) 
				{
					dist[int(pythia.event[i].pT()/2-2)]++; 
					pTSum[int(pythia.event[i].pT()/2-2)] += pythia.event[i].pT(); 
				}
			for (size_t j=0; j<N; j++)
			{
				dist[j] /= weight;
				//sqSum[j] += pow(dist[j]-dist_[j], 2) * pow(sigmaNorm*nEvent, 2); 
			}
		} 
		for (size_t j=0; j<N; j++)
		{
			pTAvg[j] = pTSum[j]/dist[j]; 
			//cs[j] += sigmaNorm*dist[j]; 
		}
	}
	for (size_t j=0; j<N; j++)
	{
		err[j] = sqrt((sqSum[j]-nEvent*nBin*cs[j])/(nEvent*nBin)/(nEvent*nBin-1)); 
		fout<<pTAvg[j]<<" "<<cs[j]/pTAvg[j]<<" "<<err[j]/pTAvg[j]<<endl;	//divided by 2 denotes (\Pi^++\Pi^-)/2, while divided by pTAvg means 1/pT*dN/dp^2
	}
}
