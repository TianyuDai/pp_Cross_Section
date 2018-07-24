// main08.cc is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates methods to emphasize generation at high pT.

#include "Pythia8/Pythia.h"
#include <fstream>
using namespace Pythia8;
ofstream fout("output08.out"); 

int main() {

  int nEvent = 100;
  int nRange = 100;
  Pythia pythia;
  Settings& settings = pythia.settings;
  Info& info = pythia.info;

  int nBin = 4;
  double pTlimit[5] = {0., 5., 10., 20., 40.};

  for (int iBin = 0; iBin < nBin; ++iBin) {
    pythia.readString("HardQCD:all = on");
    pythia.readString("SoftQCD:nonDiffractive = off");
    settings.parm("PhaseSpace:pTHatMin", pTlimit[iBin]);
    settings.parm("PhaseSpace:pTHatMax", pTlimit[iBin + 1]);

    pythia.readString("Beams:eCM = 200.");
    pythia.init();

    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;
      double pTHat  = info.pTHat();
      double weight = info.weight();
      //fout<<weight<<endl; 
      fout<<info.sigmaGen()<<" "<< info.sigmaErr()<<endl; 
    }
    fout<<endl; 

  // End of pT-bin loop.
  }

  // Done.
  return 0;
}
