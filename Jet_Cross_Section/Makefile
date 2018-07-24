# deal with pythia
PYTHIADIR   = /home/tianyu/lib/pythia8230
PYTHIAFLAGS = -I$(PYTHIADIR)/include -O2 -pedantic -W -Wall -Wshadow -fPIC
PYTHIALIBS  =  -L$(PYTHIADIR)/lib -Wl,-rpath,$(PYTHIADIR) -lpythia8 -ldl

# deal with fastjet
FASTJETDIR  = /home/tianyu/lib/FastJet3.3.0
FASTJETARGS = `$(FASTJETDIR)/bin/fastjet-config --cxxflags --libs`

# LHAPDF if you want to replace Pythia's dummy default
LHAPDFARGS =
#LHAPDFARGS = -L$(HOME)/local/lib -lLHAPDF -I$(HOME)/local/include

main72:
	g++ -o main72 main72.cc $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 
main71:
	g++ -o main71 main71.cc $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 
jet:
	g++ -o Jet Jet_cs.cpp $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 
jet_200:
	g++ -std=c++11 -o Jet_200 Jet_cs_200.cpp $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 
jet_multi: 
	g++ -std=c++11 -o Jet_multi Jet_multi.cpp $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 

