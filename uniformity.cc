// Include files
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2F.h"
#include "TBranch.h"
#include "TBranchElement.h"
//
#include "Cintex/Cintex.h"
//
#include "Event.hh"
#include "RunHeader.hh"
#include "TrackerHit.hh"
#include "CalorimeterHit.hh"
#include "DRCalorimeterHit.hh"
#include "DRTSCalorimeterHit.hh"
#include "PhotonHit.hh"
#include "TrackerHit.hh"
//
//
// standard C++ header files
//
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include<math.h>
using namespace std;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

int main(int argc, char** argv) {

    TSystem ts;
    gSystem->Load("libCintex");
    ROOT::Cintex::Cintex::Enable();
    gSystem->Load("../lib/libClassesDict");
    //  ROOT::Cintex::Cintex::SetDebug(2);

    if (argc < 3) {
        G4cout << "Program requires 2 argument: name of input file, name of output file" << G4endl;
        exit(1);
    }
    TFile* outfile = new TFile(argv[2], "RECREATE");
    TDirectory *histotop = outfile->mkdir("histos");
    histotop->cd();
    TH2F *h2 = new TH2F("uniformity", "uniformity", 100, -25., 25., 100, -25., 25.);
    string fname(argv[1]);
    ifstream infile;
    string line;
    infile.open(fname.c_str());
    if (infile.is_open()) {
        while (infile.good()) {
            getline(infile, line);
            //           cout << line << endl;
            if (string::npos != line.find("root")) {


                TFile fo(line.c_str());
                RunHeader *runh = new RunHeader();
                Event *event = new Event();
                TTree *T = (TTree*) fo.Get("T");
                T->SetBranchAddress("event.", &event);
                T->SetBranchAddress("RunHeader.", &runh);
                TBranch* frunhbranch = T->GetBranch("RunHeader.");
                TBranch* fevtbranch = T->GetBranch("event.");
                frunhbranch->GetEntry(0);
//                runh->Print();
                Int_t nevent = fevtbranch->GetEntries();
                //            G4cout << " Nr. of Events:  " << nevent << G4endl;
                double nphotons = 0.0;
                for (Int_t i = 0; i < nevent; i++) {
                    fevtbranch->GetEntry(i);

                    //-------------------------------------------
                    //
                    // Now we deal with the Hit Collections.
                    //
                    //-------------------------------------------
                    std::map<G4String, std::vector<G4VHit*> >* hcmap = event->GetHCMap();
                    std::map<G4String, std::vector<G4VHit*> >::iterator hciter;
                    for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
                        std::vector<G4VHit*> hits = (*hciter).second;
                        G4int NbHits = hits.size();
                        std::vector<std::string> y = split((*hciter).first, '_');
                        std::string Classname = y[1];
                        if (Classname == "PhotonDetector") {
                            nphotons = nphotons + double(NbHits);
                        }
                    }
                } // end loop over events
                G4ThreeVector partpos = runh->GetParticlePosition();
                //G4cout << " nr of photons per event  " << nphotons / double(nevent) << G4endl;
 //               G4cout << partpos.getY()<<"  "<< partpos.getZ()<<G4endl;
                h2->Fill(partpos.getY(), partpos.getZ(), nphotons / double(nevent));
            }
        }
    }
    G4cout << "===========================================" << G4endl;
    G4cout << " nr of bytes written:  " << outfile->Write() << G4endl;
    G4cout << "===========================================" << G4endl;
}
