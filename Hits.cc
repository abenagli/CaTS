// Include files
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
//
#include "Cintex/Cintex.h"
//
#include "Event.hh"
#include "TrackerHit.hh"
#include "CalorimeterHit.hh"
#include "DRCalorimeterHit.hh"
#include "PhotonHit.hh"
#include "TrackerHit.hh"
//
#include <vector>
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
    gSystem->Load("libClassesDict");
    ROOT::Cintex::Cintex::Enable();
    if (argc < 3) {
        G4cout << "Program requires 2 arguments:" << G4endl;
        G4cout << "Name of input file:" << G4endl;
        G4cout << "Name of output file:" << G4endl;
        exit(1);
    }
    TFile* outfile = new TFile(argv[2], "RECREATE");
    G4cout << argv[1] << "   " << argv[2] << G4endl;
    TDirectory *histotop = outfile->mkdir("histos");
    histotop->cd();

    TH1F *Gas = new TH1F("Gas", "Energy deposition", 100, 0., 10.);
    TH1F *innerAlWall = new TH1F("innerAlWall", "Energy deposition", 100, 0., 0.0001);
    TH1F *innerAuWall = new TH1F("innerAuWall", "Energy deposition", 100, 0., 0.0001);
    TH1F *innerMylarWall = new TH1F("innerMylarWall", "Energy deposition", 100, 0., 0.001);
    TH1F *innerWire = new TH1F("innerWire", "Energy deposition", 100, 0., 200.);
    TH1F *outerAlWall = new TH1F("outerAlWall", "Energy deposition", 100, 0., 0.0001);
    TH1F *outerWire = new TH1F("outerWire", "Energy deposition", 100, 0., 0.0001);
    //
    TFile fo(argv[1]);
    fo.GetListOfKeys()->Print();
    Event *event = new Event();
    TTree *T = (TTree*) fo.Get("T");
    T->SetBranchAddress("event.", &event);
    Int_t nevent = T->GetEntries();
    G4cout << " Nr. of Events:  " << nevent << G4endl;
    for (Int_t i = 0; i < nevent; i++) {
        T->GetEntry(i);
        std::map<G4String, std::vector<G4VHit*> >* hcmap = event->GetHCMap();
        std::map<G4String, std::vector<G4VHit*> >::iterator hciter;
        for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
            G4cout << " Collection name: " << (*hciter).first << G4endl;
            double sumE = 0.0; // total deposited energy
            std::vector<G4VHit*> hits = (*hciter).second;
            G4int NbHits = hits.size();
            std::vector<std::string> y = split((*hciter).first, '_');
            std::string Classname = y[1];
            G4cout << Classname << "   " << NbHits << G4endl;
            if (Classname == "Calorimeter") {
                for (G4int ii = 0; ii < NbHits; ii++) { // loop over calorimeter hits
                    CalorimeterHit* DRHit = dynamic_cast<CalorimeterHit*> (hits[ii]);
                    sumE = sumE + DRHit->GetEdep();
                    G4ThreeVector HitPosition = DRHit->GetPos();
                } // end loop over DRCalorimeterHits
                if ((*hciter).first == "GasVolume_Calorimeter_HC") Gas->Fill(sumE);
                if ((*hciter).first == "innerAlWallVolume_Calorimeter_HC")innerAlWall->Fill(sumE);
                if ((*hciter).first == "innerAuWallVolume_Calorimeter_HC")innerAuWall->Fill(sumE);
                if ((*hciter).first == "innerMylarWallVolume_Calorimeter_HC")innerMylarWall->Fill(sumE);
                if ((*hciter).first == "innerWireVolume_Calorimeter_HC")innerWire->Fill(sumE);
                if ((*hciter).first == "outerAlWallVolume_Calorimeter_HC")outerAlWall->Fill(sumE);
                if ((*hciter).first == "outerWireVolume_Calorimeter_HC")outerWire->Fill(sumE);
                G4cout << (*hciter).first << "  sumE:   " << sumE << G4endl;
            } // end (Classname == "Calorimeter") 
        } // end loop over events
        //

        //G4cout << " nr of bytes written:  " << outfile->Write() << G4endl;
    } // end loop over events 
    G4cout << " nr of bytes written:  " << outfile->Write() << G4endl;
}
