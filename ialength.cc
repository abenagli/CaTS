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
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
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
// standard C++ header files
//
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <map>

using namespace std;
string tags[4] = {"Material",
    "Physics List",
    "Particle",
    "Energy"};
string tagvalue[4];
map<string, vector<string> > FileNames;
TFile* outfile;
TH1F *hFirst;
TH1F* tempo;

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

void init(string fname) {
    ifstream infile;
    string line;
    infile.open(fname.c_str());
    if (infile.is_open()) {
        while (infile.good()) {
            getline(infile, line);
            vector<string> x = split(line, '_');
            if (x.size() == 5) {
                for (unsigned int i = 0; i < x.size() - 1; i++) {
                    tagvalue[i] = x[i];
                }
                string material = tagvalue[0];
                string physlist = tagvalue[1];
                string particle = tagvalue[2];
                if (FileNames.find(particle) == FileNames.end()) {
                    // not found
                    vector<string> tempvec;
                    tempvec.push_back(line);
                    FileNames.insert(std::pair<string, vector<string> >(particle, tempvec));
                } else {
                    // found
                    FileNames[particle].push_back(line);
                }
            }
        }
        infile.close();
    } else {
        cout << "couldn't Open infile! Exiting" << endl;
        exit(1);
    }
}

TFitResultPtr getlength(const char *filename) {
    TFile fo(filename);
    vector<string> x = split(filename, '_');
    for (unsigned int i = 0; i < x.size() - 1; i++) {
        tagvalue[i] = x[i];
    }
    vector<string> y = split(x[3], 'G');
    for (unsigned int i = 0; i < y.size() - 1; i++) {
        tagvalue[i + 3] = y[i];
    }
    double Ein = double(atof(tagvalue[3].c_str()));
    string particle = tagvalue[2];
    string energy = tagvalue[3];
    const char *p = particle.c_str();
    outfile->cd(p);
    char histoname[50];
    char histotitle[100];
    sprintf(histoname, "longi%fGeV", Ein);
    sprintf(histotitle, "longitudinal profile %fGeV", Ein);
    hFirst->Reset();
    tempo->Reset();
    //TH1F *hFirst = new TH1F("First", "Energy deposition/first cell", 100, 0., 30.);
    //TH1F* tempo = new TH1F("longi", "longitudinal profile", 2500, 0., 2500.);
    TH1F* start = new TH1F(histoname, histotitle, 2500, 0., 2500.);
    Event *event = new Event();
    TTree *T = (TTree*) fo.Get("T");
    T->SetBranchAddress("event.", &event);
    Int_t nevent = T->GetEntries();
    G4cout << " Nr. of Events:  " << nevent << G4endl;
    for (Int_t i = 0; i < nevent; i++) {
        tempo->Reset();
        T->GetEntry(i);
        std::map<G4String, std::vector<G4VHit*> >* hcmap = event->GetHCMap();
        std::map<G4String, std::vector<G4VHit*> >::iterator hciter;
        for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
            std::vector<G4VHit*> hits = (*hciter).second;
            std::vector<std::string> y = split((*hciter).first, '_');
            std::string Classname = y[1];
            if (Classname == "DRCalorimeter") {
                for (unsigned int ii = 0; ii < hits.size(); ii++) { // loop over calorimeter hits
                    DRCalorimeterHit* DRHit = dynamic_cast<DRCalorimeterHit*> (hits[ii]);
                    if (ii == 0) {
                        hFirst->Fill(DRHit->GetEdep());
                    }
                    G4ThreeVector HitPosition = DRHit->GetPos();
                    //Int_t num = 2500;          // 2500 crystal sheets
                    //Double_t cellsize = 1;     // 1 mm thickness of sheets
                    //Int_t layer = (HitPosition.z()/cellsize) + (num/2.);
                    //Double_t zp = HitPosition.z() + 1250.;
                    //
                    tempo->Fill(HitPosition.z() + 1250., DRHit->GetEdep());
                } // end loop over DRCalorimeterHits
            } // end (Classname == "DRCalorimeter")
            //
            // Now find the starting point of the shower
            // we assume that the shower is unlikely to start in the first layer
            //
            int nbins = tempo->GetNbinsX();
            double meanfirst = hFirst->GetMean();
            double sigmafirst = hFirst->GetRMS();
            double threshold = meanfirst + 2. * sigmafirst;
            int counter = 0;
            for (int jj = 0; jj < nbins; jj++) {
                if (tempo->GetBinContent(jj) > threshold) {
                    counter++;
                    if (counter > 6) {
                        start->Fill(jj - 6);
                        break;
                    }
                } else {
                    counter = 0;
                }
            } // end loop over bins
        }
    } // end loop over events
    //
    TF1 *ialength = new TF1("ialength", "[0]*TMath::Exp(-x/[1])", 0., 1000.);
    ialength->SetParameters(100., 200.0);
    ialength->SetParNames("norm", "lambda");
    TFitResultPtr fitpt = start->Fit("ialength", "S");
    G4cout << " nr of bytes written:  " << outfile->Write() << G4endl;
    return fitpt;
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
    hFirst = new TH1F("First", "Energy deposition/first cell", 100, 0., 30.);
    tempo = new TH1F("longi", "longitudinal profile", 2500, 0., 2500.);
    string fname(argv[1]);
    init(fname.c_str());
    outfile = new TFile(argv[2], "RECREATE");
    for (map<string, vector<string> >::iterator ii = FileNames.begin(); ii != FileNames.end(); ++ii) {
        cout << (*ii).first << endl;
        const char* dirname = (*ii).first.c_str();
        outfile->mkdir(dirname);
    }
    //
    vector<TFitResultPtr> fitptrvec;
    vector<double> Ein;
    for (map<string, vector<string> >::iterator ii = FileNames.begin(); ii != FileNames.end(); ++ii) {
        fitptrvec.clear();
        Ein.clear();
        vector<string>::const_iterator cii;
        for (cii = (*ii).second.begin(); cii != (*ii).second.end(); cii++) {
            vector<string> x = split(*cii, '_');
            for (unsigned int i = 0; i < x.size() - 1; i++) {
                tagvalue[i] = x[i];
            }
            vector<string> y = split(x[3], 'G');
            for (unsigned int i = 0; i < y.size() - 1; i++) {
                tagvalue[i + 3] = y[i];
            }
            Ein.push_back(double(atof(tagvalue[3].c_str())));
            const char *p;
            p = (*cii).c_str();
            fitptrvec.push_back(getlength(p));
            //	    cout << "Norm: "  << fitptr->Parameter(0) << " +/- " <<fitptr->ParError(0); 
            //cout << "  Lambda:  "  << fitptr->Parameter(1) << " +/- " <<fitptr->ParError(1) << endl; 
        }
        cout << "energies:  " << fitptrvec.size() << endl;
        Double_t energies_tvec[100];
        Double_t errenergies_tvec[100];
        Double_t lambda_tvec[100];
        Double_t errlambda_tvec[100];

        for (unsigned int i = 0; i < fitptrvec.size(); i++) {
            energies_tvec[i] = Ein[i];
            errenergies_tvec[i] = 0.001;
            lambda_tvec[i] = fitptrvec[i]->Parameter(1);
            errlambda_tvec[i] = fitptrvec[i]->ParError(1);
        }
        outfile->mkdir("top");
        outfile->cd("top");
        TCanvas *c = new TCanvas("c", "response", 200, 10, 1000, 800);
        TGraphErrors *gr_ref = new TGraphErrors(fitptrvec.size(), energies_tvec, lambda_tvec, errenergies_tvec, errlambda_tvec);
        gr_ref->SetLineColor(4);
        gr_ref->SetLineWidth(1);
        gr_ref->SetLineStyle(2);
        gr_ref->SetMarkerColor(4);
        gr_ref->SetMarkerStyle(20);
        gr_ref->SetMarkerSize(1.2);
        gr_ref->SetTitle("nuclear interacion length ");
        gr_ref->SetName("visible Energy");
        gr_ref->GetXaxis()->SetTitle("Ein [GeV]");
        gr_ref->GetYaxis()->SetTitle("IA [mm]");
        gr_ref->Draw("ACP");
        c->Write();

    }
}
