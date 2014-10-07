//Program that analyzes the profiles of hits produced with CaTS
//

// Include files
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include <TCanvas.h>
#include "TStyle.h"
#include "TFitResultPtr.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TVectorF.h"

#include "Cintex/Cintex.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "Event.hh"
#include "DRCalorimeterHit.hh"

using namespace std;

vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    return split(s, delim, elems);
}

int main(int argc, char** argv) {
    TSystem ts;
    gSystem->Load("libCintex");
    gSystem->Load("libClassesDict");
    //  ROOT::Cintex::Cintex::SetDebug(2);
    ROOT::Cintex::Cintex::Enable();
    if (argc < 2) cout << "Missing name of the file to read!" << endl;
    vector<string> pifnames;
    vector<double> pienergies;
    string tags[4] = {"Material",
        "Physics List",
        "Particle",
        "Energy"};
    char histoname[50];
    char histotitle[100];
    vector<TH1F*> vechEdepz;
    vector<TH2F*> vechEdepxy;
    vector<TH1F*> vechEdepx;
    vector<TH1F*> vechEdepy;
    vector<TH1F*> vechEdepR;
    vector<TH1F*> vechEdepR2;
    //
    vector<TH1F*> vechNCerenz;
    vector<TH2F*> vechNCerenxy;
    vector<TH1F*> vechNCerenx;
    vector<TH1F*> vechNCereny;
    vector<TH1F*> vechNCerenR;
    vector<TH1F*> vechNCerenR2;
    //
    vector<TFitResultPtr> vechEdepzfrprt;
    vector<TFitResultPtr> vechEdepxyfrptr;
    vector<TFitResultPtr> vechEdepxfrptr;
    vector<TFitResultPtr> vechEdepyfrptr;
    //
    vector<TFitResultPtr> vechNCerenzfrptr;
    vector<TFitResultPtr> vechNCerenxyfrptr;
    vector<TFitResultPtr> vechNCerenxfrptr;
    vector<TFitResultPtr> vechNCerenyfrptr;
    //
    
    //
    TFile* outfile;
    TDirectory *piontop;
    TDirectory *finaltop;


    //BEGIN INITIALIZATION

    string fname(argv[1]);
    string tagvalue[4];
    ifstream infile;
    string line;
    infile.open(fname.c_str());
    if (infile.is_open()) {
        while (infile.good()) {
            getline(infile, line);
            vector<string> x = split(line, '_');
            if (x.size() >= 4) {
                for (unsigned int i = 0; i < x.size() - 1; i++) {
                    tagvalue[i] = x[i];
                }
                vector<string> y = split(x[3], 'G');
                for (unsigned int i = 0; i < y.size() - 1; i++) {
                    tagvalue[i + 3] = y[i];
                }
                double evalue = double(atof(tagvalue[3].c_str()));
                string material = tagvalue[0];
                string physlist = tagvalue[1];
                string particle = tagvalue[2];
                string energy = tagvalue[3];
                if (particle == "pi-") {
                    pifnames.push_back(line);
                    pienergies.push_back(evalue);
                }
            }
        }
        infile.close();
        //
        // Now print out what we got:
        //
        cout << "Read data from file " << fname << endl;
        cout << "Nr. of Pion files: " << pifnames.size() << endl;
        cout << endl;
        cout << "pi energies (GeV): ";
        cout << endl;
        for (unsigned int i = 0; i < pienergies.size(); i++) {
            cout << ", " << pienergies[i];
        }
        cout << endl;
        string histofile = tagvalue[0] + "_" + tagvalue[1] + "_profilehistos.root";
        cout << "histo out file:  " << histofile << endl;
        outfile = new TFile(histofile.c_str(), "RECREATE");
        // create a subdirectory "top" in this file
        TDirectory *cdtop = outfile->mkdir("top");
        cdtop->cd(); // make the "top" directory the current directory
        finaltop = cdtop->mkdir("final");
        piontop = cdtop->mkdir("pions");
        piontop->cd();
    } else {
        cout << "couldn't Open infile! Exiting" << endl;
        exit(1);
    }
    
    int n_energies = pienergies.size();
    const int n_graphs = 5;
    string graphs[n_graphs] = {"z_mean","x_rms","y_rms","R_rms","R2_rms"};
    TVectorF * graph_vec[2*n_graphs];
    for(int i=0; i<2*n_graphs; i++)
    {
        graph_vec[i] = new TVectorF(n_energies);
    }
    
    TVectorF pienergies_tvec(n_energies);
    for(int i=0; i<n_energies; i++)
    {
       pienergies_tvec[i] = pienergies[i];
    }

    
    for (unsigned int index = 0; index < pienergies.size(); index++) {
        piontop->cd();
        double Ein = pienergies[index];
        sprintf(histoname, "Edepz%fGeV", Ein);
        sprintf(histotitle, "z-position of Energy deposition (Ein %f GeV)", Ein);
        vechEdepz.push_back(new TH1F(histoname, histotitle, 100, -1250., 1250.));
        //
        sprintf(histoname, "Edepxy%fGeV", Ein);
        sprintf(histotitle, "xy-position of Energy deposition (Ein %f GeV)", Ein);
        vechEdepxy.push_back(new TH2F(histoname, histotitle, 100, -1250., 1250., 100, -1250., 1250.));
        //
        sprintf(histoname, "Edepx%fGeV", Ein);
        sprintf(histotitle, "x-position of Energy deposition (Ein %f GeV)", Ein);
        vechEdepx.push_back(new TH1F(histoname, histotitle, 100, -1250., 1250.));
        //
        sprintf(histoname, "Edepy%fGeV", Ein);
        sprintf(histotitle, "y-position of Energy deposition (Ein %f GeV)", Ein);
        vechEdepy.push_back(new TH1F(histoname, histotitle, 100, -1250., 1250.));
        //
        sprintf(histoname, "EdepR%fGeV", Ein);
        sprintf(histotitle, "R-position of Energy deposition (Ein %f GeV)", Ein);
        vechEdepR.push_back(new TH1F(histoname, histotitle, 100, 0., 100.));
        //
        sprintf(histoname, "EdepR2%fGeV", Ein);
        sprintf(histotitle, "R2-position of Energy deposition (Ein %f GeV)", Ein);
        vechEdepR2.push_back(new TH1F(histoname, histotitle, 100, 0., 5000.));
        //
        sprintf(histoname, "NCerenz%fGeV", Ein);
        sprintf(histotitle, "z-position of Cerenkov photons (Ein %f GeV)", Ein);
        vechNCerenz.push_back(new TH1F(histoname, histotitle, 100, -1250., 1250.));
        //
        sprintf(histoname, "NCerenxy%fGeV", Ein);
        sprintf(histotitle, "xy-position of Cerenkov photons (Ein %f GeV)", Ein);
        vechNCerenxy.push_back(new TH2F(histoname, histotitle, 100, -1250., 1250., 100, -1250., 1250.));
        //
        sprintf(histoname, "NCerenx%fGeV", Ein);
        sprintf(histotitle, "x-position of Cerenkov photons (Ein %f GeV)", Ein);
        vechNCerenx.push_back(new TH1F(histoname, histotitle, 100, -1250., 1250.));
        //
        sprintf(histoname, "NCereny%fGeV", Ein);
        sprintf(histotitle, "y-position of Cerenkov photons (Ein %f GeV)", Ein);
        vechNCereny.push_back(new TH1F(histoname, histotitle, 100, -1250., 1250.));

        //
        sprintf(histoname, "NCerenR%fGeV", Ein);
        sprintf(histotitle, "R-position of Cerenkov photons (Ein %f GeV)", Ein);
        vechNCerenR.push_back(new TH1F(histoname, histotitle, 100, 0., 100.));
        //
        sprintf(histoname, "NCerenR2%fGeV", Ein);
        sprintf(histotitle, "R2-position of Cerenkov photons (Ein %f GeV)", Ein);
        vechNCerenR2.push_back(new TH1F(histoname, histotitle, 100, 0., 5000.));
        //
        TFile f(pifnames[index].c_str());
        //f.GetListOfKeys()->Print();
        Event *event = new Event();
        TTree *T = (TTree*) f.Get("T");
        T->SetBranchAddress("event.", &event);
        Int_t nevent = T->GetEntries();
        cout << endl << "Nr. of Events: " << nevent << endl;
        //loop over events
        for (Int_t i = 0; i < nevent; i++) {
            if( !(i % 20) ) cout << "Energy: " << Ein << " GeV - " << i*100/nevent << "%\r" << flush;
            T->GetEntry(i);
            map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
            //G4cout << "There are: " << hcmap->size() << "Hit Collections" << G4endl;
            map<G4String, vector<G4VHit*> >::iterator hciter;
            for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
                //G4cout << "Collection: " << (*hciter).first
                //        << "  Size:  " << (*hciter).second.size() << G4endl;
                vector<G4VHit*> hits = (*hciter).second;
                G4int NbHits = hits.size();
                vector<string> y = split((*hciter).first, '_');
                string Classname = y[1];
                if (Classname == "DRCalorimeter") {
                    for (G4int ii = 0; ii < NbHits; ii++) {
                        DRCalorimeterHit* DRHit = dynamic_cast<DRCalorimeterHit*> (hits[ii]);
                        G4ThreeVector pos = DRHit->GetPos();
                        Double_t Radiussq = 0.01 * pos.getX() * pos.getX() + 0.01 * pos.getY() * pos.getY();
                        Double_t Radius = TMath::Sqrt(Radiussq);
                        vechEdepz[index]->Fill(pos.getZ(), DRHit->GetEdep());
                        vechEdepxy[index]->Fill(pos.getX(), pos.getY(), DRHit->GetEdep());
                        vechEdepx[index]->Fill(pos.getX(), DRHit->GetEdep());
                        vechEdepy[index]->Fill(pos.getY(), DRHit->GetEdep());
                        vechEdepR[index]->Fill(Radius, DRHit->GetEdep());
                        vechEdepR2[index]->Fill(Radiussq, DRHit->GetEdep());
                        vechNCerenz[index]->Fill(pos.getZ(), DRHit->GetNCeren());
                        vechNCerenxy[index]->Fill(pos.getX(), pos.getY(), DRHit->GetNCeren());
                        vechNCerenx[index]->Fill(pos.getX(), DRHit->GetNCeren());
                        vechNCereny[index]->Fill(pos.getY(), DRHit->GetNCeren());
                        vechNCerenR[index]->Fill(Radius, DRHit->GetNCeren());
                        vechNCerenR2[index]->Fill(Radiussq, DRHit->GetNCeren());
                    }
                }
            }
        }
        cout << endl;
        // end loop over events
        Double_t scale = 100. / vechEdepz[index]->Integral();
        vechEdepz[index]->Scale(scale);
        vechEdepzfrprt.push_back(vechEdepz[index]->Fit("landau", "S"));
        Double_t scalen = 100. / vechNCerenz[index]->Integral();
        vechNCerenz[index]->Scale(scalen);
        vechNCerenzfrptr.push_back(vechNCerenz[index]->Fit("landau", "S"));
        piontop->cd();
        vechEdepz[index]->Write();
        vechEdepxy[index]->Write();
        vechEdepx[index]->Write();
        vechEdepy[index]->Write();
        vechEdepR[index]->Write();
        vechEdepR2[index]->Write();
        vechNCerenz[index]->Write();
        vechNCerenxy[index]->Write();
        vechNCerenx[index]->Write();
        vechNCereny[index]->Write();
        vechNCerenR[index]->Write();
        vechNCerenR2[index]->Write();

        (*(graph_vec[0]))[index]=vechEdepz[index]->GetMean();
        (*(graph_vec[1]))[index]=vechEdepx[index]->GetRMS();
        (*(graph_vec[2]))[index]=vechEdepy[index]->GetRMS();
        (*(graph_vec[3]))[index]=vechEdepR[index]->GetRMS();
        (*(graph_vec[4]))[index]=vechEdepR2[index]->GetRMS();

        (*(graph_vec[5]))[index]=vechNCerenz[index]->GetMean();
        (*(graph_vec[6]))[index]=vechNCerenx[index]->GetRMS();
        (*(graph_vec[7]))[index]=vechNCereny[index]->GetRMS();
        (*(graph_vec[8]))[index]=vechNCerenR[index]->GetRMS();
        (*(graph_vec[9]))[index]=vechNCerenR2[index]->GetRMS();

    }

    finaltop->cd();
    for(int i = 0; i<n_graphs;i++)
    {
        TMultiGraph * mg = new TMultiGraph( graphs[i].c_str(), graphs[i].c_str() );

        TGraph * gEdep   = new TGraph(pienergies_tvec,*(graph_vec[i]));
        gEdep->SetLineColor(1);

        TGraph * gNCeren = new TGraph(pienergies_tvec,*(graph_vec[i+5]));
        gNCeren->SetLineColor(2);

        mg->Add(gEdep);
        mg->Add(gNCeren);
        mg->Write();
    }
}



