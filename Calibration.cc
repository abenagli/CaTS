//---------------------------------------------------------------------------
//  
// routine to automate the energy correction for a dual read out calorimeter
// author: Hans Wenzel
// 
//---------------------------------------------------------------------------
#include "RunHeader.hh"
#include "Event.hh"
#include "DRTSCalorimeterHit.hh"
//
// Root header files:
//
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h" 
#include "TMath.h" 
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "Cintex/Cintex.h"
#include "TLegend.h"
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
TFile* outfile;
TH2F * banana2;
TProfile *hprof2;
vector<TH2F*> eban;
vector<TH2F*> ban;
//
// number of slices in the banana plot
//
const unsigned int nslices = 20;
TH1F * hslice[nslices];
double slicewidth;
//
// calibration factors for scintilation (S) and cerenkov (C)
//
double scorr, ccorr;

double corrf[4];

vector<string> efnames;
vector<double> eenergies;
vector<string> pifnames;
vector<double> pienergies;

TDirectory *electrontop; // directory to collect electron related histograms
TDirectory *piontop; // directory to collect pion related histograms
TDirectory *finaltop; // directory to collect summary plots and graphs
TDirectory *femtop; // directory to collect em fraction related plots

//
// raw distributions without any correction 
// and pointers to results of gaussian fit to this distributions
//
vector<TH1F*> electronEdep;
vector<TH1F*> electronCeren;
vector<TFitResultPtr> electronEdepptr;
vector<TFitResultPtr> electronCerenptr;
//
// electron and pion histograms after simple energy calibration 
// and pointers to results of gaussian fit to this distributions
//
vector<TH1F*> celectronEdep;
vector<TH1F*> celectronCeren;
vector<TH1F*> cpionEdep;
vector<TH1F*> cpionCeren;

vector<TFitResultPtr> cpionEdepptr;
vector<TFitResultPtr> cpionCerenptr;
vector<TFitResultPtr> celectronEdepptr;
vector<TFitResultPtr> celectronCerenptr;

string tags[4] = {"Material",
    "Physics List",
    "Particle",
    "Energy"};
string tagvalue[4];

//
// dual read out corrected distributions and pointers 
// to results of gaussian fit to this distributions
//
vector<TH1F*> drcpionEdep;
vector<TFitResultPtr> drcpionEdepptr;
//
// hadronic response distributions:
// 
vector<TH1F*> Pi0histo;
vector<TH1F*> Pi0hshisto;
vector<TH1F*> Pi0hchisto;
vector<TH1F*> FEMhisto;
vector<TH1F*> FEMhshisto;
vector<TH1F*> FEMhchisto;
vector<TH1F*> FEMmeanhisto;

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

Bool_t OddEv(Int_t m) {
    return m % 2;
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
                /*
                                vector<string> y = split(x[3], 'G');
                                for (unsigned int i = 0; i < y.size() - 1; i++) {
                                    tagvalue[i + 3] = y[i];
                                }
                 */
                double evalue = double(atof(tagvalue[3].c_str()));
                string material = tagvalue[0];
                string physlist = tagvalue[1];
                string particle = tagvalue[2];
                string energy = tagvalue[3];
                if (particle == "e-") {
                    efnames.push_back(line);
                    eenergies.push_back(evalue);
                }
                if (particle == "pi-") {
                    pifnames.push_back(line);
                    pienergies.push_back(evalue + 0.13957);
                }
            }
        }
        infile.close();
        //
        // Now print out what we got:
        //
        cout << "Read data from file " << fname << endl;
        cout << "Nr. of electron files: " << efnames.size() << endl;
        cout << "Nr. of Pion files: " << pifnames.size() << endl;
        cout << "e  energies (GeV): ";
        for (unsigned int i = 0; i < eenergies.size(); i++) {
            cout << eenergies[i] << ", ";
        }
        cout << endl;
        cout << "pi energies (GeV): ";
        for (unsigned int i = 0; i < pienergies.size(); i++) {
            cout << pienergies[i] << ", ";
        }
        cout << endl;
    } else {
        cout << "couldn't Open infile! Exiting" << endl;
        exit(1);
    }
    string histofile = tagvalue[0] + "_" + tagvalue[1] + "_histos_test.root";
    cout << "histo out file:  " << histofile << endl;
    outfile = new TFile(histofile.c_str(), "RECREATE");
    // create a subdirectory "top" in this file
    TDirectory *cdtop = outfile->mkdir("top");
    cdtop->cd(); // make the "top" directory the current directory
    electrontop = cdtop->mkdir("electrons");
    piontop = cdtop->mkdir("pions");
    femtop = cdtop->mkdir("fem");
    finaltop = cdtop->mkdir("final");
    finaltop->cd();

    banana2 = new TH2F("banana2", "S/Ein vs. C/S", 100, 0, 2., 100, 0, 2.);
    banana2->GetXaxis()->SetTitle("C/S");
    banana2->GetYaxis()->SetTitle("S/Ein");

    hprof2 = new TProfile("hprof2", "S/Ein vs. C/S", nslices, 0, 1.1, 0, 1.1);
    hprof2->GetXaxis()->SetTitle("C/S");
    hprof2->GetYaxis()->SetTitle("S/Ein");

    slicewidth = 1.0 / double(nslices);
    for (unsigned int ii = 0; ii < nslices; ii++) {
        char histoname[20];
        char histotitle[50];
        sprintf(histotitle, "C/S slice %d", ii);
        sprintf(histoname, "slice%d", ii);
        hslice[ii] = new TH1F(histoname, histotitle, 100, 0, 3.0);
        hslice[ii]->GetXaxis()->SetTitle("S/Ein");
        hslice[ii]->GetYaxis()->SetTitle("Evts.");
        hslice[ii]->SetFillColor(7);
        hslice[ii]->SetLineColor(4);
    }
    cdtop->cd();
    cout << " initialized" << endl;
}
//
// Dual read out correction function:
//
// Ein = S/fcorr(C/S)
// 

double fcorr(double ratio) {

    double x = ratio;
    double fc = 0;
    if (ratio < 0.4) {
        x = 0.4;
    }
    if (x > 1) {
        fc = 1.;
        return fc;
    }
    double p0 = corrf[0];
    double p1 = corrf[1];
    double p2 = corrf[2];
    double p3 = corrf[3];
    fc = p0 + p1 * x + p2 * x * x + p3 * x * x*x;

    return fc;
}

void pidist(string particle, unsigned int index) {
    char histoname[50];
    char histotitle[100];

    if (particle == "Electrons") {
        electrontop->cd();
        double Ein = eenergies[index];

        TFile f(efnames[index].c_str());
        /*
        Event *event = new Event();
        TTree *T = (TTree*) f.Get("T");
        T->SetBranchAddress("event.", &event);
        Int_t nevent = T->GetEntries();
        G4cout << "Nr. of Events: " << nevent << G4endl;
         */

        RunHeader *runh = new RunHeader();
        Event *event = new Event();
        TTree *T = (TTree*) f.Get("Events");
        T->SetBranchAddress("event.", &event);
        TBranch* fevtbranch = T->GetBranch("event.");
        TTree *Trh = (TTree*) f.Get("Runheader");
        Trh->SetBranchAddress("RunHeader.", &runh);
        TBranch* frunhbranch = Trh->GetBranch("RunHeader.");
        frunhbranch->GetEntry(0);
        runh->Print();
        Int_t nevent = T->GetEntries();
        //G4cout << "Nr. of Events: " << nevent << G4endl;






        Double_t MeanNCeren = 0.0;
        Double_t MeanEdep = 0.0;
        Double_t SumEp = 0.0;
        Double_t SumNCerenp = 0.0;

        //First loop over events
        for (Int_t i = 0; i < nevent; i++) {
            Double_t SumE = 0.0;
            Int_t SumNCeren = 0;

            T->GetEntry(i);
            map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
            map<G4String, vector<G4VHit*> >::iterator hciter;
            for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
                vector<G4VHit*> hits = (*hciter).second;
                G4int NbHits = hits.size();
                vector<string> y = split((*hciter).first, '_');
                string Classname = y[1];
                if (Classname == "DRTSCalorimeter") {
                    for (G4int ii = 0; ii < NbHits; ii++) {
                        DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
                        G4ThreeVector pos = DRHit->GetPos();
                        double z = pos.getZ();
                        int num = 100;
                        int cellsize = 25;
                        int layer = z / cellsize + num / 2;
                        //Bool_t oetl = OddEv(layer);
                        //if (oetl == kTRUE)
                        //Function to read 1 out of 3 layers
                        //double testint = (layer - 1.)/3.;
                        //double testint = (layer - 2.)/3.;
                        //double testint = (layer - 3.)/3.;
                        //Function to read 1 out of 4 layers
                        //double testint = (layer - 1.)/4.;
                        //double testint2 = fmod(testint,1);
                        //if(testint2==0)
                        //{
                        SumE = SumE + (DRHit->GetEdep() / 1000);
                        SumNCeren = SumNCeren + DRHit->GetNCeren();
                        //G4cout<<"Layer:  " <<layer<<G4endl;
                        //}
                    } // End loop over Hits
                }
            }// End loop over Hit collections

            double corrE = SumE*scorr;
            double corrC = SumNCeren*ccorr;

            SumNCerenp = SumNCerenp + corrC;
            SumEp = SumEp + corrE;

        }// End loop over Events

        MeanNCeren = SumNCerenp / double(nevent);
        MeanEdep = SumEp / double(nevent);

        G4cout << "Electron Corrected Mean Edep:  " << MeanEdep << G4endl;
        G4cout << "Electron Corrected Mean NCeren:  " << MeanNCeren << G4endl;

        electrontop->cd();

        Double_t SumNCeren2 = 0.0;
        Double_t SumE2 = 0.0;

        //Second loop over events(Variation Calculation)
        for (Int_t i = 0; i < nevent; i++) {
            Double_t SumE = 0.0;
            Int_t SumNCeren = 0;

            T->GetEntry(i);
            map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
            map<G4String, vector<G4VHit*> >::iterator hciter;
            for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
                vector<G4VHit*> hits = (*hciter).second;
                G4int NbHits = hits.size();
                vector<string> y = split((*hciter).first, '_');
                string Classname = y[1];
                if (Classname == "DRTSCalorimeter") {
                    for (G4int ii = 0; ii < NbHits; ii++) {
                        DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
                        G4ThreeVector pos = DRHit->GetPos();
                        double z = pos.getZ();
                        int num = 100;
                        int cellsize = 25;
                        int layer = z / cellsize + num / 2;
                        //Bool_t oetl = OddEv(layer);
                        //double testint = (layer - 1.)/ 3.;
                        //double testint = (layer - 2.)/3;
                        //double testint = (layer - 3.)/3.;
                        //double testint = (layer - 1.)/ 4.;
                        //double testint2 = fmod(testint,1);
                        //if(testint2==0)
                        //if (oetl == kTRUE)
                        //{
                        //cout<<"Layer:  "<<layer<<endl;
                        SumE = SumE + (DRHit->GetEdep() / 1000);
                        SumNCeren = SumNCeren + DRHit->GetNCeren();
                        //}
                    } // End loop over Hits
                }
            }// End loop over Hit collections
            double corrE = SumE*scorr;
            double corrC = SumNCeren*ccorr;

            SumE2 = SumE2 + (corrE - MeanEdep)*(corrE - MeanEdep);
            SumNCeren2 = SumNCeren2 + (corrC - MeanNCeren)*(corrC - MeanNCeren);
        }// End loop over Events

        Double_t VarEdep = SumE2 / (double(nevent) - 1.);
        Double_t VarNCeren = SumNCeren2 / (double(nevent) - 1.);

        Double_t sigmaEdep = sqrt(VarEdep);
        Double_t sigmaNCeren = sqrt(VarNCeren);

        G4cout << "Electron Corrected sigma Edep:  " << sigmaEdep << G4endl;
        G4cout << "Electron Corrected sigma NCeren:  " << sigmaNCeren << G4endl;

        //Third loop over events (Fill histograms)
        sprintf(histoname, "celectronEdep%fGeV", Ein);
        sprintf(histotitle, "electron total Energy deposition (corrected) (Ein %f GeV)", Ein);
        celectronEdep.push_back(new TH1F(histoname, histotitle, 250, MeanEdep - 5 * sigmaEdep, MeanEdep + 5 * sigmaEdep));
        //celectronEdep.push_back(new TH1F(histoname, histotitle, 250, 0, 1.1*Ein));
        celectronEdep[index]->GetXaxis()->SetTitle("E[GeV]");
        celectronEdep[index]->GetYaxis()->SetTitle("Nr of Evts.");
        celectronEdep[index]->SetFillColor(7);
        celectronEdep[index]->SetLineColor(4);

        sprintf(histoname, "celectronCeren%fGeV", Ein);
        sprintf(histotitle, "electron total cerenkov Energy deposition (corrected) (Ein %f GeV)", Ein);
        celectronCeren.push_back(new TH1F(histoname, histotitle, 500, MeanNCeren - 5 * sigmaNCeren, MeanNCeren + 5 * sigmaNCeren));
        celectronCeren[index]->GetXaxis()->SetTitle("E[GeV]");
        celectronCeren[index]->GetYaxis()->SetTitle("Nr of Evts.");
        celectronCeren[index]->SetFillColor(7);
        celectronCeren[index]->SetLineColor(4);

        sprintf(histoname, "banplot%fGeV", Ein);
        sprintf(histotitle, "S/E, C/S histogram (Ein %f GeV)", Ein);
        eban.push_back(new TH2F(histoname, histotitle, 500, 0.0, 2.0, 500, 0.0, 2.0));
        eban[index]->GetXaxis()->SetTitle("C/S");
        eban[index]->GetYaxis()->SetTitle("S/E_{in}");
        eban[index]->SetFillColor(7);
        eban[index]->SetLineColor(4);

        for (Int_t i = 0; i < nevent; i++) {
            Double_t SumE = 0.0;
            Int_t SumNCeren = 0;
            T->GetEntry(i);
            map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
            map<G4String, vector<G4VHit*> >::iterator hciter;
            for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
                vector<G4VHit*> hits = (*hciter).second;
                G4int NbHits = hits.size();
                vector<string> y = split((*hciter).first, '_');
                string Classname = y[1];
                if (Classname == "DRTSCalorimeter") {
                    for (G4int ii = 0; ii < NbHits; ii++) {
                        DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
                        G4ThreeVector pos = DRHit->GetPos();
                        double z = pos.getZ();
                        int num = 100;
                        int cellsize = 25;
                        int layer = z / cellsize + num / 2;
                        //Bool_t oetl = OddEv(layer);
                        //double testint = (layer - 1.)/ 3.;
                        //double testint = (layer -2.)/3;
                        //double testint = (layer - 3.)/3.;
                        //double testint = (layer - 1.)/ 4.;
                        //double testint2 = fmod(testint,1);
                        //if(testint2==0)
                        //if (oetl == kTRUE)
                        //{
                        SumE = SumE + (DRHit->GetEdep());
                        SumNCeren = SumNCeren + DRHit->GetNCeren();
                        //}
                    } // End loop over Hits
                }
            }// End loop over Hit collections

            double corrE = SumE * scorr / 1000.;
            double corrC = SumNCeren*ccorr;
            celectronEdep[index]->Fill(corrE);
            celectronCeren[index]->Fill(corrC);
            banana2->Fill(corrC / corrE, corrE / Ein);
            hprof2->Fill(corrC / corrE, corrE / Ein);
            eban[index]->Fill(corrC / corrE, corrE / Ein);

            int slice = (corrC / corrE) / slicewidth;
            if (slice > 19) {
                slice = 19;
            }
            hslice[slice]->Fill(corrE / Ein);
        }// End loop over Events

        electrontop->cd();
        celectronEdepptr.push_back(celectronEdep[index]->Fit("gaus", "S"));
        celectronEdep[index]->Write();
        celectronCerenptr.push_back(celectronCeren[index]->Fit("gaus", "S"));
        celectronCeren[index]->Write();
        eban[index]->Write();

    } else if (particle == "Pions") {
        double Ein = pienergies[index];
        femtop->cd();
        sprintf(histoname, "Pi0histo%fGeV", Ein);
        sprintf(histotitle, "energy fraction deposited by pi0 (Ein %f GeV)", Ein);
        Pi0histo.push_back(new TH1F(histoname, histotitle, 50, 0.0, 1.));
        Pi0histo[index]->GetXaxis()->SetTitle("frac. E");
        Pi0histo[index]->GetYaxis()->SetTitle("Nr of Evts.");
        Pi0histo[index]->SetFillColor(7);
        Pi0histo[index]->SetLineColor(4);
        //
        sprintf(histoname, "Pi0hshisto%fGeV", Ein);
        sprintf(histotitle, "hs (Pi0) (Ein %f GeV)", Ein);
        Pi0hshisto.push_back(new TH1F(histoname, histotitle, 50, 0.0, 1.));
        Pi0hshisto[index]->GetXaxis()->SetTitle("hs");
        Pi0hshisto[index]->GetYaxis()->SetTitle("Nr of Evts.");
        Pi0hshisto[index]->SetFillColor(7);
        Pi0hshisto[index]->SetLineColor(4);
        //
        sprintf(histoname, "Pi0hchisto%fGeV", Ein);
        sprintf(histotitle, "hc (Pi0) (Ein %f GeV)", Ein);
        Pi0hchisto.push_back(new TH1F(histoname, histotitle, 50, 0.0, 1.));
        Pi0hchisto[index]->GetXaxis()->SetTitle("hc");
        Pi0hchisto[index]->GetYaxis()->SetTitle("Nr of Evts.");
        Pi0hchisto[index]->SetFillColor(7);
        Pi0hchisto[index]->SetLineColor(4);
        //
        sprintf(histoname, "FEMhisto%fGeV", Ein);
        sprintf(histotitle, "EM energy fraction (Ein %f GeV)", Ein);
        FEMhisto.push_back(new TH1F(histoname, histotitle, 50, 0.0, 1.));
        FEMhisto[index]->GetXaxis()->SetTitle("frac. E");
        FEMhisto[index]->GetYaxis()->SetTitle("Nr of Evts.");
        FEMhisto[index]->SetFillColor(7);
        FEMhisto[index]->SetLineColor(4);
        //
        sprintf(histoname, "FEMhshisto%fGeV", Ein);
        sprintf(histotitle, "hs (FEM) (Ein %f GeV)", Ein);
        FEMhshisto.push_back(new TH1F(histoname, histotitle, 50, 0.0, 1.));
        FEMhshisto[index]->GetXaxis()->SetTitle("hs");
        FEMhshisto[index]->GetYaxis()->SetTitle("Nr of Evts.");
        FEMhshisto[index]->SetFillColor(7);
        FEMhshisto[index]->SetLineColor(4);
        //
        sprintf(histoname, "FEMhchisto%fGeV", Ein);
        sprintf(histotitle, "hc (FEM) (Ein %f GeV)", Ein);
        FEMhchisto.push_back(new TH1F(histoname, histotitle, 50, 0.0, 1.));
        FEMhchisto[index]->GetXaxis()->SetTitle("hc");
        FEMhchisto[index]->GetYaxis()->SetTitle("Nr of Evts.");
        FEMhchisto[index]->SetFillColor(7);
        FEMhchisto[index]->SetLineColor(4);
        //
        piontop->cd();

        TFile f(pifnames[index].c_str());
        
        RunHeader *runh = new RunHeader();
        Event *event = new Event();
        TTree *T = (TTree*) f.Get("Events");
        T->SetBranchAddress("event.", &event);
        TBranch* fevtbranch = T->GetBranch("event.");
        TTree *Trh = (TTree*) f.Get("Runheader");
        Trh->SetBranchAddress("RunHeader.", &runh);
        TBranch* frunhbranch = Trh->GetBranch("RunHeader.");
        frunhbranch->GetEntry(0);
        //runh->Print();
        //Int_t nevent = T->GetEntries();
        //G4cout << "Nr. of Events: " << nevent << G4endl;


        
        
        
        /*
        Event *event = new Event();
        TTree *T = (TTree*) f.Get("T");
        T->SetBranchAddress("event.", &event);
*/
        Int_t nevent = T->GetEntries();
        G4cout << "Nr. of Events: " << nevent << G4endl;

        Double_t MeanNCeren = 0.0;
        Double_t MeanEdep = 0.0;
        Double_t SumEp = 0.0;
        Int_t SumNCerenp = 0.0;

        //First loop over events (Mean)
        for (Int_t i = 0; i < nevent; i++) {
            Double_t SumE = 0.0;
            Int_t SumNCeren = 0;

            T->GetEntry(i);
            map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
            map<G4String, vector<G4VHit*> >::iterator hciter;
            for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
                vector<G4VHit*> hits = (*hciter).second;
                G4int NbHits = hits.size();
                vector<string> y = split((*hciter).first, '_');
                string Classname = y[1];
                if (Classname == "DRTSCalorimeter") {
                    for (G4int ii = 0; ii < NbHits; ii++) {
                        DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
                        G4ThreeVector pos = DRHit->GetPos();
                        double z = pos.getZ();
                        int num = 100;
                        int cellsize = 25;
                        int layer = z / cellsize + num / 2;
                        //Bool_t oetl = OddEv(layer);
                        //double testint = (layer - 1.)/3.; 
                        //double testint = (layer - 2.)/3.;
                        //double testint = (layer - 3.)/3.;
                        //double testint = (layer - 1.)/ 4.;
                        //double testint2 = fmod(testint,1);
                        //if(testint2==0)
                        //if (oetl == kTRUE)
                        //{
                        SumE = SumE + (DRHit->GetEdep() / 1000.);
                        SumNCeren = SumNCeren + DRHit->GetNCeren();
                        //}
                    } // End loop over Hits
                }
            }// End loop over Hit collections
            double corrE = SumE*scorr;
            double corrC = SumNCeren*ccorr;

            SumNCerenp = SumNCerenp + corrC;
            SumEp = SumEp + corrE;

        } // end loop over evts.

        MeanNCeren = SumNCerenp / double(nevent);
        MeanEdep = SumEp / double(nevent);

        G4cout << "Pion MeanEdep:  " << MeanEdep << G4endl;
        G4cout << "Pion MeanNCeren:  " << MeanNCeren << G4endl;

        Double_t SumNCeren2 = 0.0;
        Double_t SumE2 = 0.0;

        //Second loop over events (Variation calculation)
        for (Int_t i = 0; i < nevent; i++) {
            Double_t SumE = 0.0;
            Int_t SumNCeren = 0;

            T->GetEntry(i);
            map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
            map<G4String, vector<G4VHit*> >::iterator hciter;
            for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
                vector<G4VHit*> hits = (*hciter).second;
                G4int NbHits = hits.size();
                vector<string> y = split((*hciter).first, '_');
                string Classname = y[1];
                if (Classname == "DRTSCalorimeter") {
                    for (G4int ii = 0; ii < NbHits; ii++) {
                        DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
                        G4ThreeVector pos = DRHit->GetPos();
                        double z = pos.getZ();
                        int num = 100;
                        int cellsize = 25;
                        int layer = z / cellsize + num / 2;
                        //Bool_t oetl = OddEv(layer);
                        //double testint = (layer - 1.)/3.; 
                        //double testint = (layer - 2.)/3.;
                        //double testint = (layer - 3.)/3.;
                        //double testint = (layer - 1.)/4.;
                        //double testint2 = fmod(testint,1);
                        //if(testint2==0)
                        //if (oetl == kTRUE)
                        //{
                        SumE = SumE + (DRHit->GetEdep() / 1000.);
                        SumNCeren = SumNCeren + DRHit->GetNCeren();
                        //}
                    } // End loop over Hits
                }
            }// End loop over Hit collections

            double corrE = SumE*scorr;
            double corrC = SumNCeren*ccorr;

            SumE2 = SumE2 + (corrE - MeanEdep)*(corrE - MeanEdep);
            SumNCeren2 = SumNCeren2 + (corrC - MeanNCeren)*(corrC - MeanNCeren);
        } // end loop over evts.

        Double_t VarEdep = SumE2 / (double(nevent) - 1.);
        Double_t VarNCeren = SumNCeren2 / (double(nevent) - 1.);

        Double_t sigmaEdep = sqrt(VarEdep);
        Double_t sigmaNCeren = sqrt(VarNCeren);

        G4cout << "Pion sigmaEdep:  " << sigmaEdep << G4endl;
        G4cout << "PionNCeren:  " << sigmaNCeren << G4endl;

        piontop->cd();

        sprintf(histoname, "cpionEdep%fGeV", Ein);
        sprintf(histotitle, "pion total Energy deposition (corrected) (Ein %f GeV)", Ein);
        cpionEdep.push_back(new TH1F(histoname, histotitle, 250, MeanEdep - 5 * sigmaEdep, MeanEdep + 5 * sigmaEdep));
        //cpionEdep.push_back(new TH1F(histoname, histotitle, 250, 0, 1.1*Ein));
        cpionEdep[index]->GetXaxis()->SetTitle("E[GeV]");
        cpionEdep[index]->GetYaxis()->SetTitle("Nr of Evts.");
        cpionEdep[index]->SetFillColor(7);
        cpionEdep[index]->SetLineColor(4);

        sprintf(histoname, "cpionCeren%fGeV", Ein);
        sprintf(histotitle, "pion total cerenkov Energy deposition (corrected) (Ein %f GeV)", Ein);
        cpionCeren.push_back(new TH1F(histoname, histotitle, 500, MeanNCeren - 5 * sigmaNCeren, MeanNCeren + 5 * sigmaNCeren));
        cpionCeren[index]->GetXaxis()->SetTitle("E[GeV]");
        cpionCeren[index]->GetYaxis()->SetTitle("Nr of Evts.");
        cpionCeren[index]->SetFillColor(7);
        cpionCeren[index]->SetLineColor(4);

        sprintf(histoname, "banplot%fGeV", Ein);
        sprintf(histotitle, "S/E, C/S histogram (Ein %f GeV)", Ein);
        ban.push_back(new TH2F(histoname, histotitle, 500, 0.0, 2.0, 500, 0.0, 2.0));
        ban[index]->GetXaxis()->SetTitle("C/S");
        ban[index]->GetYaxis()->SetTitle("S/E_{in}");
        ban[index]->SetFillColor(7);
        ban[index]->SetLineColor(4);

        //Third loop over events (Fill histograms)
        for (Int_t i = 0; i < nevent; i++) {
            Double_t SumE = 0.0;
            Double_t SumEEM = 0.0;
            Int_t SumNCeren = 0;

            T->GetEntry(i);
            map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
            map<G4String, vector<G4VHit*> >::iterator hciter;
            for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
                vector<G4VHit*> hits = (*hciter).second;
                G4int NbHits = hits.size();
                vector<string> y = split((*hciter).first, '_');
                string Classname = y[1];
                if (Classname == "DRTSCalorimeter") {
                    for (G4int ii = 0; ii < NbHits; ii++) {
                        DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
                        G4ThreeVector pos = DRHit->GetPos();
                        double z = pos.getZ();
                        int num = 100;
                        int cellsize = 25;
                        int layer = z / cellsize + num / 2;
                        //Bool_t oetl = OddEv(layer);
                        //double testint = (layer - 1.)/3.; 
                        //double testint = (layer - 2.)/3.;
                        //double testint = (layer - 3.)/3.;
                        //double testint = (layer - 1.)/4.;
                        //double testint2 = fmod(testint,1);
                        //if(testint2==0)
                        //if (oetl == kTRUE)
                        //{
                        SumE = SumE + DRHit->GetEdep();
                        SumEEM = SumEEM + DRHit->GetEdepEM();
                        SumNCeren = SumNCeren + DRHit->GetNCeren();
                        //}
                    } // End loop over Hits
                }
            }// End loop over Hit collections

            double corrE = 0.001 * SumE * scorr;
            double corrC = SumNCeren*ccorr;
            double fem = (0.001 * SumEEM) / Ein;
            FEMhisto[index]->Fill(fem);
            /*
            double pi0fem = event->GetPi0Energy() / Ein;
            //G4cout << "pi0fem"<< pi0fem<<G4endl;
            Pi0histo[index]->Fill(pi0fem);
            double hsfem = ((corrE / Ein) - fem) / (1. - fem);
            FEMhshisto[index]->Fill(hsfem);
            double hspi0 = ((corrE / Ein) - pi0fem) / (1. - pi0fem);
            Pi0hshisto[index]->Fill(hspi0);
            double hcfem = ((corrC / Ein) - fem) / (1. - fem);
            FEMhchisto[index]->Fill(hcfem);
            double hcpi0 = ((corrC / Ein) - pi0fem) / (1. - pi0fem);
            Pi0hchisto[index]->Fill(hcpi0);
             */
            cpionEdep[index]->Fill(corrE);
            cpionCeren[index]->Fill(corrC);

            banana2->Fill(corrC / corrE, corrE / Ein);
            hprof2->Fill(corrC / corrE, corrE / Ein);
            piontop->cd();
            ban[index]->Fill(corrC / corrE, corrE / Ein);

            int slice = (corrC / corrE) / slicewidth;
            if (slice > 19) {
                slice = 19;
            }
            hslice[slice]->Fill(corrE / Ein);
        } // end loop over evts.

        piontop->cd();
        cpionEdepptr.push_back(cpionEdep[index]->Fit("gaus", "S"));
        cpionEdep[index]->Write();
        cpionCerenptr.push_back(cpionCeren[index]->Fit("gaus", "S"));
        cpionCeren[index]->Write();
        ban[index]->Write();

        femtop->cd();
        FEMhisto[index]->Write();
        G4cout << "Mean:  " << FEMhisto[index]->GetMean() << G4endl;
        Pi0histo[index]->Write();
        FEMhshisto[index]->Write();
        Pi0hshisto[index]->Write();
        FEMhchisto[index]->Write();
        Pi0hchisto[index]->Write();

        //FEMmeanhisto[index]->Fill(FEMhisto[index]->GetMean());
    }

}

void CorrectedE(unsigned int index) {
    char histoname[50];
    char histotitle[100];
    piontop->cd();
    double Ein = pienergies[index];

    TFile f(pifnames[index].c_str());
    /*
    Event *event = new Event();
    TTree *T = (TTree*) f.Get("T");
    T->SetBranchAddress("event.", &event);
    Int_t nevent = T->GetEntries();
*/
    
    
        RunHeader *runh = new RunHeader();
        Event *event = new Event();
        TTree *T = (TTree*) f.Get("Events");
        T->SetBranchAddress("event.", &event);
        TBranch* fevtbranch = T->GetBranch("event.");
        TTree *Trh = (TTree*) f.Get("Runheader");
        Trh->SetBranchAddress("RunHeader.", &runh);
        TBranch* frunhbranch = Trh->GetBranch("RunHeader.");
        frunhbranch->GetEntry(0);
      //  runh->Print();
        Int_t nevent = T->GetEntries();
    //    G4cout << "Nr. of Events: " << nevent << G4endl;


    
    
    
    
    
    
    
    Double_t MeanEdep = 0.0;
    Double_t MeanNCeren = 0.0;
    Double_t SumEp = 0.0;
    Int_t SumNCerenp = 0.0;

    //First loop over events(Get Mean)
    for (Int_t i = 0; i < nevent; i++) {
        Double_t SumE = 0.0;
        Int_t SumNCeren = 0;

        T->GetEntry(i);
        map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
        map<G4String, vector<G4VHit*> >::iterator hciter;
        for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
            vector<G4VHit*> hits = (*hciter).second;
            G4int NbHits = hits.size();
            vector<string> y = split((*hciter).first, '_');
            string Classname = y[1];
            if (Classname == "DRTSCalorimeter") {
                for (G4int ii = 0; ii < NbHits; ii++) {
                    DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
                    G4ThreeVector pos = DRHit->GetPos();
                    double z = pos.getZ();
                    int num = 100;
                    int cellsize = 25;
                    int layer = z / cellsize + num / 2;
                    //Bool_t oetl = OddEv(layer);
                    //double testint = (layer - 1.)/3.; 
                    //double testint = (layer - 2.)/3.;
                    //double testint = (layer - 3.)/3.;
                    //double testint = (layer - 1.)/4.;
                    //double testint2 = fmod(testint,1);
                    //if(testint2==0)
                    //if (oetl == kTRUE)
                    //{
                    SumE = SumE + DRHit->GetEdep();
                    SumNCeren = SumNCeren + DRHit->GetNCeren();
                    //}
                } // End loop over Hits
            }
        }// End loop over Hit collections
        double corrE = SumE * scorr / 1000.;
        double corrC = SumNCeren*ccorr;
        double ratio = corrC / corrE;
        double drcEcorr = corrE / fcorr(ratio);

        SumNCerenp = SumNCerenp + corrC;
        SumEp = SumEp + drcEcorr;

    } // end loop over evts.

    MeanNCeren = SumNCerenp / double(nevent);
    MeanEdep = SumEp / double(nevent);

    //Second loop over events (Variance Calculation)
    Int_t SumNCeren2 = 0;
    Double_t SumE2 = 0.0;
    Double_t VarEdep = 0.0;
    Double_t VarNCeren = 0.0;

    for (Int_t i = 0; i < nevent; i++) {
        Double_t SumE = 0.0;
        Int_t SumNCeren = 0;
        T->GetEntry(i);
        map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
        map<G4String, vector<G4VHit*> >::iterator hciter;
        for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
            vector<G4VHit*> hits = (*hciter).second;
            G4int NbHits = hits.size();
            vector<string> y = split((*hciter).first, '_');
            string Classname = y[1];
            if (Classname == "DRTSCalorimeter") {
                for (G4int ii = 0; ii < NbHits; ii++) {
                    DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
                    G4ThreeVector pos = DRHit->GetPos();
                    double z = pos.getZ();
                    int num = 100;
                    int cellsize = 25;
                    int layer = z / cellsize + num / 2;
                    //Bool_t oetl = OddEv(layer);
                    //double testint = (layer - 1.)/3.; 
                    //double testint = (layer - 2.)/3.;
                    //double testint = (layer - 3.)/3.;
                    //double testint = (layer - 1.)/4.;
                    //double testint2 = fmod(testint,1);
                    //if(testint2==0)
                    //if (oetl == kTRUE)
                    //{
                    SumE = SumE + DRHit->GetEdep();
                    SumNCeren = SumNCeren + DRHit->GetNCeren();
                    //}
                } // End loop over Hits
            }
        }// End loop over Hit collections

        double corrE = SumE * scorr / 1000.;
        double corrC = SumNCeren*ccorr;
        double ratio = corrC / corrE;
        double drcEcorr = corrE / fcorr(ratio);

        SumE2 = SumE2 + (drcEcorr - MeanEdep)*(drcEcorr - MeanEdep);
        SumNCeren2 = SumNCeren2 + (corrC - MeanNCeren)*(corrC - MeanNCeren);

    } // end loop over evts.
    VarEdep = SumE2 / (double(nevent) - 1);
    VarNCeren = SumNCeren2 / (double(nevent) - 1);

    Double_t sigmaEdep = sqrt(VarEdep);
    Double_t sigmaNCeren = sqrt(VarNCeren);

    piontop->cd();

    sprintf(histoname, "drcpiEdep%fGeV", Ein);
    sprintf(histotitle, "pion total Energy deposition (dr corrected) (Ein %f GeV)", Ein);
    drcpionEdep.push_back(new TH1F(histoname, histotitle, 250, MeanEdep - 5 * sigmaEdep, MeanEdep + 5 * sigmaEdep));
    //drcpionEdep.push_back(new TH1F(histoname, histotitle, 250, 0, 1.1*Ein));
    drcpionEdep[index]->GetXaxis()->SetTitle("E[GeV]");
    drcpionEdep[index]->GetYaxis()->SetTitle("Evts.");
    drcpionEdep[index]->SetFillColor(7);
    drcpionEdep[index]->SetLineColor(4);

    //Third loop over events(Fill histograms)
    for (Int_t i = 0; i < nevent; i++) {
        Double_t SumE = 0.0;
        Int_t SumNCeren = 0;
        T->GetEntry(i);
        map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
        map<G4String, vector<G4VHit*> >::iterator hciter;
        for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
            vector<G4VHit*> hits = (*hciter).second;
            G4int NbHits = hits.size();
            vector<string> y = split((*hciter).first, '_');
            string Classname = y[1];
            if (Classname == "DRTSCalorimeter") {
                for (G4int ii = 0; ii < NbHits; ii++) {
                    DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
                    G4ThreeVector pos = DRHit->GetPos();
                    double z = pos.getZ();
                    int num = 100;
                    int cellsize = 25;
                    int layer = z / cellsize + num / 2;
                    //Bool_t oetl = OddEv(layer);
                    //double testint = (layer - 1.)/3.; 
                    //double testint = (layer - 2.)/3.;
                    //double testint = (layer - 3.)/3.;
                    //double testint = (layer - 1.)/4.;
                    //double testint2 = fmod(testint,1);
                    //if(testint2==0)
                    //if (oetl == kTRUE)
                    //{
                    SumE = SumE + DRHit->GetEdep();
                    SumNCeren = SumNCeren + DRHit->GetNCeren();
                    //}
                } // End loop over Hits
            }
        }// End loop over Hit collections
        double corrE = SumE * scorr / 1000.;
        double corrC = SumNCeren*ccorr;
        double ratio = corrC / corrE;
        double drcEcorr = corrE / fcorr(ratio);
        drcpionEdep[index]->Fill(drcEcorr);

    } // end loop over evts.
    piontop->cd();
    drcpionEdepptr.push_back(drcpionEdep[index]->Fit("gaus", "S"));
    drcpionEdep[index]->Write();
    return;
}

void getEconst(unsigned int index) {
    //
    // subroutine to correct and analyze the response of electrons in a dual
    // readout calorimeter. It gets the constants to convert 
    // scint and cerenkov response to the energy of the 
    // incoming monoenergetic electron.
    // Input parameters:
    // index of electron input file 
    //
    electrontop->cd();
    char histoname[50];
    char histotitle[100];
    double Ein = eenergies[index];

    TFile f(efnames[index].c_str());
    /*
    Event *event = new Event();
    TTree *T = (TTree*) f.Get("T");
    T->SetBranchAddress("event.", &event);
*/
    
    
    
        RunHeader *runh = new RunHeader();
        Event *event = new Event();
        TTree *T = (TTree*) f.Get("Events");
        T->SetBranchAddress("event.", &event);
        TBranch* fevtbranch = T->GetBranch("event.");
        TTree *Trh = (TTree*) f.Get("Runheader");
        Trh->SetBranchAddress("RunHeader.", &runh);
        TBranch* frunhbranch = Trh->GetBranch("RunHeader.");
        frunhbranch->GetEntry(0);
      //  runh->Print();
      //  Int_t nevent = T->GetEntries();
    //    G4cout << "Nr. of Events: " << nevent << G4endl;


    
    
    
    Int_t nevent = T->GetEntries();

    Double_t MeanNCeren = 0.0;
    Double_t MeanEdep = 0.0;
    Double_t SumEp = 0.0;
    Int_t SumNCerenp = 0.0;

    //First loop over events (Mean)
    for (Int_t i = 0; i < nevent; i++) {
        Double_t SumE = 0.0;
        Int_t SumNCeren = 0;

        T->GetEntry(i);
        map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
        map<G4String, vector<G4VHit*> >::iterator hciter;
        for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
            vector<G4VHit*> hits = (*hciter).second;
            G4int NbHits = hits.size();
            vector<string> y = split((*hciter).first, '_');
            string Classname = y[1];
            if (Classname == "DRTSCalorimeter") {
                for (G4int ii = 0; ii < NbHits; ii++) {
                    DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
                    G4ThreeVector pos = DRHit->GetPos();
                    double z = pos.getZ();
                    int num = 100;
                    int cellsize = 25;
                    int layer = z / cellsize + num / 2;
                    //  Bool_t oetl = OddEv(layer);
                    //double testint = (layer - 1.)/3.; 
                    //double testint = (layer - 2.)/3.;
                    //double testint = (layer - 3.)/3.;
                    //double testint = (layer - 1.)/4.;
                    //double testint2 = fmod(testint,1);
                    //if(testint2==0)
                    //if (oetl == kTRUE)
                    //{
                    SumE = SumE + (DRHit->GetEdep() / 1000.);
                    SumNCeren = SumNCeren + DRHit->GetNCeren();
                    //}
                } // End loop over Hits
            }
        }// End loop over Hit collections
        SumEp = SumEp + SumE;
        SumNCerenp = SumNCerenp + SumNCeren;
    } // end loop over evts.

    electrontop->cd();

    MeanNCeren = SumNCerenp / double(nevent);
    MeanEdep = SumEp / double(nevent);
    G4cout << "MeanEdep: " << MeanEdep << G4endl;
    G4cout << "MeanNCeren:  " << MeanNCeren << G4endl;

    Double_t SumNCeren2 = 0.0;
    Double_t SumE2 = 0.0;

    //Second loop over events (Variation calculation)
    for (Int_t i = 0; i < nevent; i++) {
        Int_t SumNCeren = 0;
        Double_t SumE = 0.0;

        T->GetEntry(i);
        map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
        map<G4String, vector<G4VHit*> >::iterator hciter;
        for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
            vector<G4VHit*> hits = (*hciter).second;
            G4int NbHits = hits.size();
            vector<string> y = split((*hciter).first, '_');
            string Classname = y[1];
            if (Classname == "DRTSCalorimeter") {
                for (G4int ii = 0; ii < NbHits; ii++) {
                    DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
                    G4ThreeVector pos = DRHit->GetPos();
                    double z = pos.getZ();
                    int num = 100;
                    int cellsize = 25;
                    int layer = z / cellsize + num / 2;
                    //Bool_t oetl = OddEv(layer);
                    //double testint = (layer - 1.)/3.; 
                    //double testint = (layer - 2.)/3.;
                    //double testint = (layer - 3.)/3.;
                    //double testint = (layer - 1.)/4.;
                    //double testint2 = fmod(testint,1);
                    //if(testint2==0)
                    //if (oetl == kTRUE)
                    //{
                    SumNCeren = SumNCeren + DRHit->GetNCeren();
                    SumE = SumE + (DRHit->GetEdep() / 1000.);
                    //}
                } // End loop over Hits
            }
        }// End loop over Hit collections

        SumE2 = SumE2 + (SumE - MeanEdep)*(SumE - MeanEdep);
        SumNCeren2 = SumNCeren2 + (SumNCeren - MeanNCeren)*(SumNCeren - MeanNCeren);
    } // end loop over evts.

    Double_t VarEdep = SumE2 / (double(nevent) - 1.);
    Double_t VarNCeren = SumNCeren2 / (double(nevent) - 1.);

    Double_t sigmaEdep = sqrt(VarEdep);
    Double_t sigmaNCeren = sqrt(VarNCeren);

    G4cout << "sigmaEdep:  " << sigmaEdep << G4endl;
    G4cout << "sigmaNCeren:  " << sigmaNCeren << G4endl;

    sprintf(histoname, "eCeren%fGeV", eenergies[index]);
    sprintf(histotitle, "total number of C photons (Ein = %f GeV)", eenergies[index]);
    electronCeren.push_back(new TH1F(histoname, histotitle, 500, MeanNCeren - 5 * sigmaNCeren, MeanNCeren + 5 * sigmaNCeren));
    electronCeren[index]->GetXaxis()->SetTitle("Nr of cer. Photons");
    electronCeren[index]->GetYaxis()->SetTitle("Nr of Evts.");
    electronCeren[index]->SetFillColor(7);
    electronCeren[index]->SetLineColor(4);

    sprintf(histoname, "eEdep%fGeV", eenergies[index]);
    sprintf(histotitle, "uncorrected total Energy deposition electrons Ein = %f GeV", eenergies[index]);
    electronEdep.push_back(new TH1F(histoname, histotitle, 250, MeanEdep - 5 * sigmaEdep, MeanEdep + 5 * sigmaEdep));
    //electronEdep.push_back(new TH1F(histoname, histotitle, 250, 0, 1.1*Ein));
    electronEdep[index]->GetXaxis()->SetTitle("E[GeV]");
    electronEdep[index]->GetYaxis()->SetTitle("Nr of Evts.");
    electronEdep[index]->SetFillColor(7);
    electronEdep[index]->SetLineColor(4);

    //Third loop over events(Fill histograms)
    for (Int_t i = 0; i < nevent; i++) {
        Int_t SumNCeren = 0;
        Double_t SumE = 0.0;
        T->GetEntry(i);
        map<G4String, vector<G4VHit*> >* hcmap = event->GetHCMap();
        map<G4String, vector<G4VHit*> >::iterator hciter;
        for (hciter = hcmap->begin(); hciter != hcmap->end(); hciter++) {
            vector<G4VHit*> hits = (*hciter).second;
            G4int NbHits = hits.size();
            vector<string> y = split((*hciter).first, '_');
            string Classname = y[1];
            if (Classname == "DRTSCalorimeter") {
                for (G4int ii = 0; ii < NbHits; ii++) {
                    DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
                    G4ThreeVector pos = DRHit->GetPos();
                    double z = pos.getZ();
                    int num = 100;
                    int cellsize = 25;
                    int layer = z / cellsize + num / 2;
                    //Bool_t oetl = OddEv(layer);
                    //double testint = (layer - 1.)/3.; 
                    //double testint = (layer - 2.)/3.;
                    //double testint = (layer - 3.)/3.;
                    //double testint = (layer - 1.)/4.;
                    //double testint2 = fmod(testint,1);
                    //if(testint2==0)
                    //if(oetl == kTRUE)
                    //{
                    SumNCeren = SumNCeren + DRHit->GetNCeren();
                    SumE = SumE + DRHit->GetEdep();
                    //}
                } // End loop over Hits
            }
        }// End loop over Hit collections
        electronCeren[index]->Fill(SumNCeren);
        electronEdep[index]->Fill((SumE) / 1000.);

    } // end loop over evts.

    electrontop->cd();
    electronEdepptr.push_back(electronEdep[index]->Fit("gaus", "S"));
    electronEdep[index]->Write();
    electronCerenptr.push_back(electronCeren[index]->Fit("gaus", "S"));
    electronCeren[index]->Write();
    return;
}

void calE() {
    vector<double> Econst;
    const unsigned int n1 = efnames.size();
    double x[100];
    double xerr[100];
    double E[100];
    double Eerr[100];
    double C[100];
    double Cerr[100];
    for (unsigned int i = 0; i < efnames.size(); i++) {
        x[i] = eenergies[i];
        xerr[i] = 0.0;
        getEconst(i);
        E[i] = electronEdepptr[i]->Parameter(1);
        Eerr[i] = electronEdepptr[i]->ParError(1);
        C[i] = electronCerenptr[i]->Parameter(1);
        Cerr[i] = electronCerenptr[i]->ParError(1);
    }
    finaltop->cd();
    TGraphErrors *gr = new TGraphErrors(n1, x, E, xerr, Eerr);
    gr->SetLineColor(4);
    gr->SetLineWidth(4);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(10);
    gr->SetTitle("electron scint. response");
    gr->SetName("electron scint. response");
    gr->GetXaxis()->SetTitle("Ein [GeV]");
    gr->GetYaxis()->SetTitle("Emeas. [GeV]");
    TFitResultPtr fitpt = gr->Fit("pol1", "S");
    scorr = 1 / (fitpt->Parameter(1));
    gr->Write();
    TGraphErrors *grc = new TGraphErrors(n1, x, C, xerr, Cerr);
    grc->SetLineColor(3);
    grc->SetLineWidth(4);
    grc->SetMarkerColor(3);
    grc->SetMarkerStyle(11);
    grc->SetTitle("electron C response");
    grc->SetName("electron C response");
    grc->GetXaxis()->SetTitle("Ein [GeV]");
    grc->GetYaxis()->SetTitle("Nceren");
    TFitResultPtr fitptc = grc->Fit("pol1", "S");
    ccorr = 1 / (fitptc->Parameter(1));
    grc->Write();

    //  vector<double> Eresp;
    const unsigned int n2 = pifnames.size();
    Double_t Epi[100];
    Double_t Einpi[100];
    Double_t oversqrtEinpi[100];
    Double_t Eresrel[100];
    Double_t Erespi[100];
    Double_t Edrcpi[100];
    Double_t Edrcrespi[100];
    for (unsigned int i = 0; i < n2; i++) {
        pidist("Pions", i);
    }
    for (unsigned int i = 0; i < n1; i++) {
        pidist("Electrons", i);
    }

    TF1* fpow = new TF1("fpow", "1-pow(x/[0],[1]-1)", 0, 100);
    Double_t E0 = 1.;
    Double_t m = 0.85;
    fpow->SetParameter(0, E0);
    fpow->SetParameter(1, m);

    TMultiGraph *mg = new TMultiGraph();

    double y1[100];
    double x1[100];
    double ey1[100];
    double ex1[100];
    int nevents = 1000;

    for (unsigned int i = 0; i < n2; i++) {
        y1[i] = FEMhisto[i]->GetMean();
        x1[i] = pienergies[i];
        ey1[i] = FEMhisto[i]->GetRMS();
        ex1[i] = 1 / sqrt(nevents);
    }

    finaltop->cd();

    TGraph *gr1 = new TGraph(n2, x1, y1);
    gr1->SetLineColor(2);
    gr1->SetLineWidth(1);
    gr1->SetMarkerColor(4);
    gr1->SetMarkerStyle(21);
    gr1->SetName("fem");
    gr1->SetTitle("<f_{em}>");
    gr1->Fit("fpow");
    gr1->Write();
    mg->Add(gr1);

    double y2[100];
    double x2[100];
    double ey2[100];
    double ex2[100];

    for (unsigned int i = 0; i < n2; i++) {
        y2[i] = Pi0histo[i]->GetMean();
        x2[i] = pienergies[i];
        ey2[i] = Pi0histo[i]->GetRMS();
        ex2[i] = 1 / sqrt(nevents);
    }

    TGraph *gr2 = new TGraph(n2, x2, y2);
    gr2->SetLineColor(3);
    gr2->SetLineWidth(1);
    gr2->SetMarkerColor(5);
    gr2->SetMarkerStyle(21);
    gr2->SetName("fPi0");
    gr2->SetTitle("<f#pi^{0}>");
    gr2->Fit("fpow");
    gr2->Write();
    mg->Add(gr2);

    mg->SetTitle("Electromagnetic Energy Fractions");
    mg->SetName("<fem> & <fPi0>");
    mg->Write();

    TMultiGraph *mg2 = new TMultiGraph();

    TGraph *gr12 = new TGraph(n2, x1, ey1);
    gr12->SetLineColor(2);
    gr12->SetLineWidth(1);
    gr12->SetMarkerColor(4);
    gr12->SetMarkerStyle(21);
    gr12->SetName("#sigma_{f_{em}}");
    gr12->SetTitle("#sigma_{f_{em}}");
    gr12->Write();
    mg2->Add(gr12);

    TGraph *gr22 = new TGraph(n2, x2, ey2);
    gr22->SetLineColor(3);
    gr22->SetLineWidth(1);
    gr22->SetMarkerColor(5);
    gr22->SetMarkerStyle(21);
    gr22->SetName("#sigma_{f_{#pi^{0}}}");
    gr22->SetTitle("#sigma_{f_{#pi^{0}}}");
    gr22->Write();
    mg2->Add(gr22);

    mg2->SetTitle("#sigma for electromagnetic and #pi^{0} fractions");
    mg2->SetName("sigma fem & fPi0");
    mg2->Write();

    double y3[100];
    double x3[100];
    double ey3[100];
    double ex3[100];
    TCanvas *c3 = new TCanvas("c3", "hadronic response", 200, 10, 1000, 800);
    c3->SetGridx();
    c3->SetGridy();
    TMultiGraph *rmg = new TMultiGraph();
    rmg->SetMinimum(0.0);
    rmg->SetMaximum(1.1);
    rmg->SetTitle("hadronic response;Ekin [GeV];hadronic response");
    //    for (vector<TGraphErrors*>::iterator ii = RatioGraphs.begin(); ii != RatioGraphs.end(); ++ii) {
    //        rmg->Add((*ii));
    //    }

    for (unsigned int i = 0; i < n2; i++) {
        y3[i] = FEMhshisto[i]->GetMean();
        x3[i] = pienergies[i];
        ey3[i] = FEMhshisto[i]->GetRMS();
        ex3[i] = 1 / sqrt(nevents);
    }
    TGraph *hsfem = new TGraph(n2, x3, y3);
    hsfem->SetLineColor(2);
    hsfem->SetLineWidth(1);
    hsfem->SetMarkerColor(2);
    hsfem->SetMarkerStyle(21);
    hsfem->SetTitle("hs(FEM)");
    hsfem->SetName("FEMhshisto");
    rmg->Add(hsfem);
    double y4[100];
    double x4[100];
    double ey4[100];
    double ex4[100];

    for (unsigned int i = 0; i < n2; i++) {
        y4[i] = FEMhchisto[i]->GetMean();
        x4[i] = pienergies[i];
        ey4[i] = FEMhchisto[i]->GetRMS();
        ex4[i] = 1 / sqrt(nevents);
    }
    TGraph *hcfem = new TGraph(n2, x4, y4);
    hcfem->SetLineColor(3);
    hcfem->SetLineWidth(1);
    hcfem->SetMarkerColor(3);
    hcfem->SetMarkerStyle(21);
    hcfem->SetTitle("hc(FEM)");
    hcfem->SetName("FEMhchisto");
    rmg->Add(hcfem);

    double y5[100];
    double x5[100];
    double ey5[100];
    double ex5[100];

    for (unsigned int i = 0; i < n2; i++) {
        y5[i] = Pi0hshisto[i]->GetMean();
        x5[i] = pienergies[i];
        ey5[i] = Pi0hshisto[i]->GetRMS();
        ex5[i] = 1 / sqrt(nevents);
    }
    TGraph *hsPi0 = new TGraph(n2, x5, y5);
    hsPi0->SetLineColor(4);
    hsPi0->SetLineWidth(1);
    hsPi0->SetMarkerColor(4);
    hsPi0->SetMarkerStyle(21);
    hsPi0->SetTitle("hs(Pi0)");
    hsPi0->SetName("Pi0hshisto");
    rmg->Add(hsPi0);
    double y6[100];
    double x6[100];
    double ey6[100];
    double ex6[100];

    for (unsigned int i = 0; i < n2; i++) {
        y6[i] = Pi0hchisto[i]->GetMean();
        x6[i] = pienergies[i];
        ey6[i] = Pi0hchisto[i]->GetRMS();
        ex6[i] = 1 / sqrt(nevents);
    }
    TGraph *hcPi0 = new TGraph(n2, x6, y6);
    hcPi0->SetLineColor(6);
    hcPi0->SetLineWidth(1);
    hcPi0->SetMarkerColor(6);
    hcPi0->SetMarkerStyle(21);
    hcPi0->SetTitle("hc(Pi0)");
    hcPi0->SetName("Pi0hchisto");
    rmg->Add(hcPi0);
    rmg->Draw("apl");
    TLegend *leg3 = c3->BuildLegend(.6, .3, 0.85, .55);
    leg3->Draw();
    c3->Write();
    double y9[100];
    double x9[100];

    for (unsigned int i = 0; i < n1; i++) {
        y9[i] = electronCeren[i]->GetMean();
        x9[i] = eenergies[i];
    }
    TGraph *NCphot = new TGraph(n1, x9, y9);
    NCphot->SetLineColor(4);
    NCphot->SetLineWidth(3);
    NCphot->SetMarkerColor(6);
    NCphot->SetTitle("Number of Cherenkov photons");
    NCphot->SetName("NCphot");
    NCphot->GetXaxis()->SetTitle("Energy[GeV]");
    NCphot->GetYaxis()->SetTitle("Number of C photons");
    NCphot->Write();

    double y10[100];
    double x10[100];
    double ey10[100];
    double ex10[100];

    for (unsigned int i = 0; i < n1; i++) {
        y10[i] = electronEdepptr[i]->Parameter(1) / eenergies[i];
        ey10[i] = electronEdepptr[i]->ParError(1);
        x10[i] = eenergies[i];
        ex10[i] = 0.001 * x10[i];

    }
    TGraphErrors *evisein = new TGraphErrors(n1, x10, y10, ex10, ey10);
    evisein->SetLineColor(kRed);
    evisein->SetLineWidth(3);
    evisein->SetMarkerColor(kRed);
    evisein->SetMarkerStyle(10);
    evisein->SetTitle("Electron response to scintillation light as a function of Ein");
    evisein->SetName("S/E");
    evisein->GetXaxis()->SetTitle("Energy[GeV]");
    evisein->GetYaxis()->SetTitle("S/E");
    evisein->Write();

    TFitResultPtr fit_ptr = hprof2->Fit("pol3", "S");
    const unsigned int npar = fit_ptr->NPar();
    TFitResultPtr slice_ptr;
    //
    vector<double> vmeanslice;
    vector<double> vmeansliceerr;
    vector<double> veslice;
    vector<double> vesliceerr;
    for (unsigned int ii = 0; ii < nslices; ii++) {
        if (hslice[ii]->GetEntries() > 100) {
            hslice[ii]->Fit("gaus", "S");
            slice_ptr = hslice[ii]->Fit("gaus", "S");
            veslice.push_back((ii + 0.5) * slicewidth);
            vesliceerr.push_back(0.0);
            vmeanslice.push_back(slice_ptr->Parameter(1));
            vmeansliceerr.push_back(slice_ptr->ParError(1));
        }
    }
    const unsigned int nrofslices = veslice.size();
    double meanslice[100];
    double meansliceerr[100];
    double eslice[100];
    double esliceerr[100];
    for (unsigned int ii = 0; ii < nrofslices; ii++) {
        meanslice[ii] = vmeanslice[ii];
        meansliceerr[ii] = vmeansliceerr[ii];
        eslice[ii] = veslice[ii];
        esliceerr[ii] = vesliceerr[ii];
    }
    finaltop->cd();
    TGraphErrors *grslices = new TGraphErrors(nrofslices, eslice, meanslice, esliceerr, meansliceerr);
    grslices->SetName("DR correction function");
    grslices->SetFillColor(1);
    grslices->SetLineColor(4);
    grslices->SetLineWidth(0);
    grslices->SetMarkerColor(4);
    grslices->SetMarkerStyle(21);
    grslices->SetTitle("Dual Readout correction function");
    grslices->GetXaxis()->SetTitle("C/S");
    grslices->GetYaxis()->SetTitle("S/Ein");
    TF1 *pol3 = new TF1("pol3", "pol3", 0., 1.1);
    pol3->SetFillColor(19);
    pol3->SetFillStyle(0);
    pol3->SetLineColor(2);
    pol3->SetLineWidth(3);
    pol3->SetLineStyle(2);
    TFitResultPtr fitptsl = grslices->Fit(pol3, "S");
    finaltop->cd();
    grslices->Write();
    for (unsigned int i = 0; i < npar; i++) {
        corrf[i] = fitptsl->Parameter(i);
    }

    Double_t Epierr[100];
    Double_t Edrcpierr[100];
    Double_t Erespierr[100];
    Double_t Edrcrespierr[100];
    Double_t Einpierr[100];

    for (unsigned int i = 0; i < n2; i++) {
        CorrectedE(i);
        Epi[i] = cpionEdepptr[i]->Parameter(1);
        Epierr[i] = cpionEdepptr[i]->ParError(1);
        Erespi[i] = cpionEdepptr[i]->Parameter(1);
        Erespierr[i] = cpionEdepptr[i]->ParError(1);
        Edrcpi[i] = drcpionEdepptr[i]->Parameter(1);
        Edrcpierr[i] = drcpionEdepptr[i]->ParError(1);
        Edrcrespi[i] = drcpionEdepptr[i]->Parameter(2);
        Edrcrespierr[i] = drcpionEdepptr[i]->ParError(2);
        Einpi[i] = pienergies[i];
        Einpierr[i] = 0.001 * Einpi[i];
        oversqrtEinpi[i] = 1. / TMath::Sqrt(pienergies[i]);
        Eresrel[i] = (100. * Edrcrespi[i]) / Einpi[i];
    }

    TGraphErrors *grepi = new TGraphErrors(n2, Einpi, Epi, Einpierr, Epierr);
    TGraphErrors *greinpi = new TGraphErrors(n2, Einpi, Einpi, Einpierr, Einpierr);
    TGraphErrors *gredrcpi = new TGraphErrors(n2, Einpi, Edrcpi, Einpierr, Edrcpierr);

    grepi->SetLineColor(4);
    grepi->SetLineWidth(4);
    grepi->SetMarkerColor(4);
    grepi->SetMarkerStyle(10);
    grepi->SetTitle(" scint. response");
    grepi->SetName("Energy response");
    grepi->GetXaxis()->SetTitle("Ein [GeV]");
    grepi->GetYaxis()->SetTitle("Emeas. [GeV]");
    gredrcpi->SetLineColor(6);
    gredrcpi->SetLineWidth(4);
    gredrcpi->SetMarkerColor(6);
    gredrcpi->SetMarkerStyle(11);
    finaltop->cd();
    grepi->Write();
    greinpi->Write();
    gredrcpi->Write();
    //
    TGraphErrors *grerespi = new TGraphErrors(n2, Einpi, Erespi, Einpierr, Erespierr);
    grerespi->SetLineColor(4);
    grerespi->SetLineWidth(4);
    grerespi->SetMarkerColor(4);
    grerespi->SetMarkerStyle(10);
    grerespi->SetTitle(" scint. resolution");
    grerespi->SetName("Energy resolution");
    grerespi->GetYaxis()->SetTitle("energy resolution [GeV]");
    grerespi->GetXaxis()->SetTitle("Ein. [GeV]");
    grerespi->Write();
    // 
    TGraphErrors *greresdrcpi = new TGraphErrors(n2, Einpi, Edrcrespi, Einpierr, Edrcpierr);
    greresdrcpi->SetName("DR corrected Energy resolution");
    greresdrcpi->SetTitle("Energy resolution (dual read out cor.)");
    greresdrcpi->SetLineColor(6);
    greresdrcpi->SetLineWidth(4);
    greresdrcpi->SetMarkerColor(6);
    greresdrcpi->SetMarkerStyle(11);
    greresdrcpi->GetYaxis()->SetTitle("energy resolution [GeV]");
    greresdrcpi->GetXaxis()->SetTitle("Ein. [GeV]");
    greresdrcpi->Write();
    //
    TGraph *erelres = new TGraph(n2, Einpi, Eresrel);
    erelres->SetName("rel. DR corrected Energy resolution");
    erelres->SetTitle("rel. Energy resolution (dual read out cor.)");
    erelres->SetLineColor(6);
    erelres->SetLineWidth(4);
    erelres->SetMarkerColor(6);
    erelres->SetMarkerStyle(11);
    erelres->GetYaxis()->SetTitle("rel energy resolution [percent]");
    erelres->GetXaxis()->SetTitle("Ein. [GeV]");
    erelres->Write();
    //
    TGraph *erelresvssqrte = new TGraph(n2, oversqrtEinpi, Eresrel);
    erelresvssqrte ->SetName("rel. DR corrected Energy resolution vs 1/sqrt(E)");
    erelresvssqrte ->SetTitle("rel. Energy resolution (dual read out cor.) vs 1/sqrt(e)");
    erelresvssqrte ->SetLineColor(6);
    erelresvssqrte ->SetLineWidth(4);
    erelresvssqrte ->SetMarkerColor(6);
    erelresvssqrte ->SetMarkerStyle(11);
    erelresvssqrte ->GetYaxis()->SetTitle("rel energy resolution [percent]");
    erelresvssqrte ->GetXaxis()->SetTitle("1/sqrt(Ein.)");
    erelresvssqrte ->Write();

}

int main(int argc, char** argv) {
    TSystem ts;
    gSystem->Load("libCintex");
    gSystem->Load("libClassesDict");
    ROOT::Cintex::Cintex::Enable();
    if (argc < 2) G4cout << "Missing name of the file to read!" << G4endl;
    string fname(argv[1]);
    init(fname.c_str());
    calE();
    finaltop->cd();
    banana2->Write();
    hprof2->Write();
    for (unsigned int ii = 0; ii < nslices; ii++) {
        hslice[ii]->Write();
    }
}
