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
#include "TMultiGraph.h"
#include "TLegend.h"
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
map<string, vector< std::pair<double, TH1F* > > >histograms;
map<string, vector< std::pair<double, TH1F* > > >drhistograms;
map<string, vector< std::pair<double, TH1F* > > >Cerenhistograms;
map<string, vector< std::pair<double, TH1F* > > >Ratiohistograms;

vector<TGraphErrors*>EdepGraphs;
vector<TGraphErrors*>CerenGraphs;
vector<TGraphErrors*>RatioGraphs;
vector<TGraphErrors*>cEdepGraphs;
vector<TGraphErrors*>drcEdepGraphs;
vector<TGraphErrors*>cCerenGraphs;
TFile* outfile;
double scorr, ccorr;

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
//
// Dual read out correction function:
//
// Ein = S/fcorr(C/S)
// 

double fcorr(double ratio) {
    //    double corrf[4];
    double x = ratio;
    double fc = 0;
    if (ratio < 0.4) {
        x = 0.4;
    }
    if (x > 1) {
        fc = 1.;
        return fc;
    }
    // hardwired to BGO sheet calorimeter
    double p0 = -0.0442095;
    double p1 = 3.30036;
    double p2 = -3.85846;
    double p3 = 1.60686;
    //    double p0 = corrf[0];
    //    double p1 = corrf[1];
    //    double p2 = corrf[2];
    //    double p3 = corrf[3];
    fc = p0 + p1 * x + p2 * x * x + p3 * x * x*x;

    return fc;
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
        cout << "couldn't Open input file! Exiting" << endl;
        exit(1);
    }
}

void getresponse(const char *filename) {
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
    Event *event = new Event();
    TTree *T = (TTree*) fo.Get("T");
    T->SetBranchAddress("event.", &event);
    Int_t nevent = T->GetEntries();
    G4cout << " Nr. of Events:  " << nevent << G4endl;
    Double_t MeanNCeren = 0.0;
    Double_t MeanEdep = 0.0;
    Double_t drMeanEdep = 0.0;
    Double_t SumEp = 0.0;
    Double_t drSumEp = 0.0;
    Double_t SumNCerenp = 0.0;

    //First loop over events to get the mean of the distribution
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
            if (Classname == "DRCalorimeter") {
                for (G4int ii = 0; ii < NbHits; ii++) {
                    DRCalorimeterHit* DRHit = dynamic_cast<DRCalorimeterHit*> (hits[ii]);
                    SumE = SumE + (DRHit->GetEdep() / 1000);
                    SumNCeren = SumNCeren + DRHit->GetNCeren();
                } // End loop over Hits
            }
        }// End loop over Hit collections

        double corrE = SumE*scorr;
        double corrC = SumNCeren*ccorr;
        double ratio = corrC / corrE;
        double drcorrE = corrE / fcorr(ratio);
        SumNCerenp = SumNCerenp + corrC;
        SumEp = SumEp + corrE;
        drSumEp = drSumEp + drcorrE;
    }// End loop over Events

    MeanNCeren = SumNCerenp / double(nevent);
    MeanEdep = SumEp / double(nevent);
    drMeanEdep = drSumEp / double(nevent);
    Double_t SumNCeren2 = 0.0;
    Double_t SumE2 = 0.0;
    Double_t drSumE2 = 0.0;
    //
    //Second loop over events(Variation Calculation)
    //
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
            if (Classname == "DRCalorimeter") {
                for (G4int ii = 0; ii < NbHits; ii++) {
                    DRCalorimeterHit* DRHit = dynamic_cast<DRCalorimeterHit*> (hits[ii]);
                    SumE = SumE + (DRHit->GetEdep() / 1000);
                    SumNCeren = SumNCeren + DRHit->GetNCeren();
                } // End loop over Hits
            }
        }// End loop over Hit collections

        double corrE = SumE*scorr;
        double corrC = SumNCeren*ccorr;
        double ratio = corrC / corrE;
        double drcorrE = corrE / fcorr(ratio);
        SumE2 = SumE2 + (corrE - MeanEdep)*(corrE - MeanEdep);
        SumNCeren2 = SumNCeren2 + (corrC - MeanNCeren)*(corrC - MeanNCeren);
        drSumE2 = drSumE2 + (drcorrE - drMeanEdep)*(drcorrE - drMeanEdep);
    }// End loop over Events

    Double_t VarEdep = SumE2 / (double(nevent) - 1.);
    Double_t drVarEdep = drSumE2 / (double(nevent) - 1.);
    Double_t VarNCeren = SumNCeren2 / (double(nevent) - 1.);
    Double_t sigmaEdep = sqrt(VarEdep);
    Double_t drsigmaEdep = sqrt(drVarEdep);
    Double_t sigmaNCeren = sqrt(VarNCeren);
    G4cout << particle << " (" << energy << " GeV), Corrected  Mean Edep:  "
            << MeanEdep << "  +/- "
            << sigmaEdep << "  (GeV)" << G4endl;
    G4cout << particle << " (" << energy << " GeV), dr Corrected  Mean Edep:  "
            << drMeanEdep << "  +/- "
            << drsigmaEdep << "  (GeV)" << G4endl;
    G4cout << particle << " (" << energy << " GeV), Corrected Mean NCeren:  "
            << MeanNCeren << " +/- "
            << sigmaNCeren << "  (GeV)" << G4endl;

    char histoname[50];
    char histotitle[100];
    sprintf(histoname, "Edep%fGeV", Ein);
    sprintf(histotitle, "Edep corrected energy Ein = %fGeV", Ein);
    if (histograms.find(particle) == histograms.end()) {
        // not found
        std::pair<double, TH1F* >tmppair(Ein, new TH1F(histoname, histotitle, 250, MeanEdep - 5 * sigmaEdep, MeanEdep + 5 * sigmaEdep));
        vector<std::pair<double, TH1F*> >tempvec;
        tempvec.push_back(tmppair);
        histograms.insert(std::pair<string, vector<std::pair<double, TH1F*> > >(particle, tempvec));
    } else {
        // found
        histograms[particle].push_back(std::pair<double, TH1F* >(Ein, new TH1F(histoname, histotitle, 250, MeanEdep - 5 * sigmaEdep, MeanEdep + 5 * sigmaEdep)));
    }
    sprintf(histoname, "drEdep%fGeV", Ein);
    sprintf(histotitle, "Edep dr corrected energy Ein = %fGeV", Ein);
    if (drhistograms.find(particle) == drhistograms.end()) {
        // not found
        std::pair<double, TH1F* >tmppair(Ein, new TH1F(histoname, histotitle, 250, drMeanEdep - 5 * drsigmaEdep, drMeanEdep + 5 * drsigmaEdep));
        vector<std::pair<double, TH1F*> >tempvec;
        tempvec.push_back(tmppair);
        drhistograms.insert(std::pair<string, vector<std::pair<double, TH1F*> > >(particle, tempvec));
    } else {
        // found
        drhistograms[particle].push_back(std::pair<double, TH1F* >(Ein, new TH1F(histoname, histotitle, 250, drMeanEdep - 5 * drsigmaEdep, drMeanEdep + 5 * drsigmaEdep)));
    }

    sprintf(histoname, "NCeren%fGeV", Ein);
    sprintf(histotitle, "NCeren  corrected energy Ein = %fGeV", Ein);
    if (Cerenhistograms.find(particle) == Cerenhistograms.end()) {
        // not found
        std::pair<double, TH1F* >tmppair(Ein, new TH1F(histoname, histotitle, 250, MeanNCeren - 5 * sigmaNCeren, MeanNCeren + 5 * sigmaNCeren));
        vector<std::pair<double, TH1F*> >tempvec;
        tempvec.push_back(tmppair);
        Cerenhistograms.insert(std::pair<string, vector<std::pair<double, TH1F*> > >(particle, tempvec));
    } else {
        // found       
        Cerenhistograms[particle].push_back(std::pair<double, TH1F* >(Ein, new TH1F(histoname, histotitle, 250, MeanNCeren - 5 * sigmaNCeren, MeanNCeren + 5 * sigmaNCeren)));
    }

    sprintf(histoname, "Ratio%fGeV", Ein);
    sprintf(histotitle, "C/ ratioEin = %fGeV", Ein);
    if (Ratiohistograms.find(particle) == Ratiohistograms.end()) {
        // not found
        std::pair<double, TH1F* >tmppair(Ein, new TH1F(histoname, histotitle, 100, 0., 1.1));
        vector<std::pair<double, TH1F*> >tempvec;
        tempvec.push_back(tmppair);
        Ratiohistograms.insert(std::pair<string, vector<std::pair<double, TH1F*> > >(particle, tempvec));
    } else {
        // found       
        Ratiohistograms[particle].push_back(std::pair<double, TH1F* >(Ein, new TH1F(histoname, histotitle, 100, 0., 1.1)));
    }


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
            if (Classname == "DRCalorimeter") {
                for (G4int ii = 0; ii < NbHits; ii++) {
                    DRCalorimeterHit* DRHit = dynamic_cast<DRCalorimeterHit*> (hits[ii]);
                    SumE = SumE + (DRHit->GetEdep());
                    SumNCeren = SumNCeren + DRHit->GetNCeren();
                    //}
                } // End loop over Hits
            }
        }// End loop over Hit collections

        double corrE = SumE * scorr / 1000.;
        double corrC = SumNCeren*ccorr;
        double ratio = corrC / corrE;
        double drcorrE = corrE / fcorr(ratio);
        histograms[particle].back().second->Fill(corrE);
        drhistograms[particle].back().second->Fill(drcorrE);
        Cerenhistograms[particle].back().second->Fill(corrC);
        Ratiohistograms[particle].back().second->Fill(corrC / corrE);
    }// End loop over Events
    G4cout << " nr of bytes written:  " << outfile->Write() << G4endl;
    return;
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
    // correction factor from electron calibration.
    scorr = 1. / 0.997244; // For now hardwired (BGO FTFP_BERT)
    ccorr = 1. / 65553.4;
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
            getresponse(p);
        }
    }
    Double_t energies[100];
    Double_t errenergies[100];

    Double_t meanv[100];
    Double_t cmeanv[100];
    Double_t drcmeanv[100];
    Double_t errmeanv[100];
    Double_t cerrmeanv[100];
    Double_t drcerrmeanv[100];
    //
    TMultiGraph *Edepmg = new TMultiGraph();
    Edepmg->SetMinimum(0.5);
    Edepmg->SetMaximum(2.5);
    Edepmg->SetTitle("relative Energy response;Ekin [GeV];S/Ekin");
    int Color = 2;
    int Style = 20;

    for (map<string, vector< std::pair<double, TH1F* > > >::iterator ii = histograms.begin(); ii != histograms.end(); ++ii) {
        unsigned int count = 0;
        string particle = (*ii).first;
        double deltaE = 0.0;
        if (particle == "antiproton") {
            deltaE = 2. * 0.9383;
        } else if (particle == "pi+") {
            deltaE = 0.13957;
        } else if (particle == "pi-") {
            deltaE = 0.13957;
        } else if (particle == "kaon+") {
            deltaE = 0.493677;
        } else if (particle == "kaon-") {
            deltaE = 0.493677;
        } else {
            deltaE = 0.0;
        }
        for (vector< std::pair<double, TH1F* > >::iterator jj = (*ii).second.begin(); jj != (*ii).second.end(); ++jj) {
            TFitResultPtr fitpt = (*jj).second->Fit("gaus", "S");
            meanv[count] = fitpt->Parameter(1) / (*jj).first;
            cmeanv[count] = fitpt->Parameter(1) / ((*jj).first + deltaE);
            errmeanv[count] = fitpt->ParError(1) / (*jj).first;
            cerrmeanv[count] = fitpt->ParError(1) / ((*jj).first + deltaE);
            energies[count] = (*jj).first;
            errenergies[count] = 0.001;
            count++;
        }

        TGraphErrors * tg = new TGraphErrors(count, energies, meanv, errenergies, errmeanv);
        tg->SetLineColor(Color);
        tg->SetLineWidth(1);
        tg->SetMarkerColor(Color);
        tg->SetMarkerStyle(Style);
        tg->SetTitle((*ii).first.c_str());
        EdepGraphs.push_back(tg);
        TGraphErrors * ctg = new TGraphErrors(count, energies, cmeanv, errenergies, cerrmeanv);
        ctg->SetLineColor(Color);
        ctg->SetLineWidth(1);
        ctg->SetMarkerColor(Color);
        ctg->SetMarkerStyle(Style);
        ctg->SetTitle((*ii).first.c_str());
        cEdepGraphs.push_back(ctg);
        Color++;
        if (Color == 5 || Color == 10) Color++; // avoid yellow and white 
        Style++;
    }
    Color = 2;
    Style = 20;

    for (map<string, vector< std::pair<double, TH1F* > > >::iterator ii = drhistograms.begin(); ii != drhistograms.end(); ++ii) {
        unsigned int count = 0;
        string particle = (*ii).first;
        double deltaE = 0.0;
        if (particle == "antiproton") {
            deltaE = 2. * 0.9383;
        } else if (particle == "pi+") {
            deltaE = 0.13957;
        } else if (particle == "pi-") {
            deltaE = 0.13957;
        } else if (particle == "kaon+") {
            deltaE = 0.493677;
        } else if (particle == "kaon-") {
            deltaE = 0.493677;
        } else {
            deltaE = 0.0;
        }
        for (vector< std::pair<double, TH1F* > >::iterator jj = (*ii).second.begin(); jj != (*ii).second.end(); ++jj) {
            TFitResultPtr fitpt = (*jj).second->Fit("gaus", "S");
            drcmeanv[count] = fitpt->Parameter(1) / ((*jj).first + deltaE);
            drcerrmeanv[count] = fitpt->ParError(1) / ((*jj).first + deltaE);
            energies[count] = (*jj).first;
            errenergies[count] = 0.001;
            count++;
        }

        TGraphErrors * drctg = new TGraphErrors(count, energies, drcmeanv, errenergies, drcerrmeanv);
        drctg->SetLineColor(Color);
        drctg->SetLineWidth(1);
        drctg->SetMarkerColor(Color);
        drctg->SetMarkerStyle(Style);
        drctg->SetTitle((*ii).first.c_str());
        drcEdepGraphs.push_back(drctg);
        Color++;
        if (Color == 5 || Color == 10) Color++; // avoid yellow and white 
        Style++;
    }
    Color = 2;
    Style = 20;
    for (map<string, vector< std::pair<double, TH1F* > > >::iterator ii = Cerenhistograms.begin(); ii != Cerenhistograms.end(); ++ii) {
        unsigned int count = 0;
        string particle = (*ii).first;
        double deltaE = 0.0;
        if (particle == "antiproton") {
            deltaE = 2. * 0.9383;
        } else if (particle == "pi+") {
            deltaE = 0.13957;
        } else if (particle == "pi-") {
            deltaE = 0.13957;
        } else if (particle == "kaon+") {
            deltaE = 0.493677;
        } else if (particle == "kaon-") {
            deltaE = 0.493677;
        } else {
            deltaE = 0.0;
        }
        for (vector< std::pair<double, TH1F* > >::iterator jj = (*ii).second.begin(); jj != (*ii).second.end(); ++jj) {
            TFitResultPtr fitpt = (*jj).second->Fit("gaus", "S");
            meanv[count] = fitpt->Parameter(1) / (*jj).first;
            cmeanv[count] = fitpt->Parameter(1) / ((*jj).first + deltaE);
            errmeanv[count] = fitpt->ParError(1) / (*jj).first;
            cerrmeanv[count] = fitpt->ParError(1) / ((*jj).first + deltaE);
            errenergies[count] = 0.001;
            count++;
        }
        TGraphErrors * tg = new TGraphErrors(count, energies, meanv, errenergies, errmeanv);
        tg->SetLineColor(Color);
        tg->SetLineWidth(1);
        tg->SetMarkerColor(Color);
        tg->SetMarkerStyle(Style);
        tg->SetTitle((*ii).first.c_str());
        CerenGraphs.push_back(tg);
        //
        TGraphErrors * ctg = new TGraphErrors(count, energies, cmeanv, errenergies, cerrmeanv);
        ctg->SetLineColor(Color);
        ctg->SetLineWidth(1);
        ctg->SetMarkerColor(Color);
        ctg->SetMarkerStyle(Style);
        ctg->SetTitle((*ii).first.c_str());
        cCerenGraphs.push_back(ctg);

        Color++;
        if (Color == 5 || Color == 10) Color++; // avoid yellow and white 
        Style++;
    }
    Color = 2;
    Style = 20;
    for (map<string, vector< std::pair<double, TH1F* > > >::iterator ii = Ratiohistograms.begin(); ii != Ratiohistograms.end(); ++ii) {
        unsigned int count = 0;
        for (vector< std::pair<double, TH1F* > >::iterator jj = (*ii).second.begin(); jj != (*ii).second.end(); ++jj) {
            TFitResultPtr fitpt = (*jj).second->Fit("gaus", "S");
            meanv[count] = fitpt->Parameter(1);
            errmeanv[count] = fitpt->ParError(1);
            energies[count] = (*jj).first;
            errenergies[count] = 0.001;
            count++;
        }
        TGraphErrors * tg = new TGraphErrors(count, energies, meanv, errenergies, errmeanv);
        tg->SetLineColor(Color);
        tg->SetLineWidth(1);
        tg->SetMarkerColor(Color);
        tg->SetMarkerStyle(Style);
        tg->SetTitle((*ii).first.c_str());
        RatioGraphs.push_back(tg);
        Color++;
        if (Color == 5 || Color == 10) Color++; // avoid yellow and white 
        Style++;
    }
    outfile->mkdir("top");
    outfile->cd("top");
    TCanvas *c = new TCanvas("c", "response", 200, 10, 1000, 800);
    c->SetGridx();
    c->SetGridy();
    TMultiGraph *mg = new TMultiGraph();
    mg->SetMinimum(0.5);
    mg->SetMaximum(2.5);
    mg->SetTitle("relative Energy response;Ekin [GeV];S/Ekin");
    for (vector<TGraphErrors*>::iterator ii = EdepGraphs.begin(); ii != EdepGraphs.end(); ++ii) {
        mg->Add((*ii));
    }
    mg->SetMaximum(2.5);
    mg->SetMinimum(0.8);
    mg->Draw("apl");
    TLegend *leg = c->BuildLegend(.6, .6, 0.85, .85);
    leg->Draw();
    c->Write();
    //
    outfile->cd("top");
    TCanvas *c2 = new TCanvas("c2", "Cerenkov response", 200, 10, 1000, 800);
    c2->SetGridx();
    c2->SetGridy();
    TMultiGraph *cmg = new TMultiGraph();
    cmg->SetMinimum(0.2);
    cmg->SetMaximum(2.);
    cmg->SetTitle("Cerenkov relative Energy response;Ekin [GeV];C/Ekin");
    for (vector<TGraphErrors*>::iterator ii = CerenGraphs.begin(); ii != CerenGraphs.end(); ++ii) {
        cmg->Add((*ii));
    }
    cmg->Draw("apl");
    TLegend *leg2 = c2->BuildLegend(.6, .6, 0.85, .85);
    leg2->Draw();
    c2->Write();
    //    
    outfile->cd("top");
    TCanvas *cc = new TCanvas("cc", "response (invariant mass considered)", 200, 10, 1000, 800);
    cc->SetGridx();
    cc->SetGridy();
    TMultiGraph *c_mg = new TMultiGraph();
    c_mg->SetMinimum(0.5);
    c_mg->SetMaximum(2.5);
    c_mg->SetTitle("relative Energy response;Ekin [GeV];S/Ein");
    for (vector<TGraphErrors*>::iterator ii = cEdepGraphs.begin(); ii != cEdepGraphs.end(); ++ii) {
        c_mg->Add((*ii));
    }
    c_mg->SetMaximum(1.05);
    c_mg->SetMinimum(0.8);
    c_mg->Draw("apl");
    TLegend *cleg = cc->BuildLegend(.6, .6, 0.85, .85);
    cleg->Draw();
    cc->Write();
    //    
    outfile->cd("top");
    TCanvas *drcc = new TCanvas("drcc", "dr corrected response (invariant mass considered)", 200, 10, 1000, 800);
    drcc->SetGridx();
    drcc->SetGridy();
    TMultiGraph *drc_mg = new TMultiGraph();
    drc_mg->SetMinimum(0.8);
    drc_mg->SetMaximum(1.05);
    drc_mg->SetTitle("relative Energy response;Ekin [GeV];(dr corrected) S/Ein");
    for (vector<TGraphErrors*>::iterator ii = drcEdepGraphs.begin(); ii != drcEdepGraphs.end(); ++ii) {
        drc_mg->Add((*ii));
    }
    drc_mg->SetMaximum(1.1);
    drc_mg->SetMinimum(0.9);
    drc_mg->Draw("apl");
    TLegend *drcleg = drcc->BuildLegend(.6, .15, 0.85, .4);
    drcleg->Draw();
    drcc->Write();
    //
    outfile->cd("top");
    TCanvas *cc2 = new TCanvas("cc2", "Cerenkov response (invariant mass considered) ", 200, 10, 1000, 800);
    cc2->SetGridx();
    cc2->SetGridy();
    TMultiGraph *ccmg = new TMultiGraph();
    ccmg->SetMinimum(0.3);
    ccmg->SetMaximum(1.1);
    ccmg->SetTitle("Cerenkov relative Energy response;Ekin [GeV];C/Ein");
    for (vector<TGraphErrors*>::iterator ii = cCerenGraphs.begin(); ii != cCerenGraphs.end(); ++ii) {
        ccmg->Add((*ii));
    }
    ccmg->Draw("apl");
    TLegend *cleg2 = cc2->BuildLegend(.6, .6, 0.85, .85);
    cleg2->Draw();
    cc2->Write();
    //   
    outfile->cd("top");
    TCanvas *c3 = new TCanvas("c3", "NCeren/Edep  ratio", 200, 10, 1000, 800);
    c3->SetGridx();
    c3->SetGridy();
    TMultiGraph *rmg = new TMultiGraph();
    rmg->SetMinimum(0.4);
    rmg->SetMaximum(1.1);
    rmg->SetTitle("C/S ratio;Ekin [GeV];C/S");
    for (vector<TGraphErrors*>::iterator ii = RatioGraphs.begin(); ii != RatioGraphs.end(); ++ii) {
        rmg->Add((*ii));
    }
    rmg->Draw("apl");
    TLegend *leg3 = c3->BuildLegend(.6, .6, 0.85, .85);
    leg3->Draw();
    c3->Write();


}
