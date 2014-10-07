//**********************************************************************************************
//
// Script to display how different particles contribute energy shower
// To run it:  
// TFile *f = TFile::Open("../../CaTS-build/data/PbWO_FTFPBERT_pi-_100.0GeV_analysis.root");
// .L NCcont.C+
// NCcont(f);
//
// Author: Hans Wenzel
//
//**********************************************************************************************

#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include "TFile.h"
#include "THStack.h"
#include "TVirtualFitter.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPie.h"

TCanvas* NCcont(TFile* file) {
    int index = 0;
    double vals[100];
    std::string partnames[100] = {
        "deuteron",
        "triton",
        "He3",
        "proton",
        "e+",
        "e-",
        "mu+",
        "mu-",
        "pi+",
        "pi-",
        "kaon+",
        "kaon-",
        "sigma+",
        "sigma-",
        "xi-",
        "anti_sigma+",
        "anti_sigma-",
        "anti_proton",
        "anti_xi-",
        "anti_omega-"
    };
    std::string labels[100] = {
        "deuteron",
        "triton",
        "He3",
        "p",
        "e^{+}",
        "e^{-}",
        "#mu^{+}",
        "#mu^{-}",
        "#pi^{+}",
        "#pi^{-}",
        "K^{+}",
        "K^{-}",
        "#Sigma^{+}",
        "#Sigma^{-}",
        "#Xi^{-}",
        "#bar{#Sigma}^{+}",
        "#bar{#Sigma}^{-}",        
        "#bar{p}",
        "#bar{#Xi}^{-}",
        "#bar{#Omega}^{-}"
    };

    const char* comblabels[100];
    int colors[100] = {1, 3, 2, 4, 5, 6, 7, 8, 9, 46, 40, 30, 38, 41, 29, 16, 32, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
    std::string dirname = "NCerencont/";
    std::string prefix = "NCeren_";
    std::vector< std::vector<int> > grouping;
    std::vector<int> temp;
    temp.push_back(0);
    temp.push_back(1);
    temp.push_back(2);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(3);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(4);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(5);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(6);
    temp.push_back(7);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(8);
    temp.push_back(9);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(10);
    temp.push_back(11);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(12);
    temp.push_back(13);
    temp.push_back(14);    
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(15);
    temp.push_back(16);
    temp.push_back(17);
    temp.push_back(18);
    temp.push_back(19);    
    grouping.push_back(temp);
    TH1F * htmp;
    TH1F * hsum;
    std::string histoname;
    THStack *hs = new THStack("hs", "Cerenkov Photon  contribution of particles in shower");
    TLegend *legend = new TLegend(0.5, 0.5, 0.9, 0.9);
    TCanvas *chisto = new TCanvas("chisto", "Cerenkov photon contribution by particle", 700, 700);
    chisto->cd();
    for (unsigned int ii = 0; ii < grouping.size(); ii++) {
        std::string label = labels[grouping[ii][0]];
        histoname = dirname + prefix + partnames[grouping[ii][0]];
        hsum = (TH1F*) file->Get(histoname.c_str());
        for (unsigned int jj = 1; jj < grouping[ii].size(); jj++) {
            histoname = dirname + prefix + partnames[grouping[ii][jj]];
            htmp = (TH1F*) file->Get(histoname.c_str());
            hsum->Add(htmp);
            label = label + "," + labels[grouping[ii][jj]];
        }
        if (hsum != 0) {
            hsum->SetFillColor(colors[index]);
            vals[index] = hsum->GetMean();
            comblabels[index] = label.c_str();
            index++;
            hs->Add(hsum);
            legend->AddEntry(hsum, label.c_str(), "f");
        }
    }
    gStyle->SetOptStat(0);
    hs->Draw();
    legend->Draw();
    for (int i = 0; i < index; i++) {
        cout << comblabels[i] << "  Fraction:  " << vals[i] << endl;
    }
    TCanvas *cpie = new TCanvas("cpie", "TPie test", 700, 700);
    cpie->cd();
    TPie *pie1 = new TPie("pie1", "Cerenkov photon contribution by particle in shower", index, vals, colors, comblabels);
    //  TLegend *pieleg = pie1->MakeLegend();
    pie1->Draw("3d");
    legend->Draw();
    return cpie;

}
