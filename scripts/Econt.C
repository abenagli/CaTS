//**********************************************************************************************
//
// Script to display how different particles contribute energy shower
// To run it:  
// TFile *f = TFile::Open("../../CaTS-build/data/PbWO_FTFPBERT_pi-_100.0GeV_analysis.root");
// .L Econt.C+
// Econt(f);
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

TCanvas* Econt(TFile* file) {
    int index = 0;
    double vals[100];
    std::string partnames[100] = {
        "Fragment",
        "He3",
        "alpha",
        "deuteron",
        "triton",
        "neutron",
        "p_ev",
        "p_he",
        "p_sp",
        "e+",
        "e-",
        "gamma",
        "mu+",
        "mu-",
        "pi+",
        "pi-",
        "lambda",
        "sigma-",
        "sigma+",
        "kaon+",
        "kaon-",
        "kaon0S",
        "kaon0L"
    };
    std::string labels[100] = {
        "nuclear Fragments",
        "He^{3}",
        "#alpha",
        "deuteron",
        "triton",
        "neutron",
        "p_ev",
        "p_he",
        "p_sp",
        "e^{+}",
        "e^{-}",
        "#gamma",
        "#mu^{+}",
        "#mu^{-}",
        "#pi^{+}",
        "#pi^{-}",
        "#Lambda",
        "#Sigma^{-}",
        "#Sigma^{+}",
        "K^{+}",
        "K^{-}",
        "K_{S}",
        "K_{L}"
    };
    const char* comblabels[100];
    int colors[100] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 30, 40, 46, 38, 41, 29, 16, 32, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
    //std::string dirname = "Energycont/";
    std::string dirname = "histos/CalorimeterVolume_DRTSCalorimeter/";
    //std::string prefix = "E_";
    std::string prefix = "";
    std::vector< std::vector<int> > grouping;
    std::vector<int> temp;
    temp.push_back(0);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(1);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(2);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(3);
    temp.push_back(4);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(5);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(6);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(7);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(8);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(9);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(10);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(11);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(12);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(13);
    grouping.push_back(temp);
    temp.clear();
    temp.push_back(14);
    temp.push_back(15);
    temp.push_back(16);
    temp.push_back(17);
    temp.push_back(18);
    temp.push_back(19);
    temp.push_back(20);
    grouping.push_back(temp);
    /*    temp.clear();
        temp.push_back(17);
        grouping.push_back(temp);
        temp.clear();
        temp.push_back(18);
        grouping.push_back(temp);
        temp.clear();
        temp.push_back(19);
        temp.push_back(20);
        grouping.push_back(temp);
     */
    TH1F * htmp;
    TH1F * hsum;
    std::string histoname;
    THStack *hs = new THStack("hs", "Energy contribution of particles in pi- shower");\
    TLegend *legend = new TLegend(0.6, 0.4, 0.9, 0.9);
    TCanvas *chisto = new TCanvas("chisto", "Energy contribution by perticle", 700, 700);
    chisto->cd();
    //chisto->cd();
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
    //return chisto;

    TCanvas *cpie = new TCanvas("cpie", "TPie test", 700, 700);
    cpie->cd();
    TPie *pie1 = new TPie("pie1", "energy distribution by particle ", index, vals, colors, comblabels);
    TLegend *pieleg = pie1->MakeLegend();
    pie1->Draw("3d");
    pieleg->Draw();
    cout << comblabels << endl;
    //    TLegend *pielegend = new TLegend(0.65, 0.4, 0.8, 0.9);
    return cpie;

}
