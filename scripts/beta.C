//**********************************************************************************************
//
// Script to demonstarte fit o an asymetric function 
// (e.g. the fem distribution in hadronic showers
// To run it:  
// TFile *f = TFile::Open("../../CaTS-build/PbF2_FTFPBERT_pi-_50.0GeV_Hits.root");
// .L beta.C++
// beta(f);
//
// Author: Hans Wenzel
//
//**********************************************************************************************

#include <stdio.h>
#include <math.h>
#include "TVirtualFitter.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <cassert>
#include <iostream>
#include "TMath.h"
#include "TFile.h"
TCanvas* beta( TFile* f)
{
  TH1F * h1 = (TH1F*) f->Get("histos/sumEdepEMfrac");
  TCanvas *cst = new TCanvas("beta", "beta", 200, 10, 1000, 800);
  h1->SetTitle("");
  h1->GetXaxis()->SetTitle("f_em");
  h1->GetYaxis()->SetTitle("Nr of Events");
  TF1 *fbeta = new TF1("fbeta", "[3]*TMath::BetaDist(x-[2], [0], [1])", 0, 1);
  fbeta->SetParameters(3, 9,1.,0.5);
  fbeta->SetParNames("alpha","beta","deltax","Norm");
  fbeta->SetLineColor(kBlue);
  h1->Fit("fbeta","R");
  TF1 *flog = new TF1("flog", "[3]*TMath::LogNormal(x, [0], [1], [2])", 0, 5);
  flog->SetParameters(0.35, 0.25, 0.3,20.);
  flog->SetParNames("sigma","theta","Norm");
  flog->SetLineColor(kRed);
  
  h1->Fit("flog","R+");
  TLegend *legend = new TLegend(.15, .7, 0.4, .8);
  legend->AddEntry(fbeta, "displaced beta function");
  legend->AddEntry(flog, "log normal function");
  legend->Draw();

  return cst;
 
}
