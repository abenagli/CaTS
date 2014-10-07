//**********************************************************************************************
//
// Script to display beta spectrum of charged particles in a hadronic shower
// To run it:  
// TFile *f = TFile::Open("../../CaTS-build/data/PbWO_FTFPBERT_pi-_100.0GeV_analysis.root");
// .L betaspec.C+
// betaspec(f);
//
// Author: Hans Wenzel
//
//**********************************************************************************************

#include <sstream>
#include <string>
#include "TFile.h"
#include "THStack.h"
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

TCanvas* betaspec( TFile* f)
{
  TH1F * h = (TH1F*) f->Get("beta");
  h->SetFillColor(kBlue-7);
  THStack *hs = new THStack("hs","beta of charge particles in e- shower");
  TH1F * hmum = (TH1F*) f->Get("Particles/mu-");
  if (hmum!=0)
    {
      hmum->SetFillColor(kMagenta);
      hs->Add(hmum);
    }
  TH1F * hmup = (TH1F*) f->Get("Particles/mu+");
  if (hmup!=0)
    {
      hmup->SetFillColor(kOrange);
      hs->Add(hmup);
    }
 TH1F * halpha = (TH1F*) f->Get("Particles/alpha");
  if (halpha!=0)
    {
      halpha->SetFillColor(kCyan-5);
      hs->Add(halpha);
    }
  TH1F * h4 = (TH1F*) f->Get("Particles/pi+");
  if (h4!=0)
    {
      h4->SetFillColor(kGreen);
      hs->Add(h4);
    }

  TH1F * h3 = (TH1F*) f->Get("Particles/pi-");
  if (h3!=0)
    {
      h3->SetFillColor(kYellow);
      hs->Add(h3);
    }
  TH1F * h5 = (TH1F*) f->Get("Particles/proton");
  if (h5!=0)
    {
      h5->SetFillColor(kCyan);
      hs->Add(h5);
    }
  TH1F * h1 = (TH1F*) f->Get("Particles/e+");
  if (h1!=0)
    {
      h1->SetFillColor(kBlue);
      hs->Add(h1);
    }
  TH1F * h2 = (TH1F*) f->Get("Particles/e-");
  if  (h2!=0)
    {
      h2->SetFillColor(kRed);
      hs->Add(h2);
    }
  TCanvas *cst = new TCanvas("beta", "beta", 200, 10, 1000, 800);
  cst->SetLogy();
  h->SetTitle("");
  gStyle->SetOptStat(0);
  h->Draw();
  h->SetMinimum(1);
  h->GetXaxis()->SetTitle("#beta");
  h->GetYaxis()->SetTitle("Nr of charged particles");
  hs->Draw("same");
  TLegend *legend = new TLegend(0.15,0.6,0.3,0.89);
  legend->AddEntry(h,"nuclear fragments","f");
  if (h1 != 0)   legend->AddEntry(h1,"e^{+}","f");
  if (h2 != 0)   legend->AddEntry(h2,"e^{-}","f");
  if (h5 != 0)   legend->AddEntry(h5,"p","f");
  if (hmum != 0) legend->AddEntry(hmum,"#mu^{-}","f");
  if (hmup != 0) legend->AddEntry(hmup,"#mu^{+}","f");
  if (halpha!=0) legend->AddEntry(halpha,"#alpha","f");
  if (h4!=0)     legend->AddEntry(h4,"#pi^{+}","f");
  if (h3!=0)     legend->AddEntry(h3,"#pi^{-}","f");
  legend->Draw();
  return cst;

}
