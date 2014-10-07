{
  TFile *f = new TFile("p_15GeV_analysis.root");
  TCanvas *c1 = new TCanvas("c1", "c1",520,89,700,500);
  c1->SetLogx();
  c1->SetLogy();
  TH1F * h1 = (TH1F*) f->Get("scdEdxweighted");
  TH1F * h2 = (TH1F*) f->Get("scdEdxunweighted");
  TH1F * h3 = (TH1F*) f->Get("sclowdEdxunweighted");
  TH1F * h4 = (TH1F*) f->Get("sclowdEdxweighted");
  h1->Add(h1, h4, 1., 1.);
  h2->Add(h2, h3, 1., 1.);
  h1->Divide(h1, h2, 1., 1.,"b");
  h1->GetXaxis()->SetRange(100.,4000.);
  h1->SetMinimum(0.3);
  h1->SetMinimum(30.);
  h1->GetYaxis()->SetTitle("-dE/dx  [MeV/g cm^{2}]");
  h1->GetXaxis()->SetTitle("proton momentum [MeV/c]");
  Int_t ci;   // for color index setting
  ci = TColor::GetColor("#000099");
  h1->SetLineColor(ci);
  h1->SetMarkerColor(2);
  h1->SetMarkerStyle(2);
  h1->Draw("EBAR");
}
