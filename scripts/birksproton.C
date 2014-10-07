{  
  TFile *f = new TFile("p_15GeV_analysis.root");
  TCanvas *c1 = new TCanvas("c1", "c1",520,89,700,500);
  TH1F * h1 = (TH1F*) f->Get("scBirksweighted");
  TH1F * h2 = (TH1F*) f->Get("scBirksunweighted");
  TH1F * h3 = (TH1F*) f->Get("sclowBirksweighted");
  TH1F * h4 = (TH1F*) f->Get("sclowBirksunweighted");
  TH1F * h5 = (TH1F*) h1->Clone();
  TH1F * h6 = (TH1F*) h2->Clone();
  h1->Add(h1, h3, 1., 1.);
  h2->Add(h2, h4, 1., 1.);
  h1->Divide(h1, h2, 1., 1.,"b");
  h1->GetYaxis()->SetTitle("average Birks Factor");
  h1->GetXaxis()->SetTitle("proton kinetic Energy [MeV]");
  h1->SetLineColor(2);
  h1->SetLineWidth(1);
  h1->Draw();
  //TH1F * h5 = (TH1F*) f->Get("Birksweighted");
  //TH1F * h6 = (TH1F*) f->Get("Birksunweighted");
  h5->Divide(h5, h6, 1., 1.,"b");
  h5->SetLineColor(3);
  h5->SetLineWidth(1);
  h5->Draw("SAME");
  TLegend *legend = new TLegend(.5, .2, .85, .4);
 legend->AddEntry(h1,  "with low momentum protons");
 legend->AddEntry(h5, "no low momentum protons");
 legend->Draw();

}

