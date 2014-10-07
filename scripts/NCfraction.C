void NCfraction (){
  const Int_t n1= 7;
  Double_t Ein[n1]= {1.0,2.0,5.0,10.0,20.0,50.0,100.0};  
  Double_t positron[n1]={11.7462,15.475,20.1397,22.8428,25.2136,29.9286,28.2492};
  Double_t electron[n1]={37.6317,44.0028,50.2964,52.1483,54.1586,55.8085,57.0366};
  Double_t Kpm[n1]={0.00115921,0.0468801,0.154134,0.358981,0.471125,0.58034,0.513579};
  Double_t mupm[n1]={0.511014,0.343874,0.374269,0.2577276,0.250603,0.212083,0.184986};
  Double_t pip[n1]={0.939494,3.18113,5.25974,7.08843,6.66464,6.17484,5.4922};
  Double_t pim[n1]={47.1543,34.1998,20.6959,14.302,10.6919,7.98568,6.49317};
  Double_t pdthe3[n1]={2.01612,2.74311,3.07101,2.98372,2.50702,2.24199,1.96989};
  Double_t sigchi[n1]={2.7e-5,0.00740945,0.00884415,0.00267166,0.00758934,0.00681457,0.0082059};
  Double_t anti[n1]={0.,0.,0.,0.00763936,0.0349533,0.0612475,0.052239};
  TCanvas *c = new TCanvas("c", "response", 200, 10, 1000, 800);
  TMultiGraph *mg = new TMultiGraph();

  mg->SetTitle("Composition of cerenkov response in #pi^{-} showers ;Ekin [GeV];[%]");
  c->SetGrid();
  TGraph *gr_electron = new TGraph(n1,Ein,electron);
  gr_electron->SetLineColor(4);
  gr_electron->SetLineWidth(1);
  gr_electron->SetMarkerColor(4);
  gr_electron->SetMarkerStyle(20);
  gr_electron->SetMarkerSize(1.2);
  gr_electron->SetTitle("e^{-}");
  mg->Add( gr_electron);

  TGraph *gr_positron = new TGraph(n1,Ein,positron);
  gr_positron->SetLineColor(2);
  gr_positron->SetLineWidth(1);
  gr_positron->SetMarkerColor(2);
  gr_positron->SetMarkerStyle(20);
  gr_positron->SetMarkerSize(1.2);
  gr_positron->SetTitle("e^{+}");
  mg->Add( gr_positron);

  TGraph *gr_pim = new TGraph(n1,Ein,pim);
  gr_pim->SetLineColor(7);
  gr_pim->SetLineWidth(1);
  gr_pim->SetMarkerColor(7);
  gr_pim->SetMarkerStyle(22);
  gr_pim->SetMarkerSize(1.2);
  gr_pim->SetTitle("#pi^{-}");
  mg->Add( gr_pim);

  TGraph *gr_pip = new TGraph(n1,Ein,pip);
  gr_pip->SetLineColor(6);
  gr_pip->SetLineWidth(1);
  gr_pip->SetMarkerColor(6);
  gr_pip->SetMarkerStyle(21);
  gr_pip->SetMarkerSize(1.2);
  gr_pip->SetTitle("#pi^{+}");
  mg->Add( gr_pip);

  TGraph *gr_Kpm = new TGraph(n1,Ein,Kpm);
  gr_Kpm->SetLineColor(3);
  gr_Kpm->SetLineWidth(1);
  gr_Kpm->SetMarkerColor(3);
  gr_Kpm->SetMarkerStyle(20);
  gr_Kpm->SetMarkerSize(1.2);
  gr_Kpm->SetTitle("K^{+-}");
  mg->Add( gr_Kpm);

  TGraph *gr_mupm = new TGraph(n1,Ein,mupm);
  gr_mupm->SetLineColor(5);
  gr_mupm->SetLineWidth(1);
  gr_mupm->SetMarkerColor(5);
  gr_mupm->SetMarkerStyle(21);
  gr_mupm->SetMarkerSize(1.2);
  gr_mupm->SetTitle("#mu^{+-}");
  mg->Add( gr_mupm);

  TGraph *gr_pdthe3 = new TGraph(n1,Ein,pdthe3);
  gr_pdthe3->SetLineColor(8);
  gr_pdthe3->SetLineWidth(1);
  gr_pdthe3->SetMarkerColor(8);
  gr_pdthe3->SetMarkerStyle(20);
  gr_pdthe3->SetMarkerSize(1.2);
  gr_pdthe3->SetTitle("p,d,t,He3");
  mg->Add( gr_pdthe3);

  TGraph *gr_sigchi = new TGraph(n1,Ein,sigchi);
  gr_sigchi->SetLineColor(9);
  gr_sigchi->SetLineWidth(1);
  gr_sigchi->SetMarkerColor(9);
  gr_sigchi->SetMarkerStyle(23);
  gr_sigchi->SetMarkerSize(1.2);
  gr_sigchi->SetTitle("#Sigma,#Chi");
  mg->Add( gr_sigchi);

  TGraph *gr_anti = new TGraph(n1,Ein,anti);
  gr_anti->SetLineColor(46);
  gr_anti->SetLineWidth(1);
  gr_anti->SetMarkerColor(46);
  gr_anti->SetMarkerStyle(21);
  gr_anti->SetMarkerSize(1.2);
  gr_anti->SetTitle("Antiparticles");
  mg->Add( gr_anti);
  mg->SetMinimum(0.);
  mg->SetMaximum(65.);
  mg->Draw("apl");
  TLegend *leg = c->BuildLegend(.6, .6, 0.85, .85);
  leg->Draw();
}
