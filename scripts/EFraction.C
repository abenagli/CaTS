void EFraction (){
  const Int_t n1= 7;
  Double_t Ein[n1]= {1.0,2.0,5.0,10.0,20.0,50.0,100.0};  
  Double_t fragment[n1]={2.63064,2.56792,2.5284,2.21382,1.83636,1.64084};
  Double_t dthe3[n1]={0.634382,1.01121,1.63965,1.24904,0.982694,0.890259,0.791409};
  Double_t alpha[n1]={1.04224,1.55261,1.40288,1.13924,0.963478,0.854013,0.771514};
  Double_t proton[n1]={26.8095,27.4247,28.0612,24.7721,21.0491,18.8348,16.7872};
  Double_t neutron[n1]={1.26606,1.21668,1.1042,0.997219,0.856071,0.772535,0.695081};
  Double_t positron[n1]={6.88568,8.50101,10.5646,12.704,15.0721,16.6683,18.0369};
  Double_t electron[n1]={31.2011,34.0171,36.8554,40.3973,44.2021,46.8557,49.1316};
  Double_t Kpm[n1]={0.00668517,0.0942877,0.156929,0.242917,0.386265,0.4321129,0.399839};
  Double_t mupm[n1]={0.309868,0.139467,0.138826,0.0718473,0.0677658,0.052096,0.0352328};
  Double_t pip[n1]={0.904588,2.34092,3.31886,4.43833,4.4204,4.23413,3.86905};
  Double_t pim[n1]={25.959,18.5026,11.4362,8.78995,6.97546,5.44015,4.58038};
  Double_t gamma[n1]={2.30205,2.47342,2.63234,2.80419,2.97644,3.10488,3.21214};
  Double_t sigpm[n1]={0.00730525,0.0141167,0.0151327,0.0112362,0.0127214,0.0144003,0.0127149};
  TCanvas *c = new TCanvas("c", "response", 200, 10, 1000, 800);
  TMultiGraph *mg = new TMultiGraph();

  mg->SetTitle("Composition of Ionization response in #pi^{-} showers ;Ekin [GeV];[%]");
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

  TGraph *gr_dthe3 = new TGraph(n1,Ein,dthe3);
  gr_dthe3->SetLineColor(8);
  gr_dthe3->SetLineWidth(1);
  gr_dthe3->SetMarkerColor(8);
  gr_dthe3->SetMarkerStyle(20);
  gr_dthe3->SetMarkerSize(1.2);
  gr_dthe3->SetTitle("d,t,He3");
  mg->Add( gr_dthe3);

  TGraph *gr_sigpm = new TGraph(n1,Ein,sigpm);
  gr_sigpm->SetLineColor(9);
  gr_sigpm->SetLineWidth(1);
  gr_sigpm->SetMarkerColor(9);
  gr_sigpm->SetMarkerStyle(23);
  gr_sigpm->SetMarkerSize(1.2);
  gr_sigpm->SetTitle("#Sigma^{+_}");
  mg->Add( gr_sigpm);

  TGraph *gr_proton = new TGraph(n1,Ein,proton);
  gr_proton->SetLineColor(46);
  gr_proton->SetLineWidth(1);
  gr_proton->SetMarkerColor(46);
  gr_proton->SetMarkerStyle(21);
  gr_proton->SetMarkerSize(1.2);
  gr_proton->SetTitle("p");
  mg->Add( gr_proton);

  TGraph *gr_alpha = new TGraph(n1,Ein,alpha);
  gr_alpha->SetLineColor(47);
  gr_alpha->SetLineWidth(1);
  gr_alpha->SetMarkerColor(47);
  gr_alpha->SetMarkerStyle(22);
  gr_alpha->SetMarkerSize(1.2);
  gr_alpha->SetTitle("#alpha");
  mg->Add( gr_alpha);

  TGraph *gr_gamma = new TGraph(n1,Ein,gamma);
  gr_gamma->SetLineColor(48);
  gr_gamma->SetLineWidth(1);
  gr_gamma->SetMarkerColor(48);
  gr_gamma->SetMarkerStyle(23);
  gr_gamma->SetMarkerSize(1.2);
  gr_gamma->SetTitle("#gamma");
  mg->Add( gr_gamma);

  TGraph *gr_fragment = new TGraph(n1,Ein,fragment);
  gr_fragment->SetLineColor(49);
  gr_fragment->SetLineWidth(1);
  gr_fragment->SetMarkerColor(49);
  gr_fragment->SetMarkerStyle(24);
  gr_fragment->SetMarkerSize(1.2);
  gr_fragment->SetTitle("Nuclear fragments");
  mg->Add( gr_fragment);

  TGraph *gr_neutron = new TGraph(n1,Ein,neutron);
  gr_neutron->SetLineColor(50);
  gr_neutron->SetLineWidth(1);
  gr_neutron->SetMarkerColor(50);
  gr_neutron->SetMarkerStyle(25);
  gr_neutron->SetMarkerSize(1.2);
  gr_neutron->SetTitle("n");
  mg->Add( gr_neutron);

  mg->SetMinimum(0.);
  mg->SetMaximum(50.);
  mg->Draw("apl");
  TLegend *leg = c->BuildLegend(.6, .45, 0.85, .85);
  leg->Draw();
}
