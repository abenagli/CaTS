void pi0(){
  const Int_t n1= 7;
  const Int_t n2= 5;
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
  Double_t pi0[n1]={11.6125,16.6218,22.8034,29.7709,37.4228,42.726,46.6687};
  Double_t npi0[n1]={1.14,1.317,2.029,3.476,6.78,14.99,27.86};
  Double_t npiplus[n1]={1.111,1.261,1.764,3.178,5.97,14.3,26.22};
  Double_t np[n1]={5.046,10.111,24.73,43.29,73.99,169.4,305.4};
  Double_t nn[n2]={70.43,127.8,281.1,509.2,864.8};
  Double_t em[n1];
  Double_t ratio[n1];

  for (int i=0;i<n1;i++)
    {
      em[i]=electron[i]+positron[i]+gamma[i];
      ratio[i]=pi0[i]/em[i];
    }
  TCanvas *c = new TCanvas("c", "response", 200, 10, 1000, 800);
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("em and pi0 composition #pi^{-} showers ;Ekin [GeV];[%]");
  c->SetGrid();
  TGraph *gr_em = new TGraph(n1,Ein,em);
  gr_em->SetLineColor(4);
  gr_em->SetLineWidth(1);
  gr_em->SetMarkerColor(4);
  gr_em->SetMarkerStyle(20);
  gr_em->SetMarkerSize(1.2);
  gr_em->SetTitle("e^{-},e^{+},#gamma");
  mg->Add( gr_em);
  
  TGraph *gr_pi0 = new TGraph(n1,Ein,pi0);
  gr_pi0->SetLineColor(2);
  gr_pi0->SetLineWidth(1);
  gr_pi0->SetMarkerColor(2);
  gr_pi0->SetMarkerStyle(20);
  gr_pi0->SetMarkerSize(1.2);
  gr_pi0->SetTitle("#pi^{0}");
  mg->Add( gr_pi0);
  
  mg->SetMinimum(0.);
  mg->SetMaximum(100.);
  mg->Draw("apl");
  TLegend *leg = c->BuildLegend(.6, .75, 0.85, .9);
  leg->Draw();
  
  TCanvas *c2 = new TCanvas("c2", "ratio", 200, 10, 1000, 800);
  c2->SetGrid();
  TGraph *gr_ratio = new TGraph(n1,Ein,ratio);
  gr_ratio->SetTitle("#pi^{0}/em  ratio in #pi^{-} showers ;Ekin [GeV]; ");
  gr_ratio->SetLineColor(4);
  gr_ratio->SetLineWidth(1);
  gr_ratio->SetMarkerColor(4);
  gr_ratio->SetMarkerStyle(20);
  gr_ratio->SetMarkerSize(1.2);
  gr_ratio->SetMinimum(0.);
  gr_ratio->SetMaximum(1.);
  gr_ratio->Draw("apl");
  TLegend *leg2 = c2->BuildLegend(.6, .75, 0.85, .85);
  leg2->Draw();
  
  TCanvas *c3 = new TCanvas("c3", "Nr. of Pi0", 200, 10, 1000, 800);
  c3->SetGrid();
  TGraph *gr_npi0 = new TGraph(n1,Ein,npi0);
  gr_npi0->SetTitle("Nr of #pi^{0} produced  in #pi^{-} showers ;Ekin [GeV]; Nr. of #pi^{0} ");
  gr_npi0->SetLineColor(4);
  gr_npi0->SetLineWidth(1);
  gr_npi0->SetMarkerColor(4);
  gr_npi0->SetMarkerStyle(20);
  gr_npi0->SetMarkerSize(1.2);
  gr_npi0->SetMinimum(0.);
  gr_npi0->SetMaximum(30);
  gr_npi0->Draw("apl");
  TLegend *leg3 = c3->BuildLegend(.4, .8, 0.65, .85);
  leg3->Draw();
  

  TCanvas *c4 = new TCanvas("c4", "Nr. of p", 200, 10, 1000, 800);
  c4->SetGrid();
  TGraph *gr_np = new TGraph(n1,Ein,np);
  gr_np->SetTitle("Nr of p produced  in #pi^{-} showers ;Ekin [GeV]; Nr. of p ");
  gr_np->SetLineColor(4);
  gr_np->SetLineWidth(1);
  gr_np->SetMarkerColor(4);
  gr_np->SetMarkerStyle(20);
  gr_np->SetMarkerSize(1.2);
  gr_np->SetMinimum(0.);
  gr_np->SetMaximum(350);
  gr_np->Draw("apl");
  TLegend *leg4 = c4->BuildLegend(.4, .8, 0.65, .85);
  leg4->Draw();

  TCanvas *c5 = new TCanvas("c5", "Nr. of neutrons", 200, 10, 1000, 800);
  c5->SetGrid();
  TGraph *gr_nn = new TGraph(n2,Ein,nn);
  gr_nn->SetTitle("Nr of neutrons produced  in #pi^{-} showers ;Ekin [GeV]; Nr. of neutrons ");
  gr_nn->SetLineColor(4);
  gr_nn->SetLineWidth(1);
  gr_nn->SetMarkerColor(4);
  gr_nn->SetMarkerStyle(20);
  gr_nn->SetMarkerSize(1.2);
  gr_nn->SetMinimum(0.);
  gr_nn->SetMaximum(1000.);
  gr_nn->Draw("apl");
  TLegend *leg5 = c5->BuildLegend(.4, .8, 0.65, .85);
  leg5->Draw();



  TCanvas *c6 = new TCanvas("c6", "Nr. of #pi^{+}", 200, 10, 1000, 800);
  c6->SetGrid();
  TGraph *gr_npiplus = new TGraph(n1,Ein,npiplus);
  gr_npiplus->SetTitle("Nr of #pi^{+} produced  in #pi^{-} showers ;Ekin [GeV]; Nr. of #pi^{+} ");
  gr_npiplus->SetLineColor(4);
  gr_npiplus->SetLineWidth(1);
  gr_npiplus->SetMarkerColor(4);
  gr_npiplus->SetMarkerStyle(20);
  gr_npiplus->SetMarkerSize(1.2);
  gr_npiplus->SetMinimum(0.);
  gr_npiplus->SetMaximum(30);
  gr_npiplus->Draw("apl");
  TLegend *leg6 = c6->BuildLegend(.4, .8, 0.65, .85);
  leg6->Draw();
  
  
}
