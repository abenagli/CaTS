void ceren (){
  const Int_t n1= 7;
  Double_t Ein[n1]= {1.0,2.0,5.0,10.0,20.0,50.0,100.0};  
  Double_t Einerr[n1]={0.001,0.001,0.001,0.001,0.001,0.001,0.001};
  // electrons:
  Double_t Evis_em[n1]={0.9922,1.984,4.962,9.927,19.86,49.66,99.33};
  Double_t Evisrms_em[n1]={0.01098,0.01368,0.02308,0.03224,0.05164,0.09118,0.1514};
  Double_t NCeren_em[n1]={7.06e4,1.413e5,3.534e5,7.067e5,1.413e6,3.535e6,7.067e6};
  Double_t NCerenrms_em[n1]={1673,2581,4266,6464,1.627e4,2.471e4,7.132e4};
  Double_t Evis_em_rel[n1];
  Double_t Evisrms_em_rel[n1];
  // pi+:
  Double_t Evis_pip[n1-1]={0.8906,1.691,4.213,8.657,17.71,45.09};
  Double_t Evisrms_pip[n1-1]={0.1577,0.2375,0.4131,0.6335,1.003,1.982};
  Double_t NCeren_pip[n1-1]={4.098e4,7.14e4,1.822e5,4.13e5,9.072e5,2.421e6};
  Double_t NCerenrms_pip[n1-1]={2.062e4,2.987e4,6.042e4,1.07e5,1.875e5,4.029e5};
  Double_t Evis_pip_rel[n1-1];
  Double_t Evisrms_pip_rel[n1-1];
  // pi-:
  Double_t Evis_pim[n1]={0.932,1.727,4.283,8.773,17.83,45.38};
  Double_t Evisrms_pim[n1]={0.1336,0.2061,0.3954,0.5758,0.9109,1.872};
  Double_t NCeren_pim[n1-1]={4.431e4,7.562e4,1.879e5,4.183e5,9.147e5,2.467e6};
  Double_t NCerenrms_pim[n1-1]={2.037e4,2.99e4,6.034e4,1.049e5,1.829e5,3.946e5};
  Double_t Evis_pim_rel[n1-1];
  Double_t Evisrms_pim_rel[n1-1];
  // protons:
  Double_t Evis_p[n1]={0.8269,1.631,4.079,8.391,17.22,44.12};
  Double_t Evisrms_p[n1]={0.1184,0.1839,0.3645,0.5237,0.7951,1.504};
  Double_t NCeren_p[n1-1]={2.998e4,6.265e4,1.506e5,3.407e5,7.732e5,2.162e6};
  Double_t NCerenrms_p[n1-1]={1.183e4,2.437e4,5.409e4,8.433e4,1.353e5,2.81e5};
  Double_t Evis_p_rel[n1-1];
  Double_t Evisrms_p_rel[n1-1];
  // neutrons:
  Double_t Evis_n[n1]={0.7458,1.539,3.94,8.179,16.95,43.78};
  Double_t Evisrms_n[n1]={0.1191,0.2019,0.4016,0.5794,0.8836,1.631};
  Double_t NCeren_n[n1-1]={2.089e4,4.746e4,1.289e5,3.151e5,7.419e5,2.131e6};
  Double_t NCerenrms_n[n1-1]={1.023e4,2.028e4,5.029e4,8.188e4,1.329e5,2.729e5};
  Double_t Evis_n_rel[n1-1];
  Double_t Evisrms_n_rel[n1-1];
  // antiprotons:
  Double_t Evis_ap[n1]={2.339,3.242,5.901,10.32,19.11,45.92};
  Double_t Evisrms_ap[n1]={0.2306,0.2622,0.383,0.5671,0.9332,1.774};
  Double_t NCeren_ap[n1-1]={1.183e5,1.741e5,3.167e5,5.442e5,9.842e5,2.382e6};
  Double_t NCerenrms_ap[n1-1]={2.727e4,3.56e4,5.752e4,9.281e4,1.514e5,3.208e5};
  Double_t Evis_ap_rel[n1-1];
  Double_t Evisrms_ap_rel[n1-1];
  //K+
  Double_t Evis_kp[n1]={1.131,1.973,4.536,8.944,17.81,45.04};
  Double_t Evisrms_kp[n1]={0.1854,0.2474,0.4105,0.6603,1.343,2.943};
  Double_t NCeren_kp[n1-1]={5.192e4,8.689e4,1.931e5,4.023e5,8.571e5,2.331e6};
  Double_t NCerenrms_kp[n1-1]={1.868e4,3.083e4,5.636e4,9.216e4,1.652e5,3.485e5};
  Double_t Evis_kp_rel[n1-1];
  Double_t Evisrms_kp_rel[n1-1];
  //K-
  Double_t Evis_km[n1]={1.16,2.004,4.55,9.028,17.98,45.21};
  Double_t Evisrms_km[n1]={0.1896,0.2391,0.4078,0.7254,1.168,2.288};
  Double_t NCeren_km[n1-1]={4.93e4,8.454e4,1.94e5,4.267e5,8.901e5,2.376e6};
  Double_t NCerenrms_km[n1-1]={1.987e4,3.111e4,6.349e4,1.029e5,1.712e5,3.614e5};
  Double_t Evis_km_rel[n1-1];
  Double_t Evisrms_km_rel[n1-1];
  // muons (minimum ionizing/don't shower)
  const Int_t n2= 1;
  Double_t Ein_mu[n2]={10.0};  
  Double_t Einerr_mu[n2]={0.001};
  Double_t Evis_mup[n2]={3.293};
  Double_t Evisrms_mup[n2]={0.5776};
  Double_t NCeren_mup[n2]={2.795e5};
  Double_t NCerenrms_mup[n2]={3.375e4};
  Double_t Evis_mum[n2]={3.291};
  Double_t Evisrms_mum[n2]={0.5819};
  Double_t NCeren_mum[n2]={2.803e5};
  Double_t NCerenrms_mum[n2]={3.628e4};
  //
  Double_t Evis_mup_rel[n2];
  Double_t Evisrms_mup_rel[n2];
  Double_t Evis_mum_rel[n2];
  Double_t Evisrms_mum_rel[n2];
  //
  for (unsigned int i=0;i<n1;i++)
    {
	Evis_em_rel[i]=NCeren_em[i]/Evis_em[i];
	Evisrms_em_rel[i]= NCerenrms_em[i]/Ein[i]; 
    }
  for (unsigned int i=0;i<n1-1;i++)
    {
	Evis_pip_rel[i]=NCeren_pip[i]/Evis_pip[i];
	Evisrms_pip_rel[i]= NCerenrms_pip[i]/Ein[i]; 
	Evis_pim_rel[i]=NCeren_pim[i]/Evis_pim[i];
	Evisrms_pim_rel[i]= NCerenrms_pim[i]/Ein[i]; 
	Evis_p_rel[i]=NCeren_p[i]/Evis_p[i];
	Evisrms_p_rel[i]= NCerenrms_p[i]/Ein[i];
	Evis_n_rel[i]=NCeren_n[i]/Evis_n[i];
	Evisrms_n_rel[i]= NCerenrms_n[i]/Ein[i];
	Evis_ap_rel[i]=NCeren_ap[i]/Evis_ap[i];
	Evisrms_ap_rel[i]= NCerenrms_ap[i]/Ein[i];
	Evis_kp_rel[i]=NCeren_kp[i]/Evis_kp[i];
	Evisrms_kp_rel[i]= NCerenrms_kp[i]/Ein[i];
	Evis_km_rel[i]=NCeren_km[i]/Evis_km[i];
	Evisrms_km_rel[i]= NCerenrms_km[i]/Ein[i];
    }
 for (unsigned int i=0;i<n2;i++)
    {
	Evis_mup_rel[i]=NCeren_mup[i]/Evis_mup[i];
	Evisrms_mup_rel[i]= NCerenrms_mup[i]/Ein_mu[i];
 	Evis_mum_rel[i]=NCeren_mum[i]/Evis_mum[i];
	Evisrms_mum_rel[i]= NCerenrms_mum[i]/Ein_mu[i];
    }
  TCanvas *c2 = new TCanvas("c2", "relative response", 200, 10, 1000, 800);
  TMultiGraph *mg = new TMultiGraph();
  mg->SetMinimum(0);
  mg->SetMaximum(100000);
  mg->SetTitle("NCeren/Evis;Ein [GeV];NCeren/Evis");
  c2->SetGrid();
  TGraphErrors *gr_rel_em = new TGraphErrors(n1-1,Ein,Evis_em_rel,Einerr,Evisrms_em_rel);
  gr_rel_em->SetLineColor(4);
  gr_rel_em->SetLineWidth(1);
  gr_rel_em->SetMarkerColor(4);
  gr_rel_em->SetMarkerStyle(20);
  gr_rel_em->GetYaxis()->SetLimits(0.,2.5);
  gr_rel_em->SetTitle("e^{-}");
  mg->Add( gr_rel_em);
  TGraphErrors *gr_rel_mup = new TGraphErrors(n2,Ein_mu,Evis_mup_rel,Einerr_mu,Evisrms_mup_rel);
  gr_rel_mup->SetLineColor(7);
  gr_rel_mup->SetLineWidth(1);
  gr_rel_mup->SetMarkerColor(7);
  gr_rel_mup->SetMarkerStyle(20);
  gr_rel_mup->SetTitle("#mu^{+} (Evis=3.293GeV)");
  mg->Add( gr_rel_mup);
  TGraphErrors *gr_rel_mum = new TGraphErrors(n2,Ein_mu,Evis_mum_rel,Einerr_mu,Evisrms_mum_rel);
  gr_rel_mum->SetLineColor(2);
  gr_rel_mum->SetLineWidth(1);
  gr_rel_mum->SetMarkerColor(2);
  gr_rel_mum->SetMarkerStyle(20);
  gr_rel_mum->SetTitle("#mu^{-} (Evis=3.291 GeV)");
  mg->Add( gr_rel_mum);
  TGraphErrors *gr_rel_pip = new TGraphErrors(n1-1,Ein,Evis_pip_rel,Einerr,Evisrms_pip_rel);
  gr_rel_pip->SetLineColor(6);
  gr_rel_pip->SetLineWidth(1);
  gr_rel_pip->SetMarkerColor(6);
  gr_rel_pip->SetMarkerStyle(21);
  gr_rel_pip->SetTitle("#pi^{+} ");
  mg->Add( gr_rel_pip);
  TGraphErrors *gr_rel_pim = new TGraphErrors(n1-1,Ein,Evis_pim_rel,Einerr,Evisrms_pim_rel);
  gr_rel_pim->SetLineColor(13);
  gr_rel_pim->SetLineWidth(1);
  gr_rel_pim->SetMarkerColor(13);
  gr_rel_pim->SetMarkerStyle(22);
  gr_rel_pim->SetTitle("#pi^{-} ");
  mg->Add( gr_rel_pim);
  TGraphErrors *gr_rel_p = new TGraphErrors(n1-1,Ein,Evis_p_rel,Einerr,Evisrms_p_rel);
  gr_rel_p->SetLineColor(7);
  gr_rel_p->SetLineWidth(1);
  gr_rel_p->SetMarkerColor(7);
  gr_rel_p->SetMarkerStyle(20);
  gr_rel_p->SetTitle("protons");
  mg->Add( gr_rel_p);
  TGraphErrors *gr_rel_n = new TGraphErrors(n1-1,Ein,Evis_n_rel,Einerr,Evisrms_n_rel);
  gr_rel_n->SetLineColor(8);
  gr_rel_n->SetLineWidth(1);
  gr_rel_n->SetMarkerColor(8);
  gr_rel_n->SetMarkerStyle(22);
  gr_rel_n->SetTitle("neutrons");
  mg->Add( gr_rel_n);
  TGraphErrors *gr_rel_ap = new TGraphErrors(n1-1,Ein,Evis_ap_rel,Einerr,Evisrms_ap_rel);
  gr_rel_ap->SetLineColor(9);
  gr_rel_ap->SetLineWidth(1);
  gr_rel_ap->SetMarkerColor(9);
  gr_rel_ap->SetMarkerStyle(23);
  gr_rel_ap->SetTitle("antiproton ");
  TGraphErrors *gr_rel_kp = new TGraphErrors(n1-1,Ein,Evis_kp_rel,Einerr,Evisrms_kp_rel);
  gr_rel_kp->SetLineColor(11);
  gr_rel_kp->SetLineWidth(1);
  gr_rel_kp->SetMarkerColor(11);
  gr_rel_kp->SetMarkerStyle(21);
  gr_rel_kp->SetTitle("#K^{+}");
  TGraphErrors *gr_rel_km = new TGraphErrors(n1-1,Ein,Evis_km_rel,Einerr,Evisrms_km_rel);
  gr_rel_km->SetLineColor(12);
  gr_rel_km->SetLineWidth(1);
  gr_rel_km->SetMarkerColor(12);
  gr_rel_km->SetMarkerStyle(22);
  gr_rel_km->SetTitle("#K^{-}");
  mg->Add( gr_rel_ap);
  mg->Draw("apl");
  TLegend *leg = c2->BuildLegend(.6, .15, 0.85, .4);
  leg->Draw();
  pave = new TPaveText(20,90000,45,95000);
  pave->SetTextAlign(12);
  pave->AddText("Material: PbF2, Physics List: FTFPBERT");
  pave->Draw();

}
