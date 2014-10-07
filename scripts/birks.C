#include"TF1.h"
#include"TLegend.h"
#include"TAxis.h"
// dedx in [MeV*cm^2/g]
// values for the scintillator Birks' coefficient C1, with C2 = 0.
//C1_scintillator = 1.31e-2; C2_scintillator = 0.0;  // Standard Birks law.
//C1_scintillator = 8.35e-3; C2_scintillator = 0.0;  // Zeus SCSN38, lower limit (-35%) IEEE TNS Vol. 39 NO.4 (1992), 511
//C1_scintillator = 1.59e-2; C2_scintillator = 0.0;  // Pilot B (same paper), upper limit (+21%)
// BGO: kB = 6.5 mum/MeV                             // NIM A439 (2000) 158-166
// with rho =7.13 g/cm^3--> C1=4.6345e-2 C2=0.0
// CsI(Tl): kB = 1.52 mum/MeV                             // NIM A439 (2000) 158-166
// with rho = 4.51 g/cm^3--> C1=6.8552e-4 C2=0.0
// GSO(Ce): kB = 5.25 mum/MeV                             // NIM A439 (2000) 158-166
// with rho = 6.7 g/cm^3--> C1=3.5175e-3 C2=0.0
double birksf(double * x,double *p) 
{      
  double dedx = x[0];
  return 1.0 / (1.0 + p[0] * dedx + p[1] * dedx * dedx);
}
void birks()
{
 TF1 *bf4  = new TF1("birks4",birksf,0,100,2);
 bf4->SetParameters(4.6345e-2,0.0);
 bf4->SetLineWidth(2); 
 bf4->SetLineColor(7);
 bf4->GetXaxis()->SetTitle("-dE/dx  [MeV/g cm^{2}]");
 bf4->GetYaxis()->SetTitle("Birks Factor");
 bf4->Draw();
 TF1 *bf  = new TF1("birks",birksf,0,100,2);
 bf->SetParameters(1.29e-2,9.59e-6);
 bf->SetLineWidth(2); 
 bf->SetLineColor(3);
 bf->Draw("SAME");
 TF1 *bf1  = new TF1("birks1",birksf,0,100,2);
 bf1->SetParameters(8.35e-3,0.0);
 bf1->SetLineWidth(2); 
 bf1->SetLineColor(4);
 bf1->Draw("SAME");
 TF1 *bf2  = new TF1("birks2",birksf,0,100,2);
 bf2->SetParameters(1.31e-2,0.0);
 bf2->SetLineWidth(2); 
 bf2->SetLineColor(2);
 bf2->Draw("SAME");
 TF1 *bf3  = new TF1("birks3",birksf,0,100,2);
 bf3->SetParameters(1.59e-2,0.0);
 bf3->SetLineWidth(2); 
 bf3->SetLineColor(6);
 bf3->Draw("SAME");
 TF1 *bf5  = new TF1("birks5",birksf,0,100,2);
 bf5->SetParameters(6.8552e-4,0.0);
 bf5->SetLineWidth(2); 
 bf5->SetLineColor(21);
 bf5->Draw("SAME");
 TF1 *bf6  = new TF1("birks6",birksf,0,100,2);
 bf6->SetParameters(3.5175e-3,0.0);
 bf6->SetLineWidth(2); 
 bf6->SetLineColor(8);
 bf6->Draw("SAME");
 TLegend *legend = new TLegend(.44, .64, .99, .99);
 legend->AddEntry(bf,  "CMS/ATLAS");
 legend->AddEntry(bf2, "Standard Birks law.");
 legend->AddEntry(bf1, "Zeus SCSN38, lower limit (-35%) IEEE TNS Vol. 39 NO.4 (1992), 511");
 legend->AddEntry(bf3, "Zeus SCSN38, upper limit (+21%) IEEE TNS Vol. 39 NO.4 (1992), 511");
 legend->AddEntry(bf4, "BGO, NIM A439 (2000) 158-166");
 legend->AddEntry(bf5, "CsI(Tl), NIM A439 (2000) 158-166");
 legend->AddEntry(bf6, "GSO(Ge), NIM A439 (2000) 158-166");
 legend->Draw();
}
