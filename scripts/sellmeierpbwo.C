#include "math.h"
#include <iostream>
#include "TGraph.h"
using namespace std;
const double c=299792458.;           // speed of light in m/sec
const double h=4.13566743E-15;       // Planck constant in eVsec
// sellmeier coefficient from CMS TN/95-184
const double ns0 = 1.5861;
const double ns1 = 1.1062;
const double lams1 =270.63;
double lambdatoe(double lambda)
{
  // input   wavelength in nm 
  // return  energy in eV
 double   E  = (h*c)/(lambda*1.e-9);
 return E;
}
double etolambda(double E)
{
  // input  energy in eV
  // return   wavelength in nm 
 double   lambda  = (h*c)/(E*1.e-9);
 return lambda;
}
void sellmeierpbwo()
{
  for (int i =850;i>319;i--)
    {
      double lambda = double(i);
      double nsquare = 1+ ns0*ns0 +(ns1*ns1)/(1.-(lams1/lambda)*(lams1/lambda));
      double nord =  sqrt(nsquare); 
      cout << lambdatoe(lambda)<<"*eV "<< nord<<endl;

    }
}
