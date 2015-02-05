#include "MyMainFrame.hh"
#include "setTDRStyle.hh"

#include "Cintex/Cintex.h"
#include "TApplication.h"

#include <iostream>



int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();
  setTDRStyle();
  
  if( argc < 2 )
  {
    std::cout << "indicate a root file. Exiting..." << std::endl;
    exit(-1);
  }
  TFile* inFile = TFile::Open(argv[1]);
  
  
  TApplication* theApp = new TApplication("App", &argc, argv);
  // Popup the GUI...
  new MyMainFrame(gClient->GetRoot(), 200, 200, inFile);
  theApp -> Run();
  
  return 0;
}
