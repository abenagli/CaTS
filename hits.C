// Include files
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TKey.h"
//*****************************************************************************
// To run this macro in cint do:
// .include  $G4INCLUDE
// gSystem->Load("libCintex.so");
// gSystem->Load("$G4WORKDIR/tmp/$G4SYSTEM/CATS/libClassesDict.so");
// ROOT::Cintex::Cintex::Enable();
// .L hits.C++
// hits();
//*****************************************************************************
#include "include/TrackerHit.hh"
#include "include/CalorimeterHit.hh"
#include "include/DRCalorimeterHit.hh"
void hits()
{
  TFile fo("hits.root");
   
  std::vector<DRCalorimeterHit*>* hits;
  fo.GetListOfKeys()->Print();
 
  TIter next(fo.GetListOfKeys());
  TKey *key;
  double tot_en;
  while ((key=(TKey*)next()))
  {
    fo.GetObject(key->GetName(), hits);
 
    tot_en = 0;
    cout << "Collection: " << key->GetName() << endl;
    cout << "Number of hits: " << hits->size() << endl;
    for (int i=0;i!=hits->size();i++)
    {
      (*hits)[i]->Print();
    }         
  }
}
