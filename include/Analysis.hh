/* ------------------------------------------------------------------------
            |\___/|       
            )     (    
           =\     /=
             )===(
            /     \         CaTS: Calorimeter and Tracker Simulation
            |     |         Author: Hans Wenzel (Fermilab)
           /       \
           \       /
            \__  _/
              ( (
               ) )
              (_(
-------------------------------------------------------------------------*/
///////////////////////////////////////////////////////////////////////////////
// File: Analysis.hh
// Description: Analysis is a singleton class and interfaces all user
//              analysis code
///////////////////////////////////////////////////////////////////////////////
#ifndef Analysis_h
#define Analysis_h 1
#include "globals.hh"
#include "TFile.h"
#include <map>
#include <set>

class AnalysisMessenger;

class Analysis {
public:
    virtual ~Analysis();

public:
    TFile* output;
    void BeginOfRun(G4int n);
//    void BeginOfEvent(G4int n);
    void EndOfRun(G4int n);
    void EndOfEvent();
    void SetFileName(G4String);
    static Analysis* getInstance();

    TDirectory *topdir;

//    TDirectory *NCerenbypartdir;
    TDirectory *Stackingdir; // Directory for StackingAction  
    TDirectory *SDdir; // Directory for sensitive detectors
    std::map<G4String, TDirectory*> SDDirectorymap;
    std::map<G4String, TDirectory*>::iterator SDDirectoryiter;
    G4int RunNr;

private:
    Analysis();
    Long64_t fnb;
    Int_t fevts;
    std::set<G4String> SDset;
    static Analysis* instance;
  
    G4String filename;
    G4double timeslice;
    AnalysisMessenger* pMessenger;
    static Analysis* analysis;
};
#endif
