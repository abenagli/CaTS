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
#ifndef INCLUDE_ROOTIO_HH 
#define INCLUDE_ROOTIO_HH 1

// ROOT Include files
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
// CaTS Include files
#include "products.hh"
#include "Event.hh"
#include "RunHeader.hh"

class RootIOMessenger;

class RootIO {
public:
    virtual ~RootIO();

    static RootIO* GetInstance();
    void SetFileName(G4String);
    void BeginOfRun();
    void Write(Event*);
    void Write(RunHeader*);
    void Close();

protected:
    RootIO();

private:
    G4String FileName;
    bool evtinitialized;
    bool runinitialized;
    bool treeinitialized;
    TFile* fo;
    TTree* ftree;
    TTree* fhtree;
    TBranch* frunbranch;
    TBranch* fevtbranch;
    Long64_t fnb;
    RootIOMessenger* pMessenger; // pointer to the Messenger

};
#endif // INCLUDE_ROOTIO_HH
