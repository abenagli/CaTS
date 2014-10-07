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
// File: Analysis.cc
// Description: Analysis interfaces all user analysis code
///////////////////////////////////////////////////////////////////////////////

#include "Analysis.hh"
#include "AnalysisMessenger.hh"
#include "StackingAction.hh"


Analysis* Analysis::instance = 0;

Analysis::Analysis() {
    pMessenger = new AnalysisMessenger(this);
}

Analysis::~Analysis() {
}

Analysis* Analysis::getInstance() {
    if (instance == 0) instance = new Analysis();
    return instance;
}

/* 
   This member reset the histograms and it is called at the begin
   of each run; 
 */
void Analysis::BeginOfRun(G4int Run) {
    RunNr = Run;
    std::stringstream convert; // stream used for the conversion
    convert << RunNr; // insert the textual representation of 'Number' in the characters in the stream
    G4String RunNrString = convert.str();
    //std::vector<std::string> y = split(filename, '.');
    //std::string Classname = y[1];
    //G4String FileName = y[0] + "_Run" + RunNrString + "." + y[1];
    G4String FileName = filename + RunNrString;
    G4cout << "Analysis::BeginOfRun:  " << RunNr
            << "   " << RunNrString
            << "  " << FileName << G4endl;
    output = new TFile(FileName.c_str(), "RECREATE", "Histogram ROOT file");
    topdir = output->CurrentDirectory();
    output->SetCompressionLevel(1);
    SDdir = output->mkdir("SensitiveDetectors");
    Stackingdir = output->mkdir("StackingAction");
}

//void Analysis::BeginOfEvent(G4int Evt) {
//}
//  This member is called at the end of each run

void Analysis::EndOfRun(G4int) {

    output->Write("", TObject::kOverwrite);
    output->Close();
}

// This member is called at the end of every event

void Analysis::EndOfEvent() {
    StackingAction::getInstance()->EndofEvent();
    return;
}

void Analysis::SetFileName(G4String fname) {
    filename = fname;
}


