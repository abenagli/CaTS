#ifndef MyMainFrame_hh
#define MyMainFrame_hh

#include "Event.hh"
#include "RunHeader.hh"
#include "DRTSCalorimeterHit2.hh"

#include <cmath>

#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGStatusBar.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TGNumberEntry.h>
#include <TGButtonGroup.h>
#include <TGComboBox.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TGToolTip.h>
#include <TStyle.h>

class Event;
class RunHeader;



class MyMainFrame : public TGMainFrame
{
private:
  TRootEmbeddedCanvas  *fCan;
  TGStatusBar          *fStatusBar;
  Long_t eventId;
  int timeType;
  int timeSliceType;
  int minTimeSlice;
  int maxTimeSlice;
  
  std::vector<G4String>* particleList;
  std::map<G4String,bool> particleEnabled;
  
  std::vector<G4String>* processList;
  std::map<G4String,bool> processEnabled;
  
  TCanvas* c1;
  TFile* inFile;
  TTree* Tevt;
  TTree* Trh;
  Event* event;
  RunHeader* runHeader;
  TH2F* h2_yx;
  TH2F* h2_yz;
  int nBinsX;
  int nBinsY;
  int nBinsZ;
  double* xAxis;
  double* yAxis;
  double* zAxis;
  double minZaxis;
  double maxZaxis;
  double Ein;
  double Eth;
  
public:
  MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h, TFile* f);
  virtual ~MyMainFrame();
  
  TGNumberEntry* fNumber;
  TGButtonGroup* fGroup1;
  TGButtonGroup* fGroup2;
  TGNumberEntry* fNumber3_1;
  TGNumberEntry* fNumber3_2;
  TGButtonGroup* fGroup4;
  TGButtonGroup* fGroup5;
  TGNumberEntry* fNumber6_1;
  TGNumberEntry* fNumber6_2;
  TGNumberEntry* fNumber6_3;
  
  void DoExit();
  void DoDraw();
  void DoSetEventId(char*);
  void DoSetTimeType(Int_t id);
  void DoSetTimeSliceType(Int_t id);
  void DoSetMinTimeSlice(char*);
  void DoSetMaxTimeSlice(char*);
  void DoSetEth(char*);
  void DoSetZAxisMin(char*);
  void DoSetZAxisMax(char*);
  void DoSetParticleEnabled(Bool_t val);
  void DoSetProcessEnabled(Bool_t val);
  void SetStatusText(const char *txt, Int_t pi);
  void EventInfo(Int_t ev, Int_t px, Int_t py, TObject *selected);
  
  ClassDef(MyMainFrame,0); 
};

#endif
