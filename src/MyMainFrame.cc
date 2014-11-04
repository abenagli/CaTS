#include "MyMainFrame.hh"



GroupBox::GroupBox(const TGWindow* p, const char* name, const char* title) :
  TGGroupFrame(p,name)
{
  // Group frame containing combobox and number entry
  
  TGHorizontalFrame* horz = new TGHorizontalFrame(this);
  AddFrame(horz, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));
  TGLabel* label = new TGLabel(horz,title);
  horz->AddFrame(label, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

  fCombo = new TGComboBox(horz);
  horz->AddFrame(fCombo, new TGLayoutHints(kLHintsRight | kLHintsExpandY, 5, 0, 5, 5));
  fCombo->Resize(100, 20);

  fEntry = new TGNumberEntry(this);
  AddFrame(fEntry, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));
}



void MyMainFrame::DoSetEventId(char* val)
{
  eventId = atoi(val);
}


void MyMainFrame::DoSetTimeType(Int_t id)
{
  if( id ==  8 ) timeType = 0;
  if( id == 32 ) timeType = 1;
  if( id == 16 ) timeType = 2;
}


void MyMainFrame::DoSetTimeSliceType(Int_t id)
{
  if( id ==  8 ) timeSliceType = 0;
  if( id == 32 ) timeSliceType = 1;
  if( id == 16 ) timeSliceType = 2;
}


void MyMainFrame::DoSetMinTimeSlice(char* val)
{
  minTimeSlice = atoi(val);
}


void MyMainFrame::DoSetMaxTimeSlice(char* val)
{
  maxTimeSlice = atoi(val);
}


void MyMainFrame::DoDraw()
{
  std::cout << "Slot DoDraw()" << std::endl;
  
  std::cout << "eventId: " << eventId << std::endl;
  std::cout << "timeSliceType: " << timeSliceType << std::endl;
  std::cout << "minTimeSlice: " << minTimeSlice << std::endl;
  std::cout << "maxTimeSlice: " << maxTimeSlice << std::endl;
  
  // Get a particular event
  Tevt -> GetEntry(eventId);
  std::cout << ">>> reading event " << eventId << " / " << Tevt->GetEntriesFast() << std::endl;
  
  std::string timeSliceString = "";
  
  if( timeType == 0 ) timeSliceString =+ "globalTime";
  if( timeType == 1 ) timeSliceString =+ "localTime1";
  if( timeType == 2 ) timeSliceString =+ "localTime2";
  
  if( timeSliceType == 0 ) timeSliceString += "Low";
  if( timeSliceType == 1 ) timeSliceString += "Med";
  if( timeSliceType == 2 ) timeSliceString += "Hig";
  std::cout << "ts: " << timeSliceString << "   " << timeSliceType << std::endl;
  
  std::map<G4String,std::map<G4String,std::map<G4ThreeVector,std::vector<G4VHit*> > > >* HCMap = event->GetHCMap();
  std::map<G4ThreeVector,std::vector<G4VHit*> > posHitMap = (*HCMap)["AbsorberVol_DRTSCalorimeter"][timeSliceString];
  
  
  // Draw something in the canvas
  if( c1 == NULL )
  {
    c1 = fCan->GetCanvas();
    c1->SetFillColor(0);
    c1->SetGrid();
    c1 -> Divide(2,1);
  }
  
  
  // Fill histograms
  if( h2_yx != NULL ) delete h2_yx;
  h2_yx = new TH2F(Form("h2_yx"),"y-x view",20,-500.,500.,20.,-500.,500.);
  
  if( h2_yz != NULL ) delete h2_yz;  
  h2_yz = new TH2F(Form("h2_yz"),"y-z view",200,0.,3000.,20.,-500.,500.);

  for(std::map<G4ThreeVector,std::vector<G4VHit*> >::const_iterator mapIt = posHitMap.begin(); mapIt != posHitMap.end(); ++mapIt)
  {
    double x = (mapIt->first).x();
    double y = (mapIt->first).y();
    double z = (mapIt->first).z();
    std::vector<G4VHit*> hitVec = mapIt->second;
    for(unsigned int vecIt = 0; vecIt < hitVec.size(); ++vecIt)
    {
      DRTSCalorimeterHit2* aHit = dynamic_cast<DRTSCalorimeterHit2*>(hitVec.at(vecIt));
      G4int timeSlice = aHit -> GetTimeSlice();
      if( timeSlice >= minTimeSlice && timeSlice <= maxTimeSlice )
      {
        h2_yx -> Fill(x,y,aHit->GetEdep());
        h2_yz -> Fill(z,y,aHit->GetEdep());
      }
    }
  }
  
  
  // Draw histograms
  c1 -> cd(1);
  gPad -> SetLogz();
  h2_yx -> Draw("COLZ");
  
  c1 -> cd(2);
  gPad -> SetLogz();
  h2_yz -> Draw("COLZ");
  
  
  // TCanvas::Update() draws the frame, after which it can be changed
  c1->Update();
  c1->GetFrame()->SetFillColor(kWhite);
  c1->GetFrame()->SetBorderSize(12);
  // c1->Modified();
  // c1->Update();
}

void MyMainFrame::DoExit()
{
   printf("Exit application...\n");
   gApplication->Terminate(0);
}

void MyMainFrame::SetStatusText(const char *txt, Int_t pi)
{
   // Set text in status bar.
   fStatusBar->SetText(txt,pi);
}

void MyMainFrame::EventInfo(Int_t ev, Int_t px, Int_t py, TObject *selected)
{
//  Writes the event status in the status bar parts

   const char *text0, *text1, *text3;
   char text2[50];
   text0 = selected->GetTitle();
   SetStatusText(text0,0);
   text1 = selected->GetName();
   SetStatusText(text1,1);
   if (ev == kKeyPress)
      sprintf(text2, "%c", (char) px);
   else
      sprintf(text2, "%d,%d", px, py);
   SetStatusText(text2,2);
   text3 = selected->GetObjectInfo(px,py);
   SetStatusText(text3,3);
}

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h, TFile* f) :
  TGMainFrame(p, w, h),
  eventId(0),
  timeType(0),
  timeSliceType(0),
  minTimeSlice(0),
  maxTimeSlice(999),
  c1(NULL),
  inFile(f),
  event(new Event()),
  runHeader(new RunHeader()),
  h2_yx(NULL),
  h2_yz(NULL)
{
  // initialize tree
  Tevt = (TTree*)(inFile->Get("EventTree"));
  Tevt -> SetBranchAddress("Event",&event);
  Trh = (TTree*)(inFile->Get("RunTree"));
  Trh -> SetBranchAddress("RunHeader",&runHeader);
  Trh -> GetEntry(0);
  
  
  // initialize variables
  std::vector<G4String> volumeList = runHeader -> GetVolumes();
  volumeList.push_back("AllVol");
  std::vector<G4String> particleList = runHeader -> GetParticleList();
  
  std::map<G4String,G4float> timeSliceSizes;
  std::map<G4String,G4float> minTimes;
  std::map<G4String,G4float> maxTimes;
  std::map<G4String,G4int> nTimeSlices;
  timeSliceSizes["Low"] = runHeader->GetTimeSliceSizeLow();
  timeSliceSizes["Med"] = runHeader->GetTimeSliceSizeMed();
  timeSliceSizes["Hig"] = runHeader->GetTimeSliceSizeHig();
  minTimes["Low"] = runHeader->GetMinTimeLow();
  minTimes["Med"] = runHeader->GetMinTimeMed();
  minTimes["Hig"] = runHeader->GetMinTimeHig();
  maxTimes["Low"] = runHeader->GetMaxTimeLow();
  maxTimes["Med"] = runHeader->GetMaxTimeMed();
  maxTimes["Hig"] = runHeader->GetMaxTimeHig();
  
  
  
  // create the header command bar
  TGHorizontalFrame* fHor = new TGHorizontalFrame(this,1000,600,kFixedWidth);
  
  TGGroupFrame* fGroup = new TGGroupFrame(fHor,"Event");
  fHor->AddFrame(fGroup, new TGLayoutHints(kLHintsLeft | kLHintsBottom, 5, 0, 5, 5));
  TGLabel* label = new TGLabel(fGroup,"Evt. number: ");
  fGroup->AddFrame(label, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 0, 5, 5));
  fNumber = new TGNumberEntry(fGroup, 0, 9,999, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax, 0, 999999);
  (fNumber->GetNumberEntry())->Connect("TextChanged(char*)", "MyMainFrame", this, "DoSetEventId(char*)");
  fGroup -> AddFrame(fNumber, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 0, 5, 5));
  
  fGroup1 = new TGButtonGroup(fHor,"Time type:");
  new TGRadioButton(fGroup1, "global", kTextTop);
  new TGRadioButton(fGroup1, "local1", kTextCenterY);
  new TGRadioButton(fGroup1, "local2", kTextBottom);
  fGroup1->SetButton(kTextTop);
  fGroup1->Connect("Pressed(Int_t)", "MyMainFrame", this, "DoSetTimeType(Int_t)");
  fHor->AddFrame(fGroup1, new TGLayoutHints(kLHintsLeft | kLHintsBottom, 5, 0, 5, 5));
  
  fGroup2 = new TGButtonGroup(fHor,"Time slice type:");
  new TGRadioButton(fGroup2, Form("low ([%.3f-%.3f] ns, %d slices, slice size: %.3f ns)",minTimes["Low"],maxTimes["Low"],int((maxTimes["Low"]-minTimes["Low"])/timeSliceSizes["Low"]),timeSliceSizes["Low"]), kTextTop);
  new TGRadioButton(fGroup2, Form("med ([%.3f-%.3f] ns, %d slices, slice size: %.3f ns)",minTimes["Med"],maxTimes["Med"],int((maxTimes["Med"]-minTimes["Med"])/timeSliceSizes["Med"]),timeSliceSizes["Med"]), kTextCenterY);
  new TGRadioButton(fGroup2, Form("hig ([%.3f-%.3f] ns, %d slices, slice size: %.3f ns)",minTimes["Hig"],maxTimes["Hig"],int((maxTimes["Hig"]-minTimes["Hig"])/timeSliceSizes["Hig"]),timeSliceSizes["Hig"]), kTextBottom);
  fGroup2->SetButton(kTextTop);
  fGroup2->Connect("Pressed(Int_t)", "MyMainFrame", this, "DoSetTimeSliceType(Int_t)");
  fHor->AddFrame(fGroup2, new TGLayoutHints(kLHintsLeft | kLHintsBottom, 5, 0, 5, 5));
  
  TGGroupFrame* fGroup3 = new TGGroupFrame(fHor,"Min/max time:");
  fHor->AddFrame(fGroup3, new TGLayoutHints(kLHintsLeft | kLHintsBottom, 5, 0, 5, 5));
  TGLabel* label3_1 = new TGLabel(fGroup3,"min time slice (-1 = underflow):");
  fGroup3->AddFrame(label3_1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 0, 5, 5));
  fNumber3_1 = new TGNumberEntry(fGroup3, 0, 9,999, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax, -1, 999);
  (fNumber3_1->GetNumberEntry())->Connect("TextChanged(char*)", "MyMainFrame", this, "DoSetMinTimeSlice(char*)");
  fGroup3 -> AddFrame(fNumber3_1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 0, 5, 5));
  TGLabel* label3_2 = new TGLabel(fGroup3,"max time slice (nMax+1 = overflow):");
  fGroup3->AddFrame(label3_2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 0, 5, 5));
  fNumber3_2 = new TGNumberEntry(fGroup3, 999, 9,999, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax, -1, 999);
  (fNumber3_2->GetNumberEntry())->Connect("TextChanged(char*)", "MyMainFrame", this, "DoSetMaxTimeSlice(char*)");
  fGroup3 -> AddFrame(fNumber3_2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 0, 5, 5));
  
  AddFrame(fHor);
  
  
  // Create the embedded canvas
  fCan = new TRootEmbeddedCanvas(0,this,500,400);
  int wid = fCan->GetCanvasWindowId();
  TCanvas* myc = new TCanvas("MyCanvas",10,10,wid);
  fCan->AdoptCanvas(myc);
  myc->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","MyMainFrame",this,"EventInfo(Int_t,Int_t,Int_t,TObject*)");
  
  AddFrame(fCan, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX  | kLHintsExpandY,0,0,1,1));
  
  
  // status bar
  Int_t parts[] = {45, 15, 10, 30};
  fStatusBar = new TGStatusBar(this, 50, 10, kVerticalFrame);
  fStatusBar->SetParts(parts, 4);
  fStatusBar->Draw3DCorner(kFALSE);
  AddFrame(fStatusBar, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
  
  
  // Create a horizontal frame containing two buttons
  TGHorizontalFrame *hframe = new TGHorizontalFrame(this, 200, 40);
  
  TGTextButton *draw = new TGTextButton(hframe, "&Draw");
  draw->Connect("Clicked()", "MyMainFrame", this, "DoDraw()");
  hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
  TGTextButton *exit = new TGTextButton(hframe, "&Exit ");
  exit->Connect("Pressed()", "MyMainFrame", this, "DoExit()");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  AddFrame(hframe, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
  
  
  // Set a name to the main frame   
  SetWindowName("Shower x-z and y-z projections");
  MapSubwindows();
  
  // Initialize the layout algorithm via Resize()
  Resize(GetDefaultSize());
  
  // Map main frame
  MapWindow();
  
  
}


MyMainFrame::~MyMainFrame()
{
   // Clean up main frame...
   Cleanup();
   delete fCan;
}
