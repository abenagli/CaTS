// g++ -Wall -o drawShowerEvolution.exe `root-config --cflags --glibs` setTDRStyle.cc drawShowerEvolution.cpp

#include "setTDRStyle.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <map>

#include "TFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"

bool fexists(const std::string& fileName)
{
  ifstream ifile(fileName.c_str());
  return ifile;
}



int main(int argc, char** argv)
{
  if( argc < 5 ) return -1;
  
  int event = atoi(argv[1]);
  std::string histoType = std::string(argv[2]);
  std::string timeSlice = std::string(argv[3]);
  int drawPerParticle = atoi(argv[4]);
  
  setTDRStyle();
  
  std::string inFileName = "./drawShower.root";
  TFile* inFile = TFile::Open(inFileName.c_str(),"READ");
  
  std::vector<std::string> particleList;
  particleList.push_back("fragment");
  particleList.push_back("proton");
  particleList.push_back("neutron");
  particleList.push_back("e");
  particleList.push_back("gamma");
  particleList.push_back("mu");
  particleList.push_back("pi");
  particleList.push_back("meson");
  particleList.push_back("baryon");
  particleList.push_back("em");
  
  std::map<std::string,int> colorMap;
  colorMap["fragment"] = kBlack;
  colorMap["proton"] = kBlue;
  colorMap["neutron"] = kGray;
  colorMap["e"] = kRed;
  colorMap["gamma"] = kRed;
  colorMap["mu"] = kGreen;
  colorMap["pi"] = kOrange-5;
  colorMap["meson"] = kOrange-7;
  colorMap["baryon"] = kBlue+1;
  colorMap["em"] = kYellow;
  
  TH1F* h_timeSlice = (TH1F*)( inFile->Get(Form("h_timeSlice%s",timeSlice.c_str())) );
  
  TLatex* latex = NULL;
  int nSlices = h_timeSlice->GetNbinsX();
  
  
  
  if( !drawPerParticle )
  { 
    std::string fileName = Form("showerEvolution_%s_%s_evt%d.gif",histoType.c_str(),timeSlice.c_str(),event);
    if( fexists(fileName) ) remove(fileName.c_str());
    for(int slice = 1; slice <= nSlices; ++slice)
    {
      std::cout << "drawing slice " << slice << " / " << nSlices << "\r" << std::flush;
      
      TH2F* h = (TH2F*)( inFile->Get(Form("event%d/h2_%s_globalSlice%s%d_evt%d",event,histoType.c_str(),timeSlice.c_str(),slice,event)) );
      
      //h -> SetMarkerStyle(20);
      //h -> SetMarkerSize(0.5);
      
      TCanvas* c = new TCanvas("c","c");
      c -> cd();
      c -> SetLogz();
      if( histoType == "xy" )
      {
        h -> GetXaxis() -> SetTitle("x (mm)");
        h -> GetYaxis() -> SetTitle("y (mm)");
      }
      if( histoType == "zy" )
      {
        h -> GetXaxis() -> SetTitle("z (mm)");
        h -> GetYaxis() -> SetTitle("y (mm)");
      }
      h -> Draw("");
      
      if( latex != NULL ) delete latex;
      latex = new TLatex(0.13,0.96,Form("Time: %.2f ns",h_timeSlice->GetBinWidth(1)*(slice-1)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> Draw("same");
      
      c -> Print(Form("%s+",fileName.c_str()));
      delete c;
    }
  }
  
  
  
  if( drawPerParticle )
  {
    std::string fileName = Form("showerEvolution_%s_%s_evt%d_perParticle.gif",histoType.c_str(),timeSlice.c_str(),event);
    if( fexists(fileName) ) remove(fileName.c_str());
    
    TLegend* legend = new TLegend(0.70,0.65,0.90,0.90);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04);  
    
    for(int slice = 1; slice <= nSlices; ++slice)
    {
      std::cout << "drawing slice " << slice << " / " << nSlices << "\r" << std::flush;
      
      bool isFirst = true;
      
      TCanvas* c = new TCanvas("c","c");
      c -> cd();
      c -> SetLogz();
      
      for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
      {
        std::string particleName = particleList.at(partIt);
        
        TH2F* h = (TH2F*)( inFile->Get(Form("event%d/h2%s_%s_globalSlice%s%d_evt%d",event,particleName.c_str(),histoType.c_str(),timeSlice.c_str(),slice,event)) );
        
        h -> SetMarkerColor(colorMap[particleName]);
        h -> SetMarkerStyle(20);
        h -> SetMarkerSize(0.5);
        
        
        if( slice == 1 )
        {
          legend -> AddEntry(h,Form("%s",particleName.c_str()),"P");
        }
        
        
        c -> cd();
        c -> SetLogz();
        if( isFirst )
        {
          if( histoType == "xy" )
          {
            h -> GetXaxis() -> SetTitle("x (mm)");
            h -> GetYaxis() -> SetTitle("y (mm)");
          }
          if( histoType == "zy" )
          {
            h -> GetXaxis() -> SetTitle("z (mm)");
            h -> GetYaxis() -> SetTitle("y (mm)");
          }
          h -> Draw("");
          isFirst = false;
        }
        else h -> Draw("same");
        
      }
      
      if( latex != NULL ) delete latex;
      latex = new TLatex(0.13,0.96,Form("Time: %.2f ns",h_timeSlice->GetBinWidth(1)*(slice-1)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> Draw("same");
      legend -> Draw("same");
      c -> Print(Form("%s+",fileName.c_str()));
      delete c;
    }
  }
}
