#!/bin/tcsh -f                                                        
#                                                                                             
############################################################################                                 
#                                                                                                            
#
# Author: Hans wenzel
# File:  combine.csh
# wenzel@fnal.gov
#
############################################################################
# $1 = name of Material
# $2 = GEANT Physics list
# $3 = Name of Particle 
# $4 = Energy of Particle 
#--------------------------------------------------------------------------- 
if ($#argv < 4) then
    echo " script needs 4 input variables"
    exit
endif
source  /home/wenzel/root/bin/thisroot.csh
. /home/wenzel/geant4.9.6.p02-install/bin/geant4.csh
. /home/wenzel/geant4.9.6.p02-install/share/Geant4-9.6.2/geant4make/geant4make.csh
set Material   = `echo $1  | sed s/\'//g`
set PHYSLIST   = `echo $2 | sed s/\'//g` 
set PHYSLISTS  = `echo $2 | sed s/\'//g | sed s/_//g` 
set Particle   = `echo $3  | sed s/\'//g`
set Energy     = `echo $4  | sed s/\'//g`
set File = `ls /data/CaTS_newdata/${Material}_${PHYSLIST}_${Particle}_${Energy}*_0_hits.root`
#
echo $Material
echo $Particle
echo $Energy
echo $PHYSLISTS
echo $File
rm -f combine.C
/bin/cat > combine.C << EOF+
#include "TTree.h"
#include "TChain.h"
#include "TSystem.h"
#include "TInterpreter.h"

void combine()
{ 
gSystem->Load("libCintex.so");
 ROOT::Cintex::Cintex::Enable();
gSystem->Load("/home/wenzel/NetBeansProjects/CaTS-install/lib/libClassesDict.so");
TChain ch("Events");
EOF+
cat >> combine.C <<EOF+ 
`ls -1  /data/CaTS_newdata/${Material}_${PHYSLIST}_${Particle}_${Energy}*hits* | sed s'&/data/CaTS_newdata/&ch.Add("/data/CaTS_newdata/&g'| sed s'&.root&.root");&g'`
EOF+
cat >> combine.C <<EOF+ 
ch.Merge("/data/CaTS_combineddata/${Material}_${PHYSLISTS}_${Particle}_${Energy}_hits.root");
TFile* f = 0;
 TFile* newFile = 0;
 TTree* t = 0;
 TTree* newTree = 0;
 f = new TFile("${File}");
 newFile = new TFile("/data/CaTS_combineddata/${Material}_${PHYSLISTS}_${Particle}_${Energy}_hits.root", "update");

 t = (TTree*) f->Get("Runheader"); 
 newTree = t->CloneTree(0);
 TTreeCloner cloner(t, newTree, "", TTreeCloner::kNoWarnings);
 newTree->SetEntries(newTree->GetEntries() + t->GetEntries());
 cloner.Exec();
 newTree->Write();
 gROOT->ProcessLine(".q");
}
EOF+

#/bin/more combine.C
root combine.C 
/bin/date
