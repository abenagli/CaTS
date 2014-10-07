#!/bin/tcsh -f
#
############################################################################
#
#
# Author: Hans wenzel
# File:  run_CaTS.csh
# wenzel@fnal.gov
#
############################################################################
# $1 = seed
# $2 = gdml file
# $3 = name of Particle 
# $4 = Energy of Particle 
# $5 = Number of events
# $6 = x-position 
# $7 = y-position
# $8 = z-position
# $9 = x-direction 
# $10= y-direction
# $11= z-direction
# $12= GEANT Physics list
# $13= Material 
#--------------------------------------------------------------------------- 
if ($#argv < 13) then
    echo " script needs 13 input variables"
    exit
endif
echo start now
/bin/date
pwd
@ SEED         = $1 
set gdmlfile   = `echo $2  | sed s/\'//g`
set Particle   = `echo $3  | sed s/\'//g`
set Energy     = `echo $4  | sed s/\'//g`
set NRofEvents = `echo $5  | sed s/\'//g`
set xposition  = `echo $6  | sed s/\'//g`
set yposition  = `echo $7  | sed s/\'//g`
set zposition  = `echo $8  | sed s/\'//g`
set xdirection = `echo $9  | sed s/\'//g`
set ydirection = `echo $10 | sed s/\'//g`
set zdirection = `echo $11 | sed s/\'//g`
setenv PHYSLIST  `echo $12 | sed s/\'//g`
set Material   = `echo $13 | sed s/\'//g`
set FILENAME=${Material}_${PHYSLIST}_${Particle}_${Energy}GeV.root
set ANA_FILENAME=${Material}_${PHYSLIST}_${Particle}_${Energy}GeV_analysis.root
set MAC_FILENAME=${Material}_${PHYSLIST}_${Particle}_${Energy}GeV.mac
set LOG_FILENAME=${Material}_${PHYSLIST}_${Particle}_${Energy}GeV.log
echo $FILENAME
echo $PHYSLIST
echo start
/bin/date
echo Particle ${Particle}
/bin/cat >  ${MAC_FILENAME} << EOF+ 
/CaTS/Analysis/Filename data/${ANA_FILENAME}
/CaTS/RootIO/Filename  data/${FILENAME}
# Shoot along various directions
/tracking/verbose 0
/gun/particle ${Particle}
/gun/direction ${xdirection} ${ydirection} ${zdirection}
/gun/position  ${xposition}  ${yposition}  ${zposition}
/gun/energy ${Energy} GeV
/run/beamOn ${NRofEvents}
exit
EOF+

/bin/more  ${MAC_FILENAME}
rm data/${LOG_FILENAME}
./CaTS ../CaTS/gdml/${gdmlfile}  ${MAC_FILENAME} >& data/${LOG_FILENAME} &
/bin/date
