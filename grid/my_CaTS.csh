#!/bin/csh
#
############################################################################
#
#
# Author: Hans wenzel
# File:  my_CaTS.csh
# wenzel@fnal.gov
#
############################################################################
# $1  = Nr of processes
# $2  = name of Particle
# $3  = Energy of Particle
# $4  = Number of events per process
# $5  = x-position 
# $6  = y-position
# $7  = z-position
# $8  = x-direction 
# $9  = y-direction
# $10 = z-direction
# $11 = GEANT Physics list
#---------------------------------------------------------------------------

if ($#argv < 11) then
    echo " scripts needs 11 input variables"
    exit
endif

echo " submitting: "  ${1} "processes to the grid"
echo start
/bin/date
rm -f  run_CaTS_grid
cat > run_CaTS_grid << +EOF
universe = grid
globusscheduler = fngp-osg.fnal.gov/jobmanager-condor
executable = $PWD/run_CaTS_grid.csh
transfer_executable = true
transfer_output = true
transfer_error = true
transfer_executable = true
environment = "ClusterProcess=\$(Cluster)-\$(Process)"
log = ./log/CaTS_grid.log.\$(Cluster).\$(Process)
notification = NEVER
output = ./stdout/CaTS_grid.out.\$(Cluster).\$(Process)
error = ./stderr/CaTS_grid.err.\$(Cluster).\$(Process)
stream_output = false
stream_error = false
ShouldTransferFiles = YES
WhenToTransferOutput = ON_EXIT
globusrsl = (maxwalltime=90000)(jobtype=single)
Arguments = \$(Cluster) \$(Process) '${2}' '${3}' '${4}' '${5}' '${6}' '${7}' '${8}' '${9}' '${10}' '${11}' 

queue ${1}
+EOF

condor_submit run_CaTS_grid
