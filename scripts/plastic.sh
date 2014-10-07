export ENABLEOPTICAL=1
export ENABLESCINTILLATION=1
for ((i=0; i<10; i++ ))
do
    for ((j=0; j<10; j++ ))
    do
	ypos=`echo "-25.+2.5+5.*$i" | bc -l`
	zpos=`echo "-25.+2.5+5.*$j" | bc -l`
	rm -f plastic_${ypos}_${zpos}.mac 
	cat > plastic_${ypos}_${zpos}.mac << +EOF
	#/CaTS/Analysis/Filename plastic_electron_1GeV_${ypos}_${zpos}_analysis.root
	/CaTS/RootIO/Filename plastic_electron_1GeV_${ypos}_${zpos}_hits.root
        # Shoot along various directions
	/tracking/verbose 0
	/gun/direction 1. 0. 0.
	/gun/position -2 ${ypos} ${zpos} mm
	/gun/particle e-
	/gun/energy 1 GeV
	/run/beamOn 10
	exit
	+EOF
	./CaTS ../../CaTS/gdml/plastictile.gdml  plastic_${ypos}_${zpos}.mac
    done
done
