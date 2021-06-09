#!/usr/bin/bash
start=${1?Error: no start file \# given}
stop=${2?Error: no stop file \# given}
bX=${3?Error: no beamX \# given}
source /opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=/sphenix/user/shulga/tpc2019_install
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
 
for (( f=$start; f<$stop; f++ ));
do 
    let Xstart=f*1000 Xend=(f+1)*1000;
    A=$( printf '%06d' $Xstart )
    B=$( printf '%06d' $Xend )
    fname="/sphenix/user/shulga/Work/IBF/macros/detectors/sPHENIX/Files/DST_G4Hits_sHijing_0-12fm_"$A"_"$B".root" ;
    foutputname="/sphenix/user/shulga/Work/IBF/readDigitalCurrents/Files/hist_G4Hits_sHijing_0-12fm_bX"$bX"_"$A"_"$B".root" ;
    echo $fname ;
    echo $foutputname ;
    root -l -b -q /sphenix/user/shulga/Work/IBF/readDigitalCurrents/macros/Fun4All_FillDCMap.C\(1000,$Xstart,$bX,\"$fname\",\"$foutputname\"\)
done

echo all done
