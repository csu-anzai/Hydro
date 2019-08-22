#/bin/bash

ms38='/data1/SNB_2007_Data/w7_ms38_l200-20070509/EnergyFnTime'
sg='/data1/SNB_2007_Data/w7_sg_l200-20070509/EnergyFnTime'
sy319='/data1/SNB_2007_Data/w7_sy319_l200-20070509/EnergyFnTime'

ms38_out="$ms38/out.txt"
sg_out="$sg/out.txt"
sy319_out="$sy319/out.txt"

rm $ms38_out
rm $sg_out
rm $sy319_out

for m in "0" "10" "20" "30" "40" "50" "60" "70" "80" "90" "100" "110" "120" "130" "140" "150"; do
	touch $ms38_out
	energy=`./eintDensSNBinary 2 $m | grep Total | cut -d: -f2`
	currentTime=`./eintDensSNBinary 2 $m | grep Time | cut -d: -f2`
	echo $m $currentTime $energy >> $ms38_out
	echo "MS38 -- Current results: " $m $currentTime $energy
done

for m in "0" "10" "20" "30" "40" "50" "60" "70" "80" "90" "100" "110" "120" "130" "140" "150"; do
	touch $sg_out
	energy=`./eintDensSNBinary 8 $m | grep Total | cut -d: -f2`
	currentTime=`./eintDensSNBinary 8 $m | grep Time | cut -d: -f2`
	echo $m $currentTime $energy >> $sg_out
	echo "SG -- Current results: " $m $currentTime $energy
done

for m in "0" "10" "20" "30" "40" "50" "60" "70" "80" "90" "100" "110" "120" "130" "140" "150"; do
	touch $sy319_out
	energy=`./eintDensSNBinary 9 $m | grep Total | cut -d: -f2`
	currentTime=`./eintDensSNBinary 9 $m | grep Time | cut -d: -f2`
	echo $m $currentTime $energy >> $sy319_out
	echo "SY319 -- Current results: " $m $currentTime $energy
done
