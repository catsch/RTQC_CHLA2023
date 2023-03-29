#!/bin/sh

GDACDATA="/DM/GDAC/coriolis/"
WORKDATA="/home/schmechtig/TRAITEMENT_FLOTTEUR/FLOAT_RT/CHLA/DATA/"

for mission in `cat coriolis_CHLA.list`
do

	if [ ! -d ${WORKDATA}${mission} ]
	then
		cp -fr ${GDACDATA}${mission} ${WORKDATA}${mission}
	fi

	ls -1 ${WORKDATA}${mission}"/profiles/"B*.nc | grep -v D.nc > liste_all_B
 
	R $mission --vanilla < ./CHLA_RTQC_ADJ.R

	mv *.png PNG_CORIOLIS/"$mission/"

done
