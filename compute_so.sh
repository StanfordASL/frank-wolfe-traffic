#!/bin/bash
#set -x #echo on
#for city in Anaheim Massachusetts SiouxFalls NewYork
#for city in Anaheim  GoldCoast Massachusetts NewYork SiouxFalls Sydney
#for city in Anaheim Chicago Massachusetts NewYork SiouxFalls

# small: Anaheim Braess Friedrichshain Massachusetts MitteCenter Mitte-Prenzlauerberg-Friedrichshain-Center PrenzlauerbergCenter Tiergarten SiouxFalls Pigou
#medium: Berlin-Center Chicago NewYork

locations=../Differential_Pricing/Locations/
results=../Differential_Pricing/Results/test/
exe=Build/Release/Launchers/./AssignTraffic

#
#for city in Chicago NewYork
for city in SiouxFalls
do
	echo "running $city"
	for (( i=0; i <= 100; i++ ))
	do
		alpha=$(echo "0 + ( 0.01 * $i )" | bc)
		$exe -n 100 -i $locations/$city/edges.csv -od $locations/$city/od.csv -o $results/$city/2
	done 
done

#for i in */; do zip -r "${i%/}.zip" "$i"; done
