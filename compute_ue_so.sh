#!/bin/bash
#set -x #echo on
#for city in Anaheim Massachusetts SiouxFalls NewYork
#for city in Anaheim  GoldCoast Massachusetts NewYork SiouxFalls Sydney
#for city in Anaheim Chicago Massachusetts NewYork SiouxFalls

# small: Anaheim Braess Friedrichshain Massachusetts MitteCenter Mitte-Prenzlauerberg-Friedrichshain-Center PrenzlauerbergCenter Tiergarten SiouxFalls Pigou
#medium: Berlin-Center Chicago NewYork

locations=../Differential_Pricing/Locations/_medium
results=../Differential_Pricing/Results
exe=Build/Release/Launchers/./AssignTraffic

#for city in Anaheim Braess Friedrichshain Massachusetts MitteCenter Mitte-Prenzlauerberg-Friedrichshain-Center PrenzlauerbergCenter Tiergarten SiouxFalls Pigou
for city in Chicago NewYork
do
	echo "running $city"
	for (( i=0; i <= 20; i++ ))
	do
		alpha=$(echo "0 + ( 0.05 * $i )" | bc)
		$exe -n 200 -i $locations/$city/edges.csv -od $locations/$city/od.csv -o $results/$city-$alpha -obj combined_eq -ce_param $alpha
	done 
done

#for i in */; do zip -r "${i%/}.zip" "$i"; done
