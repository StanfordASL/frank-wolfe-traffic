#!/bin/bash
#set -x #echo on
#for city in Anaheim Massachusetts SiouxFalls NewYork
#for city in Anaheim  GoldCoast Massachusetts NewYork SiouxFalls Sydney
#for city in Anaheim Chicago Massachusetts NewYork SiouxFalls

# small: Anaheim Braess Friedrichshain Massachusetts MitteCenter Mitte-Prenzlauerberg-Friedrichshain-Center PrenzlauerbergCenter Tiergarten SiouxFalls Pigou
#medium: Berlin-Center Chicago NewYork

locations=../Differential_Pricing/Locations
results=../Differential_Pricing/Results
exe=Build/Release/Launchers/./AssignTraffic

#for city in Anaheim Braess Friedrichshain Massachusetts MitteCenter Mitte-Prenzlauerberg-Friedrichshain-Center PrenzlauerbergCenter Tiergarten SiouxFalls  Pigou
#for city in Anaheim Friedrichshain Massachusetts Mitte-Prenzlauerberg-Friedrichshain-Center PrenzlauerbergCenter SiouxFalls Tiergarten
for city in PrenzlauerbergCenter SiouxFalls Tiergarten  
do
	echo "running $city"
	for (( i=0; i <= 20; i++ ))
	do
		alpha=$(echo "1 + ( 0.05 * $i )" | bc)
		echo "with $alpha"
		$exe -n 100 -i $locations/$city/edges.csv -od $locations/$city/od.csv -o $results/$city-constrained-$alpha -obj sys_opt -const_param $alpha -a constrained
	done	
done

# zip all subfolders of a given folder
# for i in */; do zip -r "${i%/}.zip" "$i"; done
