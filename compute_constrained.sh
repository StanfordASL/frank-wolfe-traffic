#!/bin/bash
#set -x #echo on
#for city in Anaheim Massachusetts SiouxFalls NewYork
#for city in Anaheim  GoldCoast Massachusetts NewYork SiouxFalls Sydney
#for city in Anaheim Chicago Massachusetts NewYork SiouxFalls

# small: Anaheim Braess Friedrichshain Massachusetts MitteCenter Mitte-Prenzlauerberg-Friedrichshain-Center PrenzlauerbergCenter Tiergarten SiouxFalls Pigou
#medium: Berlin-Center Chicago NewYork

locations=../Differential_Pricing/Locations/_small
results=../Differential_Pricing/Results
exe=Build/Release/Launchers/./AssignTraffic

#for city in Anaheim Braess Friedrichshain Massachusetts MitteCenter Mitte-Prenzlauerberg-Friedrichshain-Center PrenzlauerbergCenter Tiergarten SiouxFalls Pigou
for city in Friedrichshain
do
	alpha=1000
	$exe -n 100 -i $locations/$city/edges.csv -od $locations/$city/od.csv -o $results/$city-constrained-$alpha -obj sys_opt -const_param $alpha -a constrained
	$exe -n 100 -i $locations/$city/edges.csv -od $locations/$city/od.csv -o $results/$city-so -obj sys_opt
done

#for i in */; do zip -r "${i%/}.zip" "$i"; done
