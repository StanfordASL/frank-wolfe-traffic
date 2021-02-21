#!/bin/bash
set -x #echo on
#for city in Anaheim Massachusetts SiouxFalls NewYork
#for city in Anaheim  GoldCoast Massachusetts NewYork SiouxFalls Sydney
#for city in Anaheim Chicago Massachusetts NewYork SiouxFalls

# small: Anaheim Braess Friedrichshain Massachusetts MitteCenter Mitte-Prenzlauerberg-Friedrichshain-Center PrenzlauerbergCenter Tiergarten SiouxFalls Pigou
#medium: Berlin-Center Chicago NewYork

locations=../Differential_Pricing/Locations/_small
results=../Differential_Pricing/Results
exe=Build/Debug/Launchers/./AssignTraffic

#for city in Anaheim Braess Friedrichshain Massachusetts MitteCenter Mitte-Prenzlauerberg-Friedrichshain-Center PrenzlauerbergCenter Tiergarten SiouxFalls Pigou
for city in SiouxFalls
do
	alpha=100
	echo $exe -n 100 -i $locations/$city/edges.csv -od $locations/$city/od.csv -o $results/$city-constrained-$alpha -obj sys-opt -const_param $alpha -a constrained
done

#for i in */; do zip -r "${i%/}.zip" "$i"; done
