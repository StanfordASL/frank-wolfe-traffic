#!/bin/bash
#set -x #echo on
#for city in Anaheim Massachusetts SiouxFalls NewYork
#for city in Anaheim  GoldCoast Massachusetts NewYork SiouxFalls Sydney
#for city in Anaheim Chicago Massachusetts NewYork SiouxFalls

# small: Anaheim Braess Friedrichshain Massachusetts MitteCenter Mitte-Prenzlauerberg-Friedrichshain-Center PrenzlauerbergCenter Tiergarten SiouxFalls Pigou
#medium: Berlin-Center Chicago NewYork

locations=../Differential_Pricing/Locations/Elastic_small
results=../Differential_Pricing/Results/Elastic_small
exe=Build/Debug/Launchers/./AssignTraffic

#
#for city in Chicago NewYork
$exe -n 100 -elastic -f modified_bpr -i $locations/edges.csv -od $locations/od.csv -o $results/

#for i in */; do zip -r "${i%/}.zip" "$i"; done
