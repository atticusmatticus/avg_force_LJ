#!/bin/bash

for ((i=1; i<=50; i++));
do
	if [ $i -lt 10 ]; then
		it=00$i
	else
		it=0$i
	fi
	sed -e "s/XX/$it/g" < plot.gpi > plot.$it.gpi
	gnuplot plot.$it.gpi
	rm plot.$it.gpi
done
