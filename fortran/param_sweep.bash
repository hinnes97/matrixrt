#!/bin/bash

for lw_fac in 0.1 1 10 100
do
    for sw_fac in 0.1 1 10 100
    do
	sed -i "s/^[[:blank:]]*output_file.*/  output_file = 5_band_${lw_fac}_${sw_fac}.nc/" input.nml
	sed -i "s/^[[:blank:]]*sw_fac.*/  sw_fac = ${sw_fac}/" input.nml
	sed -i "s/^[[:blank:]]*lw_fac.*/  lw_fac = ${lw_fac}/" input.nml
	./main.exe
    done
done

