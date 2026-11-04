#!/bin/bash

cd mutationWhox

gfortran -w param_mod.mod.f90 evol.f90 -o run.e
cp run.e ../
rm *.mod

cd ..

cd addhox

#gfortran -w addhoxIniEqual.f90 -o addhox.e


gfortran -w addhoxIniEqual_uplow.f90 -o addhox.e


cp addhox.e ../
rm *.mod

cd ..



cd toothmakerWhox
gfortran -O2 -w source_original.f90 -o runt.e
cp runt.e ../
rm *.mod

cd ..

cd fitnessWfit2

gfortran -O2 -w overlapping.mod.f90 occlusion.f90 -o fitness.e
cp fitness.e ../
rm *.mod

cd ..


