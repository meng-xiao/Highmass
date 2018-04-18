# Highmass, preparation of the inputs for the final combination

## to build intermediate root files that contains the final variables, dbkg_kin, ZZMass, weight

root -q -b qqzz.c

root -q -b zx.c

root -q -b data.c

## to build the 2D pdf for the bkg and import to a workspace 

runbkg.sh

## to build the 2D pdf for the sig and import to a workspace

run.sh
