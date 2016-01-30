
#!/bin/bash

run=$1

group=$2

file=$3

echo ${run} ", " ${group} ", " ${file}

rm MuCSMerger_${run}_${group}.fcl
 
cp /uboone/app/users/kalousis/larsoft/srcs/uboonecode/uboone/MuCS/MuCSMerger.fcl MuCSMerger_${run}_${group}.fcl

sed -i 's/myrun/'${run}'/g' MuCSMerger_${run}_${group}.fcl

sed -i 's/mygroup/'${group}'/g' MuCSMerger_${run}_${group}.fcl

lar -c MuCSMerger_${run}_${group}.fcl -s $3

rm MuCSMerger_${run}_${group}.fcl
