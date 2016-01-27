
#!/bin/bash

run=$1

group=$2

file=$3

echo ${run} ", " ${group} ", " ${file}

rm MuCSDT_${run}_${group}.fcl
 
cp /uboone/app/users/kalousis/larsoft/srcs/uboonecode/uboone/MuCS/MuCSDT.fcl MuCSDT_${run}_${group}.fcl

sed -i 's/myrun/'${run}'/g' MuCSDT_${run}_${group}.fcl

sed -i 's/mygroup/'${group}'/g' MuCSDT_${run}_${group}.fcl

lar -c MuCSDT_${run}_${group}.fcl -s $3

rm MuCSDT_${run}_${group}.fcl
