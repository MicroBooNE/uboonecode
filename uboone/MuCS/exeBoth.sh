
#!/bin/bash

run=$1

group=$2

file=$3

echo ${run} ", " ${group} ", " ${file}

source srcs/uboonecode/uboone/MuCS/exeMuCSDT.sh ${run} ${group} ${file}

source srcs/uboonecode/uboone/MuCS/exeMuCSMerger.sh ${run} ${group} ${file}
