#!/bin/bash

declare -a STEPS=( ${@:6:$(($#-5))} ) # -g option to declare this array global
                                      # is not present in bash version 4.1.2(1)-release
                                      # so the array has been declared in the global scope

function initialize
{
    trap 'LASTERR=$?; echo -e "\nCI MSG BEGIN\n `basename $0`: error at line ${LINENO}\n Stage: ${STEPS[STEP]}\n Task: ${TASKSTRING}\n exit status: ${LASTERR}\nCI MSG END\n"; exit ${LASTERR}' ERR

    echo "initialize $@"

    NEVENTS=$1
    STEP=$2
    LARSOFT_REFERENCE_VERSION=$3
    BASEFILENAME=$4
    EXPCODE=$5

    #declare -ga STEPS=( ${@:6:$(($#-5))} )

    declare -a SeedsFiles=(GenRandomSeeds_Ref.dat G4RandomSeeds_Ref.dat DetSimRandomSeeds_Ref.dat Reco1RandomSeeds_Ref.dat)
    if [ x${SeedsFiles[STEP-1]} != x ]; then
        xrdcp --nopbar xroot://fndca1.fnal.gov/pnfs/fnal.gov/usr/uboone/scratch/users/vito/ci_tests_inputfiles/${SeedsFiles[STEP-1]} .
    fi

    INPUT_FILE="${BASEFILENAME}_Reference_${STEPS[STEP-1]}_${LARSOFT_REFERENCE_VERSION}.root"
    if [ x"${STEPS[STEP-1]}" == xnone ]; then INPUT_FILE=""; fi

    if [ x$INPUT_FILE != x ]; then
        xrdcp --nopbar xroot://fndca1.fnal.gov/pnfs/fnal.gov/usr/uboone/scratch/users/vito/ci_tests_inputfiles/${INPUT_FILE} .
    fi

    REFERENCE_FILE="${BASEFILENAME}_Reference_${STEPS[STEP]}_${LARSOFT_REFERENCE_VERSION}.root"
    xrdcp --nopbar xroot://fndca1.fnal.gov/pnfs/fnal.gov/usr/uboone/scratch/users/vito/ci_tests_inputfiles/${REFERENCE_FILE} .

    OUTPUT_FILE="${BASEFILENAME}_Current_${STEPS[STEP]}.root"

    FHiCL_FILE="ci_test_${STEPS[STEP]}_${EXPCODE}.fcl"

    echo "Input file:  $INPUT_FILE"
    echo "Output file: $OUTPUT_FILE"
    echo "FHiCL file:  $FHiCL_FILE"
    echo
    echo -e "\nRunning\n `basename $0` $@"


}


function larsoft_data_production
{

    trap 'LASTERR=$?; echo -e "\nCI MSG BEGIN\n `basename $0`: error at line ${LINENO}\n Stage: ${STEPS[STEP]}\n Task: ${TASKSTRING}\n exit status: ${LASTERR}\nCI MSG END\n"; exit ${LASTERR}' ERR

    echo -e "\nNumber of events for ${STEPS[STEP]} step: $NEVENTS\n"
    echo lar --rethrow-default -n ${NEVENTS} -o ${OUTPUT_FILE} --config ${FHiCL_FILE} ${INPUT_FILE}
    echo

    lar --rethrow-default -n $NEVENTS -o $OUTPUT_FILE --config $FHiCL_FILE ${INPUT_FILE}

}

function compare_data_products
{
    trap 'LASTERR=$?; echo -e "\nCI MSG BEGIN\n `basename $0`: error at line ${LINENO}\n Stage: ${STEPS[STEP]}\n Task: ${TASKSTRING}\n exit status: ${LASTERR}\nCI MSG END\n"; exit ${LASTERR}' ERR

    declare -a CHECKMSG=("Check for added/removed data products" "Check for differences in the size of data products")


    if [ ${COMPAREINIT} -eq 0 ]; then

        #echo lar --rethrow-default -n $NEVENTS --config eventdump.fcl ${BASEFILENAME}_Reference_${STEPS[STEP]}_${LARSOFT_REFERENCE_VERSION}.root

        lar --rethrow-default -n $NEVENTS --config eventdump.fcl ${BASEFILENAME}_Reference_${STEPS[STEP]}_${LARSOFT_REFERENCE_VERSION}.root > ${BASEFILENAME}_Reference_${STEPS[STEP]}.dump
        OUTPUT_REFERENCE=$(cat ${BASEFILENAME}_Reference_${STEPS[STEP]}.dump | sed -e  '/PROCESS NAME/,/^\s*$/!d ; s/PROCESS NAME.*$// ; /^\s*$/d' )

        #echo lar --rethrow-default -n $NEVENTS --config eventdump.fcl ${BASEFILENAME}_Current_${STEPS[STEP]}.root

        lar --rethrow-default -n $NEVENTS --config eventdump.fcl ${BASEFILENAME}_Current_${STEPS[STEP]}.root > ${BASEFILENAME}_Current_${STEPS[STEP]}.dump
        OUTPUT_CURRENT=$(cat ${BASEFILENAME}_Current_${STEPS[STEP]}.dump | sed -e  '/PROCESS NAME/,/^\s*$/!d ; s/PROCESS NAME.*$// ; /^\s*$/d' )

        echo -e "\nCompare data products."
        echo    "Reference files for ${STEPS[STEP]} step generated using LARSOFT_VERSION $LARSOFT_REFERENCE_VERSION"
        echo -e "Current files for ${STEPS[STEP]} step generated using LARSOFT_VERSION $LARSOFT_VERSION\n"
        echo -e "\n${BASEFILENAME}_Reference_${STEPS[STEP]}.dump\n"
        echo "$OUTPUT_REFERENCE"
        echo -e "\n${BASEFILENAME}_Current_${STEPS[STEP]}.dump\n"
        echo "$OUTPUT_CURRENT"
    fi

    trap - ERR

    DIFF=$(diff  <(echo "${OUTPUT_REFERENCE}" | awk -v MYNF=$((1-$1)) 'NF{NF-=MYNF}1' | sed 's/\.//g' ) <(echo "${OUTPUT_CURRENT}" | awk -v MYNF=$((1-$1)) 'NF{NF-=MYNF}1' | sed 's/\.//g' ) )

    trap 'LASTERR=$?; echo -e "\nCI MSG BEGIN\n `basename $0`: error at line ${LINENO}\n Stage: ${STEPS[STEP]}\n Task: ${TASKSTRING}\n exit status: ${LASTERR}\nCI MSG END\n"; exit ${LASTERR}' ERR

    echo
    echo ${CHECKMSG[$1]}
    echo "difference(s)"
    echo

    if [ ${#DIFF} -gt 0 ]; then
       echo "${DIFF}"
       exitstatus $((2000+$1)) "$FUNCNAME step $1";
    else
       echo -e "none\n\n"
    fi

    COMPAREINIT=1

}

function exitstatus
{
    EXITSTATUS=$1

    echo -e "\nCI MSG BEGIN\n Stage: ${STEPS[STEP]}\n Task: ${TASKSTRING}\n exit status: ${EXITSTATUS}\nCI MSG END\n"
    if [ ${EXITSTATUS} -ne 0 ]; then
      exit ${EXITSTATUS}
    fi
}

trap 'LASTERR=$?; echo -e "\nCI MSG BEGIN\n `basename $0`: error at line ${LINENO}\n Stage: ${STEPS[STEP]}\n Task: ${TASKSTRING}\n exit status: ${LASTERR}\nCI MSG END\n"; exit ${LASTERR}' ERR

TASKSTRING="initialize"
initialize $@

exitstatus $?

TASKSTRING="larsoft_data_production"
larsoft_data_production

exitstatus $?


COMPAREINIT=0
TASKSTRING="compare_data_products 0"
compare_data_products 0 #Check for added/removed data products

exitstatus $?

TASKSTRING="compare_data_products 1"
compare_data_products 1 #Check for differences in the size of data products

exitstatus $?
