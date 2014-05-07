#!/bin/bash
#
# Lists all the FCL files whose name matches the specified pattern.
# 
# Usage:  ListFCL.sh [Pattern]
# 
# The pattern is a regex as interpreted by grep ("-e" option).
#

SCRIPTNAME="$(basename "$0")"

: ${FORMAT:="%f (in %h)"}

function help() {
	cat <<-EOH
	Looks for FHICL files in the ART FHICL search directories.
	
	Usage:  ${SCRIPTNAME} [options] [--] [Pattern]
	
	Options:
	--format=FORMAT
	    use the specified format; the default format is "%f (in %h)";
	    the placeholder are as documented in find (1) manual
	--name , -n
	    print the file name only
	--path , -p
	    print the full path
	
	EOH
} # help()


function isFlagSet() {
	local VarName="$1"
	[[ -n "${!VarName//0}" ]]
} # isFlagSet()

function STDERR() { echo "$*" >&2 ; }
function ERROR() { STDERR "ERROR: $@" ; }
function FATAL() {
	local Code="$1"
	shift
	STDERR "FATAL ERROR (${Code}): $*"
	exit $Code
} # FATAL()
function LASTFATAL() {
	local Code="$?"
	[[ "$Code" != 0 ]] && FATAL "$Code""$@"
} # LASTFATAL()

function Filter() {
	local Pattern="$1"
	if [[ -n "$Pattern" ]]; then
		grep -e "$Pattern"
	else
		cat
	fi
} # Filter()

function CleanSortKey() { sed -e 's/^[^ ]* \+//' ; }

function ApplyFormat() {
	local Format="$1"
	local Data
	while read Data ; do
		local FileName="${Data%%/*}"
		local Path="${Data#*/}"
		
		case "$Format" in
			( '' )
				echo "${FileName} (in ${Path})"
				;;
			( * )
				sed -e "s@\([^\]|\^\)%f@${FileName}@" -e "s@\([^\]|\^\)%h@${Path}@" <<< "$Format"
		esac
	done
}

################################################################################

declare -i NoMoreOptions=0
declare -a Patterns
declare -i nPatterns=0
for (( iParam = 1 ; iParam <= $# ; ++iParam )); do
	Param="${!iParam}"
	if ! isFlagSet NoMoreOptions && [[ "${Param:0:1}" == '-' ]]; then
		case "$Param" in
			( '--help' | '-h' | '-?' ) DoHelp=1  ;;
			
			### format options
			( '--name' | '-n' ) FORMAT="%f" ;;
			( '--path' | '-p' ) FORMAT="%p" ;;
			( '--format='* )    FORMAT="${Param#--*=}" ;;
			
			### other stuff
			( '-' | '--' )
				NoMoreOptions=1
				;;
			( * )
				echo "Unrecognized script option #${iParam} - '${Param}'"
				exit 1
				;;
		esac
	else
		NoMoreOptions=1
		Patterns[nPatterns++]="$Param"
	fi
done

if isFlagSet DoHelp ; then
	help
	# set the exit code (0 for help option, 1 for missing parameters)
	isFlagSet DoHelp
	exit $?
fi

# explanation:
# - start with paths in FHICL_FILE_PATH
# - split by ":"
# - find in all those directories (but not in their subdirectories),
#   and in that order, all the FCL files
# - print for each its name and the string to be presented to the user as output
# - soft them by FCL file name, preserving the relative order of files with the
#   same name from different directories
# - filter them on sort key (file name) by user's request
# - remove the sort key (file name) from the output
tr ':' "\n" <<< "$FHICL_FILE_PATH" | xargs -I SEARCHPATH find SEARCHPATH -maxdepth 1 -name "*.fcl" -printf "%f ${FORMAT}\n" 2> /dev/null | sort -s -k1,1 -u | Filter "^[^ ]*${Patterns[0]}[^ ]* " | CleanSortKey
