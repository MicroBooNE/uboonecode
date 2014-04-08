#! /bin/bash

#------------------------------------------------------------------
#
# Purpose: Generate a sam wrapper fcl file.
#
# Usage: make_sam_wrapper.sh <job-fcl> <project-url> <consumer-process-id>
#
# Notes:
#
# 1.  This script is invoked with three positional arguments, all 
#     of which are required.
#
# 2.  The wrapper fcl file is written to standard output.  This
#     output should normally be captured to a file.
#
# Created: H. Greenlee, 30-Jul-2013
#
#------------------------------------------------------------------

jobfcl=''
prjurl=''
cpid=''

if [ $# -ne 3 ]; then
  echo "Usage: make_sam_wrapper.sh <job-fcl> <project-url> <consumer-process-id>" >&2
  exit 1
fi
jobfcl=$1
prjurl=$2
cpid=$3

cat <<EOF
#include "${jobfcl}"

services.user.IFDH:
{
  IFDH_BASE_URI: "http://samweb.fnal.gov:8480/sam/uboone/api"
}

services.user.CatalogInterface:
{
  service_provider: "IFCatalogInterface"
  webURI: "${prjurl}"
}

services.user.FileTransfer:
{
  service_provider: "IFFileTransfer"
}

source.fileNames: [ "${cpid}" ]

EOF
