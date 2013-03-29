#!/bin/bash

# Move to project dir, e.g. /home/zzz/model
# Only works from model or selfcheck directory
if [ -L $0 ] ; then
    ME=$(readlink $0)
else
    ME=$0
fi
DIR=$(dirname $ME)
if [ "$DIR" = "." ]; then
    cd .. 
fi

# Save paths
modelDir="$PWD"
selfcheckDir="$modelDir/selfcheck"
buildDir="$modelDir/build"

# List of output files to compare:
testOutputPrefix="$selfcheckDir/output/new_000000000001_"
referenceOutputPrefix="$selfcheckDir/output/reference_000000000001_"

fileSuffixList="eta_out.nc psi_out.nc u_out.nc v_out.nc"

# Delete previous test results
for fileSuffix in ${fileSuffixList} 
do
rm -f $testOutputPrefix$fileSuffix
done

# Execute model binary from selfcheck dir
cd $selfcheckDir
time $buildDir/./model

# Compare results
set -x
for fileSuffix in ${fileSuffixList} 
do
cdo diff $referenceOutputPrefix$fileSuffix $testOutputPrefix$fileSuffix
done

