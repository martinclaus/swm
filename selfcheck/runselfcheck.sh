#!/bin/bash

# arguments
# $1 test directory

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
testDir="$selfcheckDir/$1"
buildDir="$modelDir/build"

# List of output files to compare:
testOutputPrefix="$testDir/output/new_000000000001_"
referenceOutputPrefix="$testDir/output/reference_000000000001_"
referenceDatasetList=$(ls $referenceOutputPrefix*)

# Delete previous test results
rm -f $testOutputPrefix*

# Execute model binary from test dir
cd $testDir
time $buildDir/./model

# Compare results
set -x
for dataset in $referenceDatasetList
do
cdo diff $dataset $(echo $dataset | sed 's/reference/new/')
done
