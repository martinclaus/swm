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
binDir="$modelDir/bin"

# List of output files to compare:
testOutputPrefix="$testDir/output/new_000000000001_"
referenceOutputPrefix="$testDir/output/reference_000000000001_"
referenceDatasetList=$(ls $referenceOutputPrefix*)

# Delete previous test results
rm -f $testOutputPrefix*

# Execute model binary from test dir
cd $testDir
time $binDir/model

# Compare results
if [ $? -eq 0 ]
  then
    for dataset in $referenceDatasetList
    do
      CDOOUTPUT=$(cdo diff $dataset ${dataset/reference/new} 2>&1)
      N_MATCH=$(echo $CDOOUTPUT | awk '/differ/ {print $(NF-4)}')
      [[ $N_MATCH = "" ]] && N_MATCH=0
      if [[ $N_MATCH != 0 ]]
      then
          echo -e "\e[00;31mDifferences in ${dataset##$referenceOutputPrefix}\e[0m"
          echo "$CDOOUTPUT"
      else
	        echo -e "\e[00;32mNo Difference in ${dataset##$referenceOutputPrefix}\e[0m"
      fi
    done
fi
