#!/bin/bash

#Constants
COMPILE_PATH="build"
RUN_PATH="/$COMPILE_PATH/src/filters"
TEST_PATH="/$COMPILE_PATH/test/unit_tests"

#INPUT OPTIONS
RUN_OPTIONS=$1

#If the COMPILE doesn't exist
if [ ! -d $COMPILE_PATH ]; then
	mkdir $COMPILE_PATH;
fi

#Make the file
cd $COMPILE_PATH;
cmake ../; make -j4 #2>> error_output.txt;
cd ../;

#Run the code
if [ $RUN_OPTIONS = 1 ]
then
	.$RUN_PATH;
fi
