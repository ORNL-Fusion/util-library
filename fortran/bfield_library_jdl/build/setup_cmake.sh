#!/bin/bash
#
# setup_cmake chooses the cmake build command for a machine
#
#

BUILD_TYPE=Release
VERBOSE_BUILD=0

MACHINE_ID=`uname -n`

echo Building for machine $MACHINE_ID
echo

if [ $# -eq 1 ]
then
	BUILD_TYPE=Debug
fi
echo cmake configured to generate a $BUILD_TYPE build.
if [ "$BUILD_TYPE" == Debug ]
then
	echo "    cmake may be reconfigured to generate a Release build by running this script with no arguments or using the command"	
	echo
	echo "    cmake -DCMAKE_BUILD_TYPE=Release"
else
	echo "    cmake may be reconfigured to generate a Debug build by running this script with a Debug argument or using the command"
	echo
	echo "    cmake -DCMAKE_BUILD_TYPE=Debug"
fi

rm -rf CMakeFiles CMakeCache.txt

echo 
if [ "$MACHINE_ID" == "megabucky" ]
then
   # megabucky 
   # module gcc should be loaded
   # cmake by default will take old gcc from /usr/bin so define manually
    cmake -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
	  -DCMAKE_Fortran_COMPILER=gfortran     \
	  -DCMAKE_C_COMPILER=gcc                \
	  -DCMAKE_CXX_COMPILER=g++              \
	  ..
elif [ "$MACHINE_ID" == "fusion2.ornl.gov" ]
then
    # fusion2
    #
    cmake -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
	  ..
elif [ "$MACHINE_ID" == "fusion3" ]
then
    # fusion3
    #
    cmake -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
	  ..
elif [ "$MACHINE_ID" == "ultrabucky" ]
then
    # fusion3
    #
    cmake -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
	  ..
else
    echo $MACHINE_ID is not supported by this script.
    echo Please add your machine.
    echo
fi

if [ $VERBOSE_BUILD -eq 1 ]
then
    make VERBOSE=$VERBOSE_BUILD
else
    make
fi

