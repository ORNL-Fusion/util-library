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
	echo "    cmake may be reconfigured to generate a Release build by running this script with no arguments or using the commmand"	
	echo
	echo "    cmake -DCMAKE_BUILD_TYPE=Release"
else
	echo "    cmake may be reconfigured to generate a Debug build by running this script with a Debug argument or using the commmand"
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
	        -DCMAKE_Fortran_COMPILER=gfortran         \
	        -DCMAKE_C_COMPILER=gcc                          \
	        -DCMAKE_CXX_COMPILER=g++                    \
	        ..
fi

if [$VERBOSE_BUILD==1]; then
    make VERBOSE=$VERBOSE_BUILD
else
    make
fi
make install



#

#
#   cmake -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
#               
#               
#              
#elif [ "MACHINE_ID" == "fusion2"]; then
#
#     # fusion2
#
#   cmake -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE  
#fi
#
#else
# MACHINE_ID is new and unknown. Inform the user how to add support for this new machine.
#	echo $MACHINE_ID not suported by this script.
#	echo To support this machine, add a new elif statement of the form
#	echo
#	echo elif [ \$MACHINE_ID == \"$MACHINE_ID\" ]
#	echo then
#	echo "    cmake -DVARIABLE=value ...
#	echo
#fi


