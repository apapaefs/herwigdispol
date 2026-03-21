#!/bin/bash
# script to put the HERWIG and THEPEG paths into the simple makefiles,
# for everything else the user is on their own!
# loop over all the files

OPENLOOPSPREFIX=""
THEPEGINCLUDE="-I/home/apapaefs/Projects/Herwig/Herwig-pol-full-python3-rivet4/include"
GSLINCLUDE="-I/home/apapaefs/Projects/Herwig/Herwig-pol-full-python3-rivet4/include"
FASTJETINCLUDE="-I/home/apapaefs/Projects/Herwig/Herwig-pol-full-python3-rivet4/include"
HERWIGINCLUDE=-I/home/apapaefs/Projects/Herwig/Herwig-pol-full-python3-rivet4/include/
HERWIGINSTALL="/home/apapaefs/Projects/Herwig/Herwig-pol-full-python3-rivet4"
FC="gfortran"
FCLIBS=" -L/opt/intel/oneapi/tcm/1.4/lib/../lib -L/opt/intel/oneapi/umf/1.0/lib/../lib -L/opt/intel/oneapi/mpi/2021.17/lib/../lib -L/opt/intel/oneapi/compiler/2025.3/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/gcc/x86_64-linux-gnu/9/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/9/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/opt/intel/oneapi/tcm/1.4/lib -L/opt/intel/oneapi/umf/1.0/lib -L/opt/intel/oneapi/mpi/2021.17/lib -L/opt/intel/oneapi/compiler/2025.3/lib -L/usr/lib/gcc/x86_64-linux-gnu/9/../../.. -lgfortran -lm -lquadmath"
CXX="g++ -std=c++17"
CXXFLAGS="-O2 -DBOOST_UBLAS_NDEBUG -Wno-deprecated-declarations -Wno-deprecated-copy"
LDFLAGS="" 
SHARED_FLAG="-shared " 

for i in *
do
# if a directory
  if [ -d $i ]; then
# check input files exists
      file=$i/Makefile.in
      if [ -e $file ]; then
	  file2=`echo $file | sed s!\.in!!`
	  echo 'Making ' $file2
	  sed "s!THEPEGINCLUDE *\=!THEPEGINCLUDE=$THEPEGINCLUDE!" < $file | \
	  sed "s!OPENLOOPSPREFIX *\=!OPENLOOPSPREFIX=$OPENLOOPSPREFIX!" | \
	  sed "s!FC *\=!FC=$FC!" | \
	  sed "s!FCLIBS *\=!FCLIBS=$FCLIBS!" | \
          sed "s!CXX *\=!CXX=$CXX!" | \
          sed "s!SHARED_FLAG *\=!SHARED_FLAG=$SHARED_FLAG!" | \
          sed "s!LDFLAGS *\=!LDFLAGS=$LDFLAGS!" | \
          sed "s!CXXFLAGS *\=!CXXFLAGS=$CXXFLAGS!" | \
	  sed "s!HERWIGINCLUDE *\=!HERWIGINCLUDE=$HERWIGINCLUDE!" | \
          sed "s!HERWIGINSTALL *\=!HERWIGINSTALL=$HERWIGINSTALL!" | \
	  sed "s!GSLINCLUDE *\=!GSLINCLUDE=$GSLINCLUDE!" | \
	  sed "s!FASTJETINCLUDE *\=!FASTJETINCLUDE=$FASTJETINCLUDE!" > $file2
      fi
  fi
done
