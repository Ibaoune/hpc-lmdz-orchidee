#!/bin/bash

export OMP_STACKSIZE=800M  # Max stack size for open MP
export OMP_NUM_THREADS=$2  # Number of openMP processes (threads)
ulimit -s unlimited        # Allowing illimted memory for ...
mkdir -p WORK$$ ; mv *_0*nc WORK$$ # saving old files

dir_tools=`echo $0 | sed -e s:\`basename $0\`::`
modipsl=`cd $dir_tools ; ( pwd | sed -e 's:modipsl.*.$:modipsl:' )`
if [ ! -d $modipsl/modeles/LMDZ5 ]  ; then
   echo Dans $0
   echo Probleme pour indentifier l arborescence modipsl
   exit
fi

# On recupere le chemin pour mpirun dans arch.fcm
mpirun=`head -1 $modipsl/modeles/LMDZ5/arch.fcm | awk ' { print $2 } ' | sed -e 's/mpif90/mpirun/'`
$mpirun  -np $1 $3

#############################################################################
# Reconstructing netcdf files from fragments comming from parallel execution
# file_0000.nc file_0001.nc ... > file.nc
# Using the rebuild commend available on ..../modipsl/bin/rebuild
#############################################################################

for f in `ls *_0000.nc` ; do
   file=`echo $f | sed -e 's/_0000.nc//'` ; echo Rebuild for $file
   if [ -f $file.nc ] ; then
       echo The $file.nc file already exists
       echo remove it and run 
       echo $modipsl/bin/rebuild -o $file.nc ${file}_0*.nc
       echo manually
   else
      $modipsl/bin/rebuild -o $file.nc ${file}_0*.nc
      if [ -f $file.nc ] ; then \rm -f ${file}_0*.nc ; fi
   fi
done
