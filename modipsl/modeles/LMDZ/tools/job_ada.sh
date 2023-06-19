#!/bin/bash


##############################
# Gestion des chemins d'acces
##############################
localdir=`pwd`
dir_tools=`echo $0 | sed -e s:\`basename $0\`::`
modipsl=`cd $dir_tools ; ( pwd | sed -e 's:modipsl.*.$:modipsl:' )`
if [ ! -d $modipsl/modeles/LMDZ5 ]  ; then
    echo Dans $0
    echo Probleme pour indentifier l arborescence modipsl
    exit
fi


##############################
# Lien local pour l'executable
##############################
if [ -f gcm.e ] ; then cp gcm.e gcm.$$ ; fi
ln -sf $3 gcm.e


###################
# Creation d'un job
###################
cat <<eod>| tmp
# @ job_type = parallel
# @ job_name = NOMSIMU
# @ output   = \$(job_name).\$(jobid)
# @ error    = \$(job_name).\$(jobid)
# @ total_tasks = $1
# @ parallel_threads = $2
# @ wall_clock_limit = 00:30:00
# @ queue

ulimit -s unlimited
export OMP_STACKSIZE=800M
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/smplocal/pub/NetCDF/4.1.3/lib:/smplocal/pub/HDF5/1.8.9/seq/lib

cd $localdir
\rm -f hist*

poe  ./gcm.e

for type in day hf mth ; do
   \rm -f \$type.nc
   $modipsl/bin/rebuild -o \$type.nc hist\${type}_*nc
done

eod

###################
# soumission du job
###################
llsubmit tmp

