#!/bin/bash
#set -vx
# Default options
#
#
#
date
#### 1  Set up the compiling options
#### Define some directories
submitdir=$( pwd )
modipsl=$submitdir/../..
arch_path=$submitdir/ARCH 

#### Set default options 
# Optimization mode
# optmode=prod/dev/debug
optmode=prod
# Resolution of LMDZ in regular mode
#resol_atm=144x142x79
resol_atm=384x220x79
# Also compile ce0l subprogram to LMDZ (y/n)
ce0l=y

# fcm_arch
fcm_arch=default
# Default values to be overritten
parallel=mpi_omp
export fcm_arch parallel xios
full_flag=""
full_xios=""
full_lmdz=""
full_orch=""

# Default netcdf_lib is used for XIOS but can be change by argument
netcdf_lib=""

# Output text file for compilation of each component
datestr=`LC_ALL=C date +"%Y%m%dT%H%M"`
outfile=$submitdir/out_compile_lmdzor.$datestr
echo > $outfile
echo; echo "Text output from compilation will be stored in file out_compile_lmdzor.$datestr"; echo 

#### Read arguments 
# Loop over all arguments to modify default set up
while (($# > 0)) ; do
    case $1 in
	"-h") cat <<end_help
########################################################################
# Usage of the script compile_lmdzor.sh
#
########################################################################

./compile_config [Options] 

Options: -full, -resol_atm, -debug, -dev, -prod(default)

Example 1: Default compilation of LMDZ-ORCHIDEE with XIOS and IOIPSL
./compile_lmdzor.sh

Example 2: Compile in debug mode
./compile_lmdzor.sh -debug

Example 3: Default compilation with full recompilation of all components
./compile_lmdzor.sh -full

Example 4: Compilation of LMDZ in regular lat-lon for other dimensions (default: 144x142x79) 
The dimension can be changed to any other 3d dimension.
./compile_lmdzor.sh -resol_atm 96x95x79

end_help
exit;;
	"-parallel") 	    parallel=$2 ; shift ; shift ;; 
	"-arch")	    fcm_arch="$2" ; shift ; shift ;;
	"-xios") 	    xios="$2" ; shift ; shift ;;
	"-debug")	    optmode=debug ; shift ;;
	"-dev") 	    optmode=dev ; shift ;;
	"-prod")	    optmode=prod ; shift ;;
	"-resol_atm")	    resol_atm=$2 ; shift ; shift ;;
	"CE0L")             ce0l=y ; shift ;;
	"-full")            full_flag="-full"; full_xios="--full" ; full_lmdz="-full" ; full_orch="-full" ; shift ;;
	"-full_xios")       full_xios="--full" ; shift ;;  # Note only full_xios is using double dash: --full
	"-full_lmdz")       full_lmdz="-full"  ; shift ;;
	"-full_orch")       full_orch="-full"  ; shift ;;
	"-netcdf_lib_seq")  netcdf_lib="--netcdf_lib netcdf4_seq"; shift ;;
	*)	            echo "unknown option "$1" , exiting..." ; exit 
    esac
done

echo "Following options are set in current compiling:" >> $outfile 
echo "   resol_atm=${resol_atm}" >> $outfile 
echo "   ce0l=${ce0l}" >> $outfile 
echo "   optmode = $optmode, parallel = $parallel, fcm_arch = $fcm_arch " >> $outfile 
echo "   full_flag=$full_flag, full_xios=$full_xios, full_lmdz=$full_lmdz, full_orch=$full_orch" >> $outfile 
echo >> $outfile

### Read host dependent default values
### These variables will not be changed if they were set as argument
###./host.sh $host
# Later : Following lines should be set in host.sh file
# begin host.sh
if [ $fcm_arch == default ] ; then
    # Find out current host and source specific paths and commands for the host
    case $( hostname -s ) in
	jean-zay*)
	    fcm_arch=X64_JEANZAY;;
        irene170|irene171|irene190|irene191|irene192|irene193)
            fcm_arch=X64_IRENE;;
        irene172|irene173|irene194|irene195)
            fcm_arch=X64_IRENE-AMD;;
	asterix*|obelix*)
	    fcm_arch=ifort_LSCE
            if [ "${parallel}" == "mpi_omp" ] ; then 
		echo "Warning!! Currently at Obelix hybrid mode is not possible."
		echo "          Option -parallel mpi_omp is now changed to -parallel mpi"
		echo ""
		parallel=mpi
            else
		echo "You used -parallel" $parallel
            fi ;;
	ciclad*|climserv*|loholt*|camelot*)
	    fcm_arch=ifort_CICLAD;;
	*)
	    echo Current host is not known. You must use option -arch to specify which architecuture files to use.
	    echo Exit now.
	    exit
    esac
fi

# Set a link to arch.env if arch-${fcm_arch}.env file exist for current fcm_arch.
# The link arch.env is also set in config.card and will be used by libIGCM to ensure the same running environnement.
if [ -f ARCH/arch-${fcm_arch}.env ] ; then
    echo >> $outfile
    echo "The file ARCH/arch-${fcm_arch}.env will now be sourced with modules needed for compilation for all components."
    echo "Note that this new environement might be kept after compilation." 
    echo "If this is the case, source again your personal environment after compilation. "
    echo " Personal module list before sourcing of ARCH/arch.env file:"    >> $outfile  
    module list   >> $outfile 2>&1 

    # Make a link to this file, to be used also in config.card
    rm -f ARCH/arch.env
    ln -s arch-${fcm_arch}.env ARCH/arch.env

    # Source the file
    source ARCH/arch.env   >> $outfile 2>&1 
    echo >> $outfile 
    echo " New module list after sourcing of ARCH/arch.env file:"    >> $outfile  
    module list   >> $outfile 2>&1 
fi

#### 2 Do the compilation 
## 2.1 Compile ioipsl
cd $modipsl/modeles/IOIPSL
echo; echo "NOW COMPILE IOIPSL"
echo >> $outfile ; echo " NOW COMPILE IOIPSL"   >> $outfile 

# Check if compilation with fcm is included in IOIPSL
if [ ! -f makeioipsl_fcm ] ; then
  echo "The file makeiopsl_fcm do not exist. Probably this is not the right version of IOIPSL."
  echo "Version of IOIPSL with fcm compilation is needed to compile with this script. Stop now."
  exit
fi

echo ./makeioipsl_fcm -$optmode -arch ${fcm_arch} -arch_path $arch_path -j 8 $full_flag   >> $outfile 
     ./makeioipsl_fcm -$optmode -arch ${fcm_arch} -arch_path $arch_path -j 8 $full_flag   >> $outfile 2>&1
# Test if compiling succeded 
if [[ $? != 0 ]] ; then 
    echo "THERE IS A PROBLEM IN IOIPSL COMPILATION - STOP"
    exit
fi

## 2.2 Compile xios
cd $modipsl/modeles/XIOS
echo; echo "NOW COMPILE XIOS"
echo >> $outfile ; echo " NOW COMPILE XIOS"   >> $outfile 
echo ./make_xios --$optmode --arch $fcm_arch --arch_path $arch_path ${netcdf_lib} --job 4 $full_xios   >> $outfile 
     ./make_xios --$optmode --arch $fcm_arch --arch_path $arch_path ${netcdf_lib} --job 4 $full_xios   >> $outfile 2>&1
# Test if compiling succeded 
if [[ $? != 0 ]] ; then 
    echo "THERE IS A PROBLEM IN XIOS COMPILATION - STOP"
    exit
fi
# Move executables to modipsl/bin
if [ -f $modipsl/modeles/XIOS/bin/xios_server.exe ] ; then 
    mv $modipsl/modeles/XIOS/bin/xios_server.exe $modipsl/bin/xios_server_${optmode}.exe
else
    echo "THERE IS A PROBLEM IN XIOS COMPILATION EXECUTABLE MISSING - STOP"
    exit
fi


## 2.3 Compile orchidee
cd $modipsl/modeles/ORCHIDEE
echo; echo "NOW COMPILE ORCHIDEE"
echo >> $outfile ; echo " NOW COMPILE ORCHIDEE"   >> $outfile 


# Check if the compilation of ORCHIDEE was previsouly interupted prematured. 
# In that case, the files fcm.bld.lock exist. If the file exist, it is here removed to avoid makeorchidee_fcm 
# to ask question and wait for interactivly answer from the user.
if [ -f build/fcm.bld.lock ] ; then
    echo >> $outfile
    echo " WARNING!! The file ORCHIDEE/build/fcm.bld.lock exist. "   >> $outfile
    echo "           This means that the compilation is either currently on going in another terminal was previous interupted before the end."   >> $outfile
    echo "           The file will now be removed. "   >> $outfile
    echo >> $outfile
    rm -f build/fcm.bld.lock
fi

echo 	./makeorchidee_fcm -j 8 -xios -parallel $parallel -$optmode -arch ${fcm_arch} -arch_path $arch_path $full_orch -driver   >> $outfile 
        ./makeorchidee_fcm -j 8 -xios -parallel $parallel -$optmode -arch ${fcm_arch} -arch_path $arch_path $full_orch    >> $outfile 2>&1
# Test if compiling finished
if [[ $? != 0 ]] ; then
    echo "THERE IS A PROBLEM IN ORCHIDEE COMPILATION - STOP"
    exit
fi

# Rename executables to contain $optmode
if [ -f $modipsl/bin/orchidee_ol ] ; then mv $modipsl/bin/orchidee_ol  $modipsl/bin/orchidee_ol_${optmode} ; fi
if [ -f $modipsl/bin/orchideedriver ] ; then mv $modipsl/bin/orchideedriver  $modipsl/bin/orchideedriver_${optmode} ; fi


## 2.4 Compile LMDZ for regular latlon configuration
cd $modipsl/modeles/LMDZ
# Compile LMDZ regular lat-lon executable
echo; echo "NOW COMPILE LMDZ REGULAR LAT-LON MODE. Resolution = ${resol_atm}"
echo >> $outfile ; echo " NOW COMPILE LMDZ REGULAR LAT-LON MODE. Resolution = ${resol_atm}"   >> $outfile 

# Check if the compilation of LMDZ was previsouly interupted prematured. 
# In that case, the files dimension.h and .lock exist. If the files exist, they are here removed to avoid makelmdz_fcm 
# to ask question and wait for interactivly answer from the user.
if [ -f libf/grid/dimensions.h ] || [ -f .lock ] ; then
  echo >> $outfile
  echo " WARNING!! The file LMDZ/libf/grid/dimension.h and/or the LMDZ/.lock exist. "   >> $outfile
  echo "           This means that the compilation is either currently on going in another terminal was previous interupted before the end."   >> $outfile
  echo "           The files dimension.h and .lock will now be removed. "   >> $outfile
  echo >> $outfile
  rm -f libf/grid/dimensions.h
  rm -f .lock
fi

echo ./makelmdz_fcm -d ${resol_atm} -p lmd -rrtm true -$optmode -mem -parallel $parallel -io xios -v orchidee2.1 -arch $fcm_arch -arch_path $arch_path -j 8 $full_lmdz gcm    >> $outfile 
./makelmdz_fcm -d ${resol_atm} -p lmd -rrtm true -$optmode -mem -parallel $parallel -io xios -v orchidee2.1 -arch $fcm_arch -arch_path $arch_path -j 8 $full_lmdz gcm    >> $outfile 2>&1
# Test if compiling finished
if [[ $? != 0 ]] ; then
    echo "THERE IS A PROBLEM IN LMDZ REGULAR LATLON COMPILATION - STOP"
    exit
fi

# Move executables to modipsl/bin folder
echo >> $outfile
echo "Move gcm.e executable to modipsl/bin" >> $outfile
echo ls -lrt $modipsl/modeles/LMDZ/bin   >> $outfile
ls -lrt $modipsl/modeles/LMDZ/bin  >> $outfile
echo >> $outfile

# Find executable suffix
if [ $parallel == seq ] || [ $parallel == none ] ; then
    suffix=_${resol_atm}_phylmd_seq_orch.e
else
    suffix=_${resol_atm}_phylmd_para_mem_orch.e
fi
echo gcm suffix = $suffix

if [ -f $modipsl/modeles/LMDZ/bin/gcm${suffix} ] ;  then 
    mv $modipsl/modeles/LMDZ/bin/gcm${suffix} $modipsl/bin/gcm_${resol_atm}_${optmode}.e
else
    echo "ERROR gcm${suffix} executable does not exist." 
    echo "THERE IS A PROBLEM IN LMDZ REGULAR LATLON COMPILATION - STOP"
    exit
fi  


# Compile ce0l initialization program for LMDZ regular lat-lon exectuable
if [ $ce0l == y ] ; then
    echo; echo "NOW COMPILE CE0L OF LMDZ. Resolution = ${resol_atm}"
    echo >> $outfile ; echo " NOW COMPILE CE0L OF LMDZ. Resolution = ${resol_atm}"   >> $outfile 
    
    echo ./makelmdz_fcm -d ${resol_atm} -p lmd -rrtm true -$optmode -io xios -arch $fcm_arch -arch_path $arch_path -j 8 $full_lmdz ce0l    >> $outfile 
    ./makelmdz_fcm -d ${resol_atm} -p lmd -rrtm true -$optmode -io xios -arch $fcm_arch -arch_path $arch_path -j 8 $full_lmdz ce0l    >> $outfile 2>&1
    # Test if compiling finished
    if [[ $? != 0 ]] ; then
	echo "THERE IS A PROBLEM IN CE0L (LMDZ) REGULAR LATLON COMPILATION - STOP"
	exit
    fi

    
    # Move executables to modipsl/bin folder
    echo >> $outfile
    echo "Move ce0l.e executable to modipsl/bin" >> $outfile
    echo ls -lrt $modipsl/modeles/LMDZ/bin   >> $outfile
    ls -lrt $modipsl/modeles/LMDZ/bin  >> $outfile
    echo >> $outfile
    suffix=_${resol_atm}_phylmd_seq.e
    echo ce0l suffix = $suffix

    if [ -f $modipsl/modeles/LMDZ/bin/ce0l${suffix} ] ;  then 
	mv $modipsl/modeles/LMDZ/bin/ce0l${suffix} $modipsl/bin/ce0l_${resol_atm}_${optmode}.e
    else
	echo "ERROR ce0l${suffix} executable does not exist." 
	echo "THERE IS A PROBLEM IN CE0L (LMDZ) REGULAR LATLON COMPILATION - STOP"
	exit
    fi  
fi


echo >>$outfile ; echo "ALL COMPILING FINISHED" >> $outfile
echo ls -lrt modipsl/bin >> $outfile
ls -lrt $modipsl/bin >> $outfile

echo; echo "ALL COMPILING FINISHED" ; echo
echo "Executables are found in modipsl/bin"
echo "Check that executable names correspond with the name set in config.card before launching the job"
echo ls -lrt modipsl/bin
ls -lrt $modipsl/bin

date

exit
