#!/bin/ksh

#**************************************************************
# Author: Sebastien Denvil
# Contact: Sebastien.Denvil__at__ipsl.jussieu.fr
# $Revision:: 1569                                     $ Revision of last commit
# $Author:: aclsce                                     $ Author of last commit
# $Date:: 2022-12-15 11:08:40 +0100 (خ, 15 دجنبر#$ Date of last commit
# IPSL (2006)
#  This software is governed by the CeCILL licence see libIGCM/libIGCM_CeCILL.LIC
#
#**************************************************************

#=========================================================
# The documentation of this file can be automatically generated
# if you use the prefix #D- for comments to be extracted.
# Extract with command: cat lib* | grep "^#D-" | cut -c "4-"
#=========================================================

#D-#==================================================
#D-LibIGCM_sys for IPSL ESPRI MESO cluster: spirit and spiritx
#D-#==================================================
#D-
#D- This ksh library if a layer under some usefull
#D-environment variables and shell commands.
#D-All those definitions depend on host particularities.
#D-It manages a stack mechanism and test validity of operations.
#D-All function described bellow must be prefixed by IGCM_sys.

#====================================================
# libIGCM_sys PARAMETERS
#====================================================

#====================================================
# set DEBUG_sys to true to output calls of function
typeset -r DEBUG_sys=${DEBUG_sys:=true}

#====================================================
# Turn in dry run mode ? (sys_Put_Rest, sys_Put_Out, sys_Get)
if ([ X${DRYRUN} != X ] && [ ${DRYRUN} = 4 ]) ; then
  typeset -r DRYRUN=0
  typeset -r DRYRUN_DEBUG=4
else
  typeset -r DRYRUN=${DRYRUN:=0}
  typeset -r DRYRUN_DEBUG=0
fi

# YOU MUST COMPILE YOUR EXE FILES FOR DRYRUN MODE !
# -------------------------------------------------------------------------------------
# | DRYRUN=  |  Date computations, | sys_Get    |  Exe    | sys_Put_Out; sys_Put_Rest |
# |          |  Cp/Exe/param/files |            |         |                           |
# |          |  Chmod Qsub         |            |         |                           |
# -------------------------------------------------------------------------------------
# |    0     |       yes           |    yes     |  yes    |      yes                  |
# -------------------------------------------------------------------------------------
# |    1     |       yes           |    yes     |  yes    |      no                   |
# -------------------------------------------------------------------------------------
# |    2     |       yes           |    yes     |  no     |      no                   |
# -------------------------------------------------------------------------------------
# |    3     |       yes           |    no      |  no     |      no                   |
# -------------------------------------------------------------------------------------
# |    4     |       yes           |    yes     |  creation of RUN_DIR and Job_debug  |
# -------------------------------------------------------------------------------------

#=====================================================
# Global Variables :
#=====================================================
# Language : "fr" or "en"
typeset -r MYLANG="fr"

#=====================================================
# Host user names project maxCpuTime
# $hostname ou hostname
typeset HOST=${HOST:=$( hostname )}
# $username ou whoami
typeset LOGIN=${LOGIN:=$( whoami )}
# $hostname of the MASTER job
typeset MASTER=${MASTER:=$( hostname )}
# project name
typeset PROJECT=NONE
# jobWarningDelay in seconds
#typeset jobWarningDelay=${PBS_WALLTIME}

#D-
#D-#==================================================
#D-Program used in libIGCM
#D-#==================================================

# Submit command
typeset SUBMIT=${SUBMIT:=sbatch}
# rsync with path
typeset -r RSYNC=/usr/bin/rsync
# RSYNC_opt args to rsync
typeset -r RSYNC_opt="-va"
# ie storage filesystem
typeset -r STOREHOST=${MASTER}

#====================================================
# Access to module command
#====================================================
. /etc/profile.d/modules.sh

#====================================================
# Set environment tools (ferret, nco, cdo, rebuild, ...)
#====================================================
if [ X${TaskType} = Xcomputing ] ; then
    IGCM_debug_Print 1 "Modules will be loaded later in IGCM_sys_activ_variables."
else
    if [ $CENTER == spirit ] ; then
        # At spirit
        . /home/igcmg/MachineEnvironment/mesoipsl/atlas_env_mesoipsl
    elif [ $CENTER == spiritx ] ; then
        # At spiritx
	. /ciclad-home/igcmg/MachineEnvironment/mesoipsl/atlas_env_mesoipsl
    fi
fi

[ ! X${TaskType} = Xchecking ] && IGCM_debug_Print 1 "List of loaded modules:"
[ ! X${TaskType} = Xchecking ] && module list

#====================================================
# Host specific DIRECTORIES
#====================================================
# For rebuild
if [ $CENTER == spirit ] ; then
    export PATH=${PATH}:/home/igcmg/rebuild/spirit/modipsl/modeles/IOIPSL/bin/
else
    export PATH=${PATH}:/ciclad-home/igcmg/rebuild/spirit/modipsl/modeles/IOIPSL/bin/
fi
#====================================================
#- MirrorlibIGCM for frontend
typeset -r MirrorlibIGCM=${MirrorlibIGCM:=false}

#====================================================
#- libIGCM_POST for frontend
typeset -r libIGCM_POST=${libIGCM}

#====================================================
#- R_EXE   (==> BIN_DIR = ${MODIPSL}/bin )
typeset -r R_EXE="${MODIPSL}/bin"

#====================================================
#- SUBMIT_DIR : submission dir
if [ X${SLURM_SUMBIT_DIR} != X ] ; then
  typeset -x SUBMIT_DIR=${SUBMIT_DIR:=${SLURM_SUBMIT_DIR}}
else
  typeset -x SUBMIT_DIR=${SUBMIT_DIR:=${PWD}}
fi

#====================================================
#- IN
typeset -r R_IN=${R_IN:=/projsu/igcmg/IGCM}

#====================================================
#- RUN_DIR_PATH : Temporary working directory

if [ $CENTER == spirit ] ; then
    RUN_DIR_PATH=${RUN_DIR_PATH:=/scratchu/${LOGIN}/RUN_DIR}
else
    RUN_DIR_PATH=${RUN_DIR_PATH:=/scratchx/${LOGIN}/RUN_DIR}
fi
if [ X${SLURM_JOBID} != X ] ; then
    typeset -r RUN_DIR_PATH=${RUN_DIR_PATH}/${SLURM_JOBID}_${$}
else
    typeset -r RUN_DIR_PATH=${RUN_DIR_PATH}/${$}
fi

#====================================================
#- OUTCOMMAND_PATH : tmp place to store command lines standard error and outputs
# Note that typeset -r can not be set at spirit because it sets the path as read-only
typeset OUTCOMMAND_PATH=/tmp

#====================================================
#- HOST_MPIRUN_COMMAND
typeset -r HOST_MPIRUN_COMMAND=${HOST_MPIRUN_COMMAND:="time mpirun"}

#====================================================
#- Max number of arguments passed to nco operator or demigration command
UNIX_MAX_LIMIT=360

#====================================================
#- set PackDefault to false on meso-ipsl 
PackDefault=false

#====================================================
#- Default number of MPI task for IPSL coupled model
#- required for backward compatibility
#-
DEFAULT_NUM_PROC_OCE=5
DEFAULT_NUM_PROC_CPL=1
(( DEFAULT_NUM_PROC_ATM = BATCH_NUM_PROC_TOT - DEFAULT_NUM_PROC_OCE - DEFAULT_NUM_PROC_CPL ))

#D-#==================================================
#D-function IGCM_sys_defineArchives
#D-* Purpose:
#D-* Define ARCHIVE : Dedicated to large files
#D-* Define STORAGE : Dedicated to small/medium files
#D-* Define R_OUT   : Output tree located on ARCHIVE
#D-* Define R_FIG   : Output tree located on STORAGE hosting figures (monitoring and atlas, and/or small files)
#D-* Define R_BUF   : Output tree located on STORAGE hosting files waiting for rebuild or pack processes
#D-* if SpaceName=TEST nothing special will hapen
#D-* Examples:
#D-
function IGCM_sys_defineArchives {
  IGCM_debug_PushStack "IGCM_sys_defineArchives"

  if [ ! X${config_UserChoices_ARCHIVE} = X ]; then
    #====================================================
    #- ARCHIVE (dedicated to large files)
    ARCHIVE=${config_UserChoices_ARCHIVE}
  else
    #====================================================
    #- ARCHIVE (dedicated to large files)
    if [ $CENTER == spirit ] ; then
	# At spirit
   	ARCHIVE=${ARCHIVE:=/data/${LOGIN}}
    else
	# At spiritx
	ARCHIVE=${ARCHIVE:=/homedata/${LOGIN}}
    fi
  fi

  if [ ! X${config_UserChoices_STORAGE} = X ]; then
    #====================================================
    #- STORAGE (dedicated to small/medium files)
    STORAGE=${config_UserChoices_STORAGE}
  else
    #====================================================
    #- STORAGE (dedicated to small/medium files)
    STORAGE=${ARCHIVE}
  fi

  # ON MESO-IPSL NO SPECIAL CASE WHEN X${config_UserChoices_SpaceName} = XTEST

  #====================================================
  #- R_OUT
  R_OUT=${ARCHIVE}/IGCM_OUT

  #====================================================
  #- R_FIG (hosting figures : monitoring and atlas, and/or small files)
  R_FIG=${STORAGE}/IGCM_OUT

  #====================================================
  #- CMIP6 (hosting CMIP6 files produced by XIOS2 and configured by dr2xml)
  CMIP6_BUF=${STORAGE}/IGCM_OUT

  #====================================================
  #- R_BUF (ONLY FOR double copy an scratch)
  R_BUF=${STORAGE}/IGCM_OUT

  IGCM_debug_Print 1 "R_OUT has been defined = ${R_OUT}"
  IGCM_debug_Print 1 "R_BUF has been defined = ${R_BUF}"
  IGCM_debug_Print 1 "R_FIG has been defined = ${R_FIG}"

  IGCM_debug_PopStack "IGCM_sys_defineArchives"
}

#D-#==================================================
#D-function IGCM_sys_RshArchive
#D-* Purpose: Archive rsh command
#D-* Examples:
#D-
function IGCM_sys_RshArchive {
  IGCM_debug_PushStack "IGCM_sys_RshArchive" $@
  /bin/ksh <<-EOF
    ${@}
EOF
  status=$?
  if [ ${status} -gt 0 ] ; then
    IGCM_debug_Print 2 "IGCM_sys_RshArchive : command failed error code ${status}"
    IGCM_debug_Exit "IGCM_sys_RshArchive"
  fi
  IGCM_debug_PopStack "IGCM_sys_RshArchive"
}

#D-#==================================================
#D-function IGCM_sys_RshArchive_NoError
#D-* Purpose: Archive rsh command, without error
#D-*          used only in monitoring.job
#D-* Examples:
#D-
function IGCM_sys_RshArchive_NoError {
  IGCM_debug_PushStack "IGCM_sys_RshArchive_NoError" $@
  /bin/ksh <<-EOF
    ${@} 2> /dev/null
EOF
  IGCM_debug_PopStack "IGCM_sys_RshArchive_NoError"
}

#D-#==================================================
#D-function IGCM_sys_MkdirArchive
#D-* Purpose: Mkdir on Archive
#D-* Examples:
#D-
function IGCM_sys_MkdirArchive {
  IGCM_debug_PushStack "IGCM_sys_MkdirArchive" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_MkdirArchive :" $@
  fi
  #- creation de repertoire sur le serveur fichier
  if [ ! -d ${1} ]; then
    \mkdir -p $1
    status=$?

    if [ ${status} -gt 0 ] ; then
      IGCM_debug_Print 2 "IGCM_sys_MkdirArchive : mkdir failed error code ${status}"
      IGCM_debug_Exit "IGCM_sys_MkdirArchive"
    fi
  fi
  IGCM_debug_PopStack "IGCM_sys_MkdirArchive"
}

#D-#==================================================
#D-function IGCM_sys_TestDirArchive
#D-* Purpose: Test Directory that must exists on Archive
#D-* Examples:
#D-
function IGCM_sys_TestDirArchive {
  IGCM_debug_PushStack "IGCM_sys_TestDirArchive" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_TestDirArchive :" $@
  fi
  typeset ExistFlag
  ExistFlag=$( [ -d $1 ] && echo 0 || echo 1 )
  IGCM_debug_PopStack "IGCM_sys_TestDirArchive"
  return ${ExistFlag}
}

#D-#==================================================
#D-function IGCM_sys_IsFileArchived
#D-* Purpose: Test file that must NOT EXISTS on Archive based on filename only
#D-* Examples:
#D-
function IGCM_sys_IsFileArchived {
  IGCM_debug_PushStack "IGCM_sys_IsFileArchived" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_IsFileArchived :" $@
  fi
  typeset IsArchivedFlag
  # Never archived for this system
  IsArchivedFlag=1
  IGCM_debug_PopStack "IGCM_sys_IsFileArchived"

  return ${IsArchivedFlag}
}

#D-#==================================================
#D-function IGCM_sys_TestFileArchive
#D-* Purpose: Test file that must NOT EXISTS on Archive
#D-* Examples:
#D-
function IGCM_sys_TestFileArchive {
  IGCM_debug_PushStack "IGCM_sys_TestFileArchive" $@
  typeset ExistFlag
  ExistFlag=$( [ -f $1 ] && echo 0 || echo 1 )
  IGCM_debug_PopStack "IGCM_sys_TestFileArchive"

  return ${ExistFlag}
}

#D-#==================================================
#D-function IGCM_sys_CountFileArchive
#D-* Purpose: Count files on Archive filesystem
#D-* Examples:
#D-
function IGCM_sys_CountFileArchive {
  IGCM_debug_PushStack "IGCM_sys_CountFileArchive" $@
  ls ${@} 2>/dev/null | wc -l
  if [ $? -gt 0 ] ; then
    echo "IGCM_sys_CountFileArchive : erreur."
  fi
  IGCM_debug_PopStack "IGCM_sys_CountFileArchive"
}

#D-#==================================================
#D-function IGCM_sys_Tree
#D-* Purpose: Tree directories with files on ${ARCHIVE}
#D-* Examples: IGCM_sys_Tree ${R_IN} ${R_OUT}
#D-
function IGCM_sys_Tree {
  IGCM_debug_PushStack "IGCM_sys_Tree" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_Tree :" $@
  fi

  \ls -lR ${@}

  IGCM_debug_PopStack "IGCM_sys_Tree"
}

#D-#==================================================
#D-function IGCM_sys_Qsub
#D-* Purpose: Qsub new job
#D-* Examples:
#D-
function IGCM_sys_Qsub {
  IGCM_debug_PushStack "IGCM_sys_Qsub" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_Qsub :" $@
  fi
  typeset options status
  options="-o ${SUBMIT_DIR}/${Script_Output} -e ${SUBMIT_DIR}/${Script_Output}"
  /usr/bin/time sbatch ${options} $1 > ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
  status=$?

  cat ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
  if [ ${status} -gt 0 ] ; then
    IGCM_debug_Print 2 "IGCM_sys_Qsub ${options} $1 : error code ${status}"
    IGCM_debug_Exit "IGCM_sys_Qsub"
  else
    JobID=$( gawk {'print $4'} ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ )
  fi
  IGCM_debug_PopStack "IGCM_sys_Qsub"
}

#D-#==================================================
#D-function IGCM_sys_QsubPost
#D-* Purpose: Qsub new job on scalaire
#D-* Examples:
#D-
function IGCM_sys_QsubPost {
  IGCM_debug_PushStack "IGCM_sys_QsubPost" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_QsubPost :" $@
  fi
  typeset options status
  options="-o ${POST_DIR}/${Script_Post_Output}.out  -e ${POST_DIR}/${Script_Post_Output}.out"
  sbatch ${options} ${libIGCM_POST}/$1.job > ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
  status=$?

  cat ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
  if [ ${status} -gt 0 ] ; then
    IGCM_debug_Print 2 "IGCM_sys_QsubPost ${options} ${libIGCM_POST}/$1.job : error code ${status}"
    IGCM_debug_Exit "IGCM_sys_QsubPost"
  else
    JobID=$( gawk {'print $4'} ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ )
  fi
  IGCM_debug_PopStack "IGCM_sys_QsubPost"
}

#D-*************************
#D- File transfer functions
#D-*************************
#D-

#D-#==================================================
#D-function IGCM_sys_RmRunDir
#D-* Purpose: rm tmpdir (dummy function most of the time batch
#D-                      scheduler will do the job)
#D-* Examples:
#D-
function IGCM_sys_RmRunDir {
  IGCM_debug_PushStack "IGCM_sys_RmRunDir" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_RmRunDir :" $@
  fi

  typeset status

  echo rm $@ > ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
  \rm $@ >> ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
  status=$?

  if [ ${status} -gt 0 ] ; then
    IGCM_debug_Print 1 "IGCM_sys_RmRunDir : rm error code is ${status}."
    cat ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
    IGCM_debug_Exit "IGCM_sys_RmRunDir"
  else
    \rm ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
  fi
  IGCM_debug_PopStack "IGCM_sys_RmRunDir"
}

#D-#==================================================
#D-function IGCM_sys_Put_Dir
#D-* Purpose: Copy a complete directory on $(ARCHIVE)
#D-* Examples:
#D-
function IGCM_sys_Put_Dir {
  IGCM_debug_PushStack "IGCM_sys_Put_Dir" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_Put_Dir :" $@
  fi
  if [ $DRYRUN = 0 ]; then
    if [ ! -d ${1} ] ; then
      echo "WARNING : IGCM_sys_Put_Dir ${1} DOES NOT EXIST ."
      IGCM_debug_PopStack "IGCM_sys_Put_Dir"
      return
    fi

    typeset status

    # Only if we use rsync
    #IGCM_sys_TestDirArchive $( dirname $2 )
    #
    #USUAL WAY
    \cp -r $1 $2 > ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
    status=$?

    if [ ${status} -gt 0 ] ; then
      IGCM_debug_Print 2 "IGCM_sys_Put_Dir : cp failed error code ${status}"
      cat ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
      IGCM_debug_Exit "IGCM_sys_Put_Dir"
    else
      \rm ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
    fi
  fi
  IGCM_debug_PopStack "IGCM_sys_Put_Dir"
}

#D-#==================================================
#D-function IGCM_sys_Get_Dir
#D-* Purpose: Copy a complete directory from ${ARCHIVE}
#D-* Examples:
#D-
function IGCM_sys_Get_Dir {
  IGCM_debug_PushStack "IGCM_sys_Get_Dir" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_Get_Dir :" $@
  fi
  if [ $DRYRUN = 0 ]; then
    typeset NB_ESSAI DELAI status i
    # number of tentative
    NB_ESSAI=3
    # time delay between tentative
    DELAI=2

    #
    # USUAL WAY
    i=0
    while [ $i -lt $NB_ESSAI ] ; do
      \cp -ur $1 $2 >> ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
      status=$?
      if [ ${status} -gt 0 ] ; then
        IGCM_debug_Print 2 "IGCM_sys_Get_Dir : cp failed error code ${status} ${i}/${NB_ESSAI}"
        IGCM_debug_Print 2 "IGCM_sys_Get_Dir : sleep ${DELAI} seconds and try again."
        sleep $DELAI
      else
        break
      fi
      (( i = i + 1 ))
    done

    if [ ${status} -gt 0 ] ; then
      IGCM_debug_Print 2 "IGCM_sys_Get_Dir : cp failed error code ${status}"
      cat ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
      IGCM_debug_Exit "IGCM_sys_Get_Dir"
    else
      \rm ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
    fi
  fi
  IGCM_debug_PopStack "IGCM_sys_Get_Dir"
}

#D-#==================================================
#D-function IGCM_sys_Put_Rest
#D-* Purpose: Put computied restarts on ${ARCHIVE}.
#D-           File and target directory must exist.
#D-* Examples:
#D-
function IGCM_sys_Put_Rest {
  IGCM_debug_PushStack "IGCM_sys_Put_Rest" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_Put_Rest :" $@
  fi
  if [ $DRYRUN = 0 ]; then
    if [ ! -f ${1} ] ; then
      echo "ERROR : IGCM_sys_Put_Rest ${1} DOES NOT EXIST ."
      IGCM_debug_Exit "IGCM_sys_Put_Rest"
    fi

    typeset status
    #
    # USUAL WAY
    \cp $1 $2 > ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
    status=$?

#       #RSYNC WITH NETWORK SSH CALL
#       echo ${RSYNC} ${RSYNC_opt} -e ssh ${RUN_DIR}/$1 ${STOREHOST}:${2} > ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
#       ${RSYNC} ${RSYNC_opt} -e ssh ${RUN_DIR}/$1 ${STOREHOST}:${2} >> ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1

#       #RSYNC WITH NFS USE
#       echo ${RSYNC} ${RSYNC_opt} ${RUN_DIR}/$1 ${2} > ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
#       ${RSYNC} ${RSYNC_opt} ${RUN_DIR}/$1 ${2} >> ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1

#       status=$?
#       IGCM_sys_Rsync_out $status

#       ${libIGCM}/libIGCM_sys/IGCM_analyse_rsync_out.awk ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
#       (( status=status+$? ))

    if [ ${status} -gt 0 ] ; then
      IGCM_debug_Print 2 "IGCM_sys_Put_Rest : cp failed error code ${status}"
      [ -f ${1} ] && ls -l ${1}
      [ -f ${2} ] && ls -l ${2}
      [ -f ${2}/${1} ] && ls -l ${2}/${1}
      cat ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
      IGCM_debug_Exit "IGCM_sys_Put_Rest"
    else

      if [ X${JobType} = XRUN ] ; then
        [ -f ${2} ] && IGCM_sys_Chmod 444 ${2}
        [ -f ${2}/${1} ] && IGCM_sys_Chmod 444 ${2}/${1}
      fi

      \rm ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
    fi
  fi
  IGCM_debug_PopStack "IGCM_sys_Put_Rest"
}

#D-#==================================================
#D-function IGCM_sys_Put_Out
#D-* Purpose: Copy a file on ${ARCHIVE} after having chmod it in readonly
#D-* Examples:
#D-
function IGCM_sys_Put_Out {
  IGCM_debug_PushStack "IGCM_sys_Put_Out" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_Put_Out :" $@
  fi

  typeset status

  if [ $DRYRUN = 0 ]; then
    if [ ! -f ${1} ] ; then
      echo "WARNING : IGCM_sys_Put_Out ${1} DOES NOT EXIST ."
      IGCM_debug_PopStack "IGCM_sys_Put_Out"
      return 1
    fi
    #
    IGCM_sys_MkdirArchive $( dirname $2 )
    #
    if [ X${JobType} = XRUN ] ; then
      if [ X${3} = X ] ; then
        IGCM_sys_Chmod 444 ${1}
      fi
    fi
    #

    echo ${RSYNC} ${RSYNC_opt} $1 $2 > ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
    ${RSYNC} ${RSYNC_opt} $1 $2 >> ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
    status=$?
    IGCM_sys_Rsync_out $status

    ${libIGCM}/libIGCM_sys/IGCM_analyse_rsync_out.awk ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
    (( status=status+$? ))

    if [ ${status} -gt 0 ] ; then
      IGCM_debug_Print 2 "IGCM_sys_Put_Out : rsync failed error code ${status}"
      cat ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
      IGCM_debug_Exit "IGCM_sys_Put_Out"
    else
      \rm ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
    fi
  fi
  IGCM_debug_PopStack "IGCM_sys_Put_Out"
  return 0
}

#D-#==================================================
#D-function IGCM_sys_Get
#D-* Purpose: Get a file from ${ARCHIVE}
#D-* Examples: IGCM_sys_Get myfile /destpath/myfile_with_PREFIX
#D-            IGCM_sys_Get /l Array_contain_myfiles /destpath/
function IGCM_sys_Get {
  IGCM_debug_PushStack "IGCM_sys_Get" $@

  typeset DEST status dm_liste

  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_Get :" $@
  fi
  if [ $DRYRUN -le 2 ]; then
    if [ X${1} = X'/l' ] ; then
      eval set +A dm_liste \${${2}}
    else
      dm_liste=${1}
    fi
    eval DEST=\${${#}}

    # test if the (first) file is present in the old computation :
    IGCM_sys_TestFileArchive ${dm_liste[0]}
    status=$?
    if [ ${status} -gt 0 ] ; then
      echo "IGCM_sys_Get, ERROR : regular file ${dm_liste[0]} DOES NOT EXIST ."
      IGCM_debug_Exit "IGCM_sys_Get"
      return
    fi

    #USUAL WAY
    \cp ${dm_liste[*]} ${DEST} > ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
    status=$?

    if [ ${status} -gt 0 ] ; then
      IGCM_debug_Print 2 "IGCM_sys_Get : cp failed error code ${status}"
      cat ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
      IGCM_debug_Exit "IGCM_sys_Get"
    else
      \rm ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
    fi
  fi
  IGCM_debug_PopStack "IGCM_sys_Get"
}

#D-#==================================================
#D-function IGCM_sys_GetDate_Monitoring
#D-* Purpose: get the last year for which the monitoring has been computed
#D-* Examples:
#D-
function IGCM_sys_GetDate_Monitoring {
  IGCM_debug_PushStack "IGCM_sys_GetDate_Monitoring" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_GetDate_Monitoring :" $@
  fi

  eval ${2}=$( cdo showyear ${1} 2> /dev/null | gawk '{ print $NF }' )

  IGCM_debug_PopStack "IGCM_sys_GetDate_Monitoring"
}

#D-#==================================================
#D-function IGCM_sys_Dods_Rm
#D-* Purpose: DO NOTHING ! Put ${ARCHIVE} files on DODS internet protocole.
#D-* Examples:
#D-
function IGCM_sys_Dods_Rm {
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_Dods_Rm :" $@
  fi
  return 0
}

#D-#==================================================
#D-function IGCM_sys_Dods_Cp
#D-* Purpose: Copy $(ARCHIVE) files on DODS internet protocole.
#D-* Examples:
#D-
function IGCM_sys_Dods_Cp {
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_Dods_Cp :" $@
  fi
  return 0
}

#D-#==================================================
#D-function IGCM_sys_Put_Dods
#D-* Purpose: Put ${ARCHIVE} files on DODS internet protocole. Dummy function here
#D-* Examples:
#D-
function IGCM_sys_Put_Dods {
  IGCM_debug_PushStack "IGCM_sys_Put_Dods" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_Put_Dods :" $@
  fi
  IGCM_debug_PopStack "IGCM_sys_Put_Dods"
}

##############################################################
# REBUILD OPERATOR

#D-#==================================================
#D-function IGCM_sys_sync
#D-* Purpose: flush buffer on disk (dummy function on Ciclad)
#D-* Examples:
#D-
function IGCM_sys_sync {
  IGCM_debug_PushStack "IGCM_sys_sync" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_sync :" $@
    echo "Dummy call, let the system do that."
  fi
  IGCM_debug_PopStack "IGCM_sys_sync"
}

############################################################
# Activate Running Environnment Variables

#D-#==================================================
#D-function IGCM_sys_activ_variables
#D-* Purpose: set environement variables prior to execution
#D-* Examples:
#D-
function IGCM_sys_activ_variables {
  IGCM_debug_PushStack "IGCM_sys_activ_variables"
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_activ_variables"
  fi

# --------------------------------------------------------------------
#D- MODULE specifications
# --------------------------------------------------------------------
  # Source the file EnvFile containing all module load needed to run the model.
  # EnvFile can be specified in config.card. If this is the case and if the file exists, 
  # it will be used. Otherwise a default file will be used.
  if [ ! X${config_UserChoices_EnvFile} = X ] && [ -f ${config_UserChoices_EnvFile} ] ; then
      # EnvFile is set in config.card and the file exists
      IGCM_debug_Print 1 "EnvFile set in config.card will be used"
      EnvFile=${config_UserChoices_EnvFile}
  else
      IGCM_debug_Print 1 "IGCM_sys_active_variables : Default modules will be used"
      if [ $CENTER == spirit ] ; then
          # At spirit
	  EnvFile=/home/igcmg/MachineEnvironment/mesoipsl/atlas_env_mesoipsl
      else
          # At spirit
	  EnvFile=/ciclad-home/igcmg/MachineEnvironment/mesoipsl/atlas_env_mesoipsl
      fi
  fi
  IGCM_debug_Print 1 "IGCM_sys_active_variables : Following EnvFile will be used :" 
  IGCM_debug_Print 1 ${EnvFile}
  . ${EnvFile}


  IGCM_debug_Print 1 "IGCM_sys_active_variables : Now loaded modules for meso-ipsl . "
  module list
  

# --------------------------------------------------------------------
#D- MPI specifications
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#D- Other specifications
# --------------------------------------------------------------------

  ulimit -s unlimited

  IGCM_debug_PopStack "IGCM_sys_activ_variables"
}

############################################################
# Desactivate Running Environnment Variables

#D-#==================================================
#D-function IGCM_sys_desactiv_variables
#D-* Purpose: unset environement variables after execution
#D-* Examples:
#D-
function IGCM_sys_desactiv_variables {
  IGCM_debug_PushStack "IGCM_sys_desactiv_variables"
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_desactiv_variables"
  fi
# --------------------------------------------------------------------
#D- MPI specifications
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#D- Other specifications
# --------------------------------------------------------------------

  IGCM_debug_PopStack "IGCM_sys_desactiv_variables"
}

############################################################
# Update job headers to be used by the scheduler

#D-#==================================================
#D-function IGCM_sys_updateHeaders
#D-* Purpose: Update job headers to be used by the scheduler
#D-* Examples: IGCM_sys_updateHeaders /path/to/Job_MYEXP
#D-
function IGCM_sys_updateHeaders {
  IGCM_debug_PushStack "IGCM_sys_updateHeaders"
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_updateHeaders"
  fi
  typeset file
  file=$1

  if [ ${executionType} -eq 1 ] ; then
    # MPMD + MPI
    sed -e "/::openMPthreads::/d"                  \
        -e "s/::JobNumProcTot::/${coreNumber}/"    \
      ${file} > ${file}.tmp

  elif [ ${executionType} -eq 2 ] ; then
    # MPMD + MPI + OMP
      (( nodeNumber = coreNumber / NB_CORE_PER_NODE ))
      [ $(( ${coreNumber} % ${NB_CORE_PER_NODE} )) -ne 0 ] && (( nodeNumber = nodeNumber + 1 ))
      sed -e "/::openMPthreads::/d"                 \
        -e "s/::JobNumProcTot::/${mpiTasks}/"       \
	-e "/ntasks/i\#SBATCH --nodes=${nodeNumber}"\
        -e "/ntasks/i\#SBATCH --exclusive"          \
         ${file} > ${file}.tmp
  elif [ ${executionType} -eq 3 ] ; then
    # SPMD + MPI/OMP
    sed -e "s/::openMPthreads::/${openMPthreads}/" \
        -e "s/::JobNumProcTot::/${mpiTasks}/"      \
      ${file} > ${file}.tmp

  elif [ ${executionType} -eq 4 ] ; then
    # SPMD + MPI only
    sed -e "s/::JobNumProcTot::/${mpiTasks}/"      \
        -e "/::openMPthreads::/d"                  \
      ${file} > ${file}.tmp

  elif [ ${executionType} -eq 5 ] ; then
    # SPMD + OMP only
    sed -e "s/::openMPthreads::/${openMPthreads}/" \
        -e "/::JobNumProcTot::/d"                  \
      ${file} > ${file}.tmp

  elif [ ${executionType} -eq 6 ] ; then
    # SEQUENTIAL THEN
    sed -e "s/::JobNumProcTot::/1/"                \
        -e "/::openMPthreads::/d"                  \
      ${file} > ${file}.tmp

  fi

  IGCM_sys_Mv ${file}.tmp ${file}

  IGCM_debug_PopStack "IGCM_sys_updateHeaders"
}

############################################################
# Build MPI/OMP scripts run file (dummy function)

#D-#==================================================
#D-function IGCM_sys_build_run_file
#D-* Purpose: build run file (deprecated)
#D-* Examples:
#D-
function IGCM_sys_build_run_file {
  IGCM_debug_Print 3 " dummy function : IGCM_sys_build_run_file "
}

############################################################
# Build MPI/OMP scripts

#D-#==================================================
#D-function IGCM_sys_build_execution_scripts
#D-* Purpose: build execution scripts to be launch by ${HOST_MPIRUN_COMMAND}
#D-* Examples:
#D-
function IGCM_sys_build_execution_scripts
{
  IGCM_debug_PushStack "IGCM_sys_build_execution_scripts" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_build_execution_scripts " $@
  fi
  IGCM_debug_Print 1 "executionType= ${executionType}"
  EXECUTION=${HOST_MPIRUN_COMMAND}

  if ( ${OK_PARA_MPMD} ) ; then

    # MPMD mode
    # 1 MPI only : executionType=1
    # 2 MPI/OpenMP : executionType=2

    if [ -f run_file ] ; then
      IGCM_sys_Rm -f run_file
    fi
    touch run_file

    # case 1 : Only MPI (MPMD)
    if  ( ! ${OK_PARA_OMP} ) ; then
      # Build run_file
      current_core=0
      # First loop on the components for the coupler ie oasis (only if oasis3)
      # the coupler ie oasis3 must be the first one
      for comp in ${config_ListOfComponents[*]} ; do

        eval ExeNameIn=\${config_Executable_${comp}[0]}
        eval ExeNameOut=\${config_Executable_${comp}[1]}

        # Only if the component has an executable
        if ( [ "X${ExeNameOut}" != X\"\" ] ) ; then

#          eval comp_proc_mpi_loc=\${${comp}_PROC_MPI}
#          (( end_core = ${current_core} + ${comp_proc_mpi_loc} - 1 ))
#          echo "${current_core}-${end_core} ./${ExeNameOut}" >> run_file
#          (( current_core = ${end_core} + 1 ))
          eval comp_proc_mpi_loc=\${${comp}_PROC_MPI}
          eval comp_proc_omp_loc=\${${comp}_PROC_OMP}
          echo "-np ${comp_proc_mpi_loc} ./${ExeNameOut}" >> run_file

	fi
      done

      EXECUTION="${HOST_MPIRUN_COMMAND} --app ./run_file"

      IGCM_sys_Chmod u+x run_file
      if ( $DEBUG_sys ) ; then
        echo "run_file contains : "
        cat run_file
      fi

    else 

    # 2 MPI/OpenMP : executionType=2

    # MPI-OpenMP (MPMD)
    # export SLURM_HOSTFILE=./hostlist 
    # srun --cpu-bind=none --distribution=arbitrary --multi-prog ./run_file
    # example of  hostlist file :
    # r3i3n33
    # r3i3n33
    # ...
    # example of run_file :
    # 0-70 ./prog_lmdz.x.sh %o %t
    # 71-430 ./prog_opa.xx.sh %o %t
    # 431-431 ./prog_xios.x.sh %o %t
    # examples of prog_file :
    # prog_lmdz.x.sh :
    # (( init = 0 + $1 ))
    # (( index = init * 10 ))
    # (( slot = index % 40 ))
    # taskset -c $slot-$((slot + 10 - 1)) ./script_lmdz.x.ksh
    # that will become
    # taskset -c 0-9 ./script_lmdz.x.ksh
    # ...
    # with script_lmdz.x.ksh
    # export OMP_STACKSIZE=3g
    # export OMP_PLACES=cores
    # export OMP_NUM_THREADS=10
    # ./lmdz.x > out_lmdz.x.out.${SLURM_PROCID} 2>out_lmdz.x.err.${SLURM_PROCID}

    #  Hosts treatment
    _bkIFS=$IFS;
	IFS=$'\n'; set -f
	listnodes=($(< <( scontrol show hostnames $SLURM_JOB_NODELIST )))
	IFS=$_bkIFS; set +f
	rm -f hostlist

     # Loop on the components to build run_file and script_exec files
	rank=0
	current_core=0
        current_core_mpi=0

      for comp in ${config_ListOfComponents[*]} ; do

        eval ExeNameIn=\${config_Executable_${comp}[0]}
        eval ExeNameOut=\${config_Executable_${comp}[1]}

        # Not possible if oasis has an executable (i.e old version of oasis3)
        if ( [ "X${ExeNameOut}" != X\"\" ] &&  [ "X${comp}" = "XCPL" ] ) ; then
          IGCM_debug_Exit "ERROR MPMD with hybrid MPI-OpenMP is not available with oasis3 version"
          IGCM_debug_Print 2 "Only available with oasis3-MCT version coupler"
          IGCM_debug_Verif_Exit
        fi

        # Only if we really have an executable for the component :
        if [ "X${ExeNameOut}" != X\"\" ] ; then

          eval comp_proc_mpi_loc=\${${comp}_PROC_MPI}
          eval comp_proc_omp_loc=\${${comp}_PROC_OMP}
          eval comp_proc_nod_loc=\${${comp}_PROC_NOD}
     

          # Build script files

          echo "#!/bin/ksh" > script_${ExeNameOut}.ksh
          echo ""  >> script_${ExeNameOut}.ksh
          if [ ${comp_proc_omp_loc} -gt 1 ] ; then
          # Check if the number of threads is correct
	      case ${comp_proc_omp_loc} in
		  2|4|5|10|20)
		      IGCM_debug_Print 1 "You run ${ExeNameOut} on ${comp_proc_omp_loc} OMP threads"
		      ;;
		  *)
		      IGCM_debug_Exit "ERROR with OMP parameters !"
		      IGCM_debug_Print 2 "${comp_proc_omp_loc} is not possible as number of OMP threads"
		      IGCM_debug_Print 2 "Only 2,4,5,10,20 as number of OMP threads are possible "
		      IGCM_debug_Verif_Exit
		      ;;
	      esac
            echo "export OMP_STACKSIZE=3g"  >> script_${ExeNameOut}.ksh
            echo "export OMP_PLACES=cores"  >> script_${ExeNameOut}.ksh
            echo "export OMP_NUM_THREADS=${comp_proc_omp_loc}"  >> script_${ExeNameOut}.ksh

          fi

          # to have out/err per process on different files
          echo "./${ExeNameOut} > out_${ExeNameOut}.out.\${SLURM_PROCID} 2>out_${ExeNameOut}.err.\${SLURM_PROCID}" >> script_${ExeNameOut}.ksh

          IGCM_sys_Chmod u+x script_${ExeNameOut}.ksh

        # Build run_file
        # Only if the component has an executable
        if ( [ "X${ExeNameOut}" != X\"\" ] ) ; then

          eval comp_proc_mpi_loc=\${${comp}_PROC_MPI}
          (( end_core = ${current_core_mpi} + ${comp_proc_mpi_loc} - 1 ))
          echo "${current_core_mpi}-${end_core} ./prog_${ExeNameOut}.sh %o %t" >> run_file
          (( current_core_mpi = ${end_core} + 1 ))
        fi

          if [ ${comp_proc_nod_loc} -gt 1 ] ; then
              (( offset_comp_proc_loc =  NB_CORE_PER_NODE / (comp_proc_mpi_loc / comp_proc_nod_loc) ))
          else
              (( offset_comp_proc_loc =  comp_proc_omp_loc ))
          fi

          # Build configuration file
 
	      echo "#!/bin/sh" > prog_${ExeNameOut}.sh
          echo "(( init = $current_core + \$1 ))" >> prog_${ExeNameOut}.sh
          echo "(( index = init * $comp_proc_omp_loc ))" >> prog_${ExeNameOut}.sh
          echo "(( slot = index % 40 ))" >> prog_${ExeNameOut}.sh
	      echo "echo ${ExeNameOut} taskset -c \$slot"-"\$((slot + $comp_proc_omp_loc - 1))" >> prog_${ExeNameOut}.sh
          echo "taskset -c \$slot"-"\$((slot + $comp_proc_omp_loc - 1)) ./script_${ExeNameOut}.ksh" >> prog_${ExeNameOut}.sh

          IGCM_sys_Chmod u+x prog_${ExeNameOut}.sh

        # Build hostlist file

          for nb_proc_mpi in `seq 0 $(($comp_proc_mpi_loc-1))`; do
             (( index_host = current_core / NB_CORE_PER_NODE ))
             host_value=${listnodes[${index_host}]}
               echo "$host_value" >> hostlist
             (( current_core = current_core + offset_comp_proc_loc ))
          done
        fi
      done

      ## variable added to stop after 60s instead of 600s by default. 
      ## This is used when no error comes from executables and when something stopped an executable without notice. 
      export SLURM_WAIT=60

      EXECUTION="${HOST_MPIRUN_COMMAND} --cpu-bind=none --distribution=arbitrary --multi-prog ./run_file"

      IGCM_sys_Chmod u+x run_file
      if ( $DEBUG_sys ) ; then
        echo "run_file contains : "
        cat run_file
      fi

    fi # if ${OK_PARA_MPMD}

  else
  # Only one executable (SPMD mode):  executionType=3, 4, 5 and 6

    for comp in ${config_ListOfComponents[*]} ; do

      # Only if we really have an executable for the component :
      eval ExeNameOut=\${config_Executable_${comp}[1]}
      if ( [ "X${ExeNameOut}" != X\"\" ] && [ "X${ExeNameOut}" != "Xinca.dat" ] ) ; then

        # Build script files

        echo "#!/bin/ksh" > script_${ExeNameOut}.ksh
        echo ""  >> script_${ExeNameOut}.ksh
        IGCM_sys_Chmod u+x script_${ExeNameOut}.ksh

        if ( ${OK_PARA_OMP} ) ; then
            eval comp_proc_omp_loc=\${${comp}_PROC_OMP}
	    # Check if the number of threads is correct
	    case ${comp_proc_omp_loc} in
		2|4|5|10|20)
		    IGCM_debug_Print 1 "You run ${ExeNameOut} on ${comp_proc_omp_loc} OMP threads"
		    ;;
		*)
		    IGCM_debug_Exit "ERROR with OMP parameters !"
		    IGCM_debug_Print 2 "${comp_proc_omp_loc} is not possible as number of OMP threads"
		    IGCM_debug_Print 2 "Only 2,4,5,10,20 as number of OMP threads are possible "
		    IGCM_debug_Verif_Exit
		    ;;
	    esac 
            echo ""  >> script_${ExeNameOut}.ksh
            echo "export OMP_STACKSIZE=3g"  >> script_${ExeNameOut}.ksh
            echo "export OMP_PLACES=cores"  >> script_${ExeNameOut}.ksh
            echo "OMP_NUM_THREADS=${comp_proc_omp_loc}" >> script_${ExeNameOut}.ksh
        fi

        eval comp_proc_mpi_loc=\${${comp}_PROC_MPI}
        
        # To have out/err per process on different files
        echo "./${ExeNameOut} > out_${ExeNameOut}.out.\${SLURM_PROCID} 2>out_${ExeNameOut}.err.\${SLURM_PROCID}"  >> script_${ExeNameOut}.ksh
        EXECUTION="${HOST_MPIRUN_COMMAND} ./script_${ExeNameOut}.ksh"

        IGCM_debug_Print 1 "sys Jean-Zay : script_${ExeNameOut}.ksh contains"
        cat script_${ExeNameOut}.ksh

      fi

    done

  fi # ${OK_PARA_MPMD}

  IGCM_debug_Print 1 "sys meso-ipsl : execution command is "
  IGCM_debug_Print 1 "$EXECUTION"

  IGCM_debug_PopStack "IGCM_sys_build_execution_scripts"
}

#D-#==================================================
#D-function IGCM_sys_check_path
#D-* Purpose: check that RUN_DIR_PATH that will be removed on some machine
#D-* do not point to an important use directory. Stop immediately in that case.
#D-* Examples:
#D-
function IGCM_sys_check_path {
  IGCM_debug_PushStack "IGCM_sys_check_path"
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_check_path"
  fi

  if ( [ X${RUN_DIR_PATH} = X${HOME} ] || [ X${RUN_DIR_PATH} = X${ARCHIVE} ] ) ; then
    IGCM_debug_Print 1 "Variable RUN_DIR_PATH is pointing to an important directory : ${RUN_DIR_PATH}"
    IGCM_debug_Print 1 "Please check the RUN_DIR_PATH definition in your Job : Job_${config_UserChoices_JobName}"
    IGCM_debug_Exit "This will stop the job"
  fi
  IGCM_debug_PopStack "IGCM_sys_check_path"
}

#D-#==================================================
#D-function IGCM_sys_check_quota. Dummy call here
#D-* Purpose: check user quota. Stop the simulation if quota above 90%
#D-* Examples:
#D-
function IGCM_sys_check_quota {
  IGCM_debug_PushStack "IGCM_sys_check_quota"
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_check_quota"
  fi

  # TO BE IMPLEMENTED

  IGCM_debug_PopStack "IGCM_sys_check_quota"
}

#D-#==================================================
#D-function IGCM_sys_projectAccounting
#D-* Purpose: store project accounting information in a file
#D-* Examples:
#D-
function IGCM_sys_projectAccounting {
  IGCM_debug_PushStack "IGCM_sys_projectAccounting"
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_check_quota"
  fi

  touch $1

  IGCM_debug_PopStack "IGCM_sys_projectAccounting"
}

#D-#==================================================
#D-function IGCM_sys_getJobSchedulerID
#D-* Purpose: Get the job ID during execution
#D-* Examples: IGCM_sys_getJobSchedulerID jobSchedulerID
#D-
function IGCM_sys_getJobSchedulerID {
  IGCM_debug_PushStack "IGCM_sys_getJobSchedulerID"
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_getJobSchedulerID"
  fi

  eval ${1}=$( echo ${PBS_JOBID} | awk -F. '{print $1}' )

  IGCM_debug_PopStack "IGCM_sys_getJobSchedulerID"
}

#D-#==================================================
#D-function IGCM_sys_GetJobID
#D-* Purpose: Get the job ID from the JobName
#D-* Examples: IGCM_sys_GetJobID ${JobName} ${TargetUsr} JobID
#D-
function IGCM_sys_GetJobID {
  IGCM_debug_PushStack "IGCM_sys_GetJobID"
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_GetJobID"
  fi

  # With -f option, the full job name is given in the last column
  ID="$( qstat -u $2 | grep -w $1 | gawk '-F ' '{print $10}' )"

  eval ${3}=${ID}
  IGCM_debug_PopStack "IGCM_sys_GetJobID"
}

#D-#==================================================
#D-function IGCM_sys_CountJobInQueue
#D-* Purpose: Count number of users job
#D-* Examples: IGCM_sys_CountJobInQueue ${JobName} NbRun
#D-
function IGCM_sys_CountJobInQueue {
  IGCM_debug_PushStack "IGCM_sys_CountJobInQueue"
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_CountJobInQueue"
  fi
  IGCM_debug_PopStack "IGCM_sys_CountJobInQueue"
}

#D-#==================================================
#D-function IGCM_sys_atlas
#D-* Purpose: encapsulate atlas call so as to manage error code and curie specificity
#D-* Examples:
#D-
function IGCM_sys_atlas {
  IGCM_debug_PushStack "IGCM_sys_atlas" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_atlas :" $@
  fi

  typeset status

  \atlas $@ > ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$ 2>&1
  status=$?
  if [ ${status} -gt 0 ] ; then
    echo "IGCM_sys_atlas : error code ${status}"
    cat ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
    \rm ${OUTCOMMAND_PATH}/out_command_${LOGIN}.$$
    IGCM_debug_PopStack "IGCM_sys_atlas"
    return 1
  else
    IGCM_debug_PopStack "IGCM_sys_atlas"
    return 0
  fi

  IGCM_debug_PopStack "IGCM_sys_atlas"
}

#D-#==================================================
#D-function IGCM_sys_rebuild_nemo
#D-* Purpose: rebuild nemo parallel files with general rebuild
#D-* Examples:
#D-

function IGCM_sys_rebuild_nemo {
  IGCM_debug_PushStack "IGCM_sys_rebuild_nemo" $@
  if ( $DEBUG_sys ) ; then
    echo "IGCM_sys_rebuild_nemo :" $@
  fi

  nemo_generic_restart_file_name_out=${1}
  nemo_extension_out=${3}
  shift ; shift ; shift

  IGCM_sys_rebuild -o ${nemo_generic_restart_file_name_out}.${nemo_extension_out} $@

  IGCM_debug_PopStack "IGCM_sys_rebuild_nemo"
}



