#!/bin/ksh

#D--------------------------------------------------------------------==

function ENSEMBLE_initialize
{
    IGCM_debug_PushStack "ENSEMBLE_initialize"

    # config.card : 
    # -------------
    IGCM_card_DefineArrayFromSection ${SUBMIT_DIR}/config.card UserChoices
    typeset option
    for option in ${config_UserChoices[*]} ; do
	IGCM_card_DefineVariableFromOption ${SUBMIT_DIR}/config.card UserChoices ${option}
    done
    
    DateBegin=$( IGCM_date_ConvertFormatToGregorian ${config_UserChoices_DateBegin} )
    IGCM_date_GetYearMonth ${DateBegin} year month
    
    IGCM_card_DefineArrayFromSection ${SUBMIT_DIR}/config.card SubJobPost
    typeset option
    for option in ${config_SubJobPost[*]} ; do
	IGCM_card_DefineVariableFromOption ${SUBMIT_DIR}/config.card SubJobPost ${option}
    done
    
    # Subjobs configuration : 
    # -----------------------

    # 1) UserChoices section of your card

    # This section groups standard options to be passed to ORCHIDEE.
    IGCM_card_DefineArrayFromSection ${SUBMIT_DIR}/ensemble.card UserChoices
    set +A ensemble_UserChoices -- ${ensemble_UserChoices[*]} > /dev/null 2>&1
    typeset option
    for option in ${ensemble_UserChoices[*]} ; do
	IGCM_card_DefineVariableFromOption ${SUBMIT_DIR}/ensemble.card UserChoices ${option}
    done

    IGCM_card_DefineArrayFromSection ${SUBMIT_DIR}/config.card SRF
    typeset option
    for option in ${config_SRF[*]} ; do
	IGCM_card_DefineVariableFromOption ${SUBMIT_DIR}/config.card SRF ${option}
    done
    if [ X"${ensemble_UserChoices_ok_stomate}" = "Xy" ] ; then
	IGCM_card_DefineArrayFromSection ${SUBMIT_DIR}/config.card SBG
	typeset option
	for option in ${config_SBG[*]} ; do
	    IGCM_card_DefineVariableFromOption ${SUBMIT_DIR}/config.card SBG ${option}
	done
    fi
    
    # 2) CONFIG section :

    # This section groups specific options for ENSEMBLE configuration :
    # - Forcing files path
    IGCM_card_DefineVariableFromOption ${SUBMIT_DIR}/ensemble.card CONFIG ForcingPath

    # - List of the sites/information to be treated :
    #   * Number of PFTs
    IGCM_card_DefineVariableFromOption ${SUBMIT_DIR}/ensemble.card CONFIG NbPFTs
    NbPFTs=${ensemble_CONFIG_NbPFTs}

    #   * Number of physical parameters on each sites per PFTs
    IGCM_card_DefineVariableFromOption ${SUBMIT_DIR}/ensemble.card CONFIG NbSitesParam
    #   * ORCHIDEE name for physical parameters on each sites
    IGCM_card_DefineArrayFromOption ${SUBMIT_DIR}/ensemble.card CONFIG NameSitesParam
    #   * Sites descriptions

    IGCM_card_DefineArrayFromOption ${SUBMIT_DIR}/ensemble.card CONFIG Groups

    if [ X${ensemble_CONFIG_Groups[0]} == X"Option" ] ; then
	typeset name=ensemble_CONFIG_Groups
	eval unset ${name}
	eval ${name}[0]=${NULL_STR}
	set +A ${name} -- Sites
    fi
    for group in ${ensemble_CONFIG_Groups[*]} ; do
	IGCM_card_DefineArrayFromOption ${SUBMIT_DIR}/ensemble.card CONFIG $group
    done

    first_group=${ensemble_CONFIG_Groups[*]}
    echo 'first group:' ${first_group}

    NumInfosBySite=4
    (( iphys = 0 ))
    while [ $iphys -lt ${ensemble_CONFIG_NbSitesParam} ] ; do
	
	# Detect if physical parameter has 
        #  array description or line per PFT description in the parameter file
	SearchParam=""
	ORCHIDEE_grep ${ensemble_CONFIG_NameSitesParam[${iphys}]}
	echo 'Param:' ${ensemble_CONFIG_NameSitesParam[${iphys}]}
	echo 'SearchParam:' ${SearchParam}
	echo 'nbPFTs:' $NbPFTs
	case ${SearchParam} in
	    undefined)
		echo "ERROR : parameter ${ensemble_CONFIG_NameSitesParam[${iphys}]} doesn't exist in run.def parameter file !"
		echo "You must correct NameSiteParam option in your card file or add new parameter in run.def file"
		echo "We must STOP here."
		exit 1
		;;
	    value)
		(( NumInfosBySite = NumInfosBySite + 1 ))
		;;
	    line)
		(( NumInfosBySite = NumInfosBySite + NbPFTs ))
		;;
	    vector)
		(( NumInfosBySite = NumInfosBySite + NbPFTs ))
		;;
	esac

	(( iphys = iphys + 1 ))
    done

    echo 'NumInfosBySite: ' $NumInfosBySite


    # Num of Columns : name, filename, length, + ${NbPFTs} PFT
    ColumnName=0
    ColumnFilename=1
    ColumnYear=2
    ColumnLength=3



    # 3) SPINUP section : files path and spinup parameters
    # This section groups specific options for SPINUP :
    IGCM_card_DefineArrayFromSection ${SUBMIT_DIR}/ensemble.card SPINUP
    set +A ensemble_SPINUP -- ${ensemble_SPINUP[*]} > /dev/null 2>&1
    typeset option
    for option in ${ensemble_SPINUP[*]} ; do
	IGCM_card_DefineVariableFromOption ${SUBMIT_DIR}/ensemble.card SPINUP ${option}
    done

    # 4) Subjob run.def parametrisation : SubJobParams in your card
    # This section give flags to put in run.def of Subjobs.
    IGCM_card_DefineArrayFromSection ${SUBMIT_DIR}/ensemble.card SubJobParams
    set +A ensemble_SubJobParams -- ${ensemble_SubJobParams[*]} > /dev/null 2>&1
    eval first_option=${ensemble_SubJobParams[0]} > /dev/null 2>&1
    if [ X${first_option} != X"Error:" ] ; then
	typeset option
	for option in ${ensemble_SubJobParams[*]} ; do
	    IGCM_card_DefineVariableFromOption ${SUBMIT_DIR}/ensemble.card SubJobParams ${option}
	done
    fi

    IGCM_debug_PopStack "ENSEMBLE_initialize"
}

# Function to sed parameter files in Subjob dirs
# For run.def :
function ORCHIDEE_def
{
    IGCM_debug_PushStack "ORCHIDEE_def"
    
    typeset ORCHIDEE_def_KEY ORCHIDEE_default_KEY

    ORCHIDEE_def_KEY=$( grep -e "^${1}\ *=" ${New_SUBMIT_DIR}/PARAM/run.def | gawk -- 'BEGIN {FS="="} {print $2}')
    if [ -f ${New_SUBMIT_DIR}/PARAM/orchidee.default ] ; then
	ORCHIDEE_default_KEY=$( grep -e "^${1}\ *=" ${New_SUBMIT_DIR}/PARAM/orchidee.default | gawk -- 'BEGIN {FS="="} {print $2}')
    fi

    if [ X"${ORCHIDEE_def_KEY}" != X ] ; then
	IGCM_debug_Print 2 "ORCHIDEE : ${1} has already been set in def file."
	if [ -f ${New_SUBMIT_DIR}/PARAM/orchidee.default ] ; then
	    IGCM_debug_Print 2 " default value : ${ORCHIDEE_default_KEY}"
	fi
	IGCM_debug_Print 2 "ORCHIDEE : ${1} has already been set in run.def file."
	IGCM_debug_Print 2 " default value : ${ORCHIDEE_default_KEY}"
	IGCM_debug_Print 2 " script value : ${2}"
	IGCM_debug_Print 2 " USER value : ${ORCHIDEE_def_KEY}"
	IGCM_debug_Print 2 " We will NOT set in again !"
    else
	echo "${1}= ${2}" >> ${New_SUBMIT_DIR}/PARAM/run.def
	echo "" >> ${New_SUBMIT_DIR}/PARAM/run.def
	echo "ORCHIDEE_def : ${1} ${2}"
    fi

    IGCM_debug_PopStack "ORCHIDEE_def"
    return $RET
}

# Function to grep parameter files in Subjob dirs
# They search for the argmument value in the file.
# The result of the function is given in ${SearchParam} global variable :
# * "line" if the parameter is a vector defined line per line with the form 
#   PARAM__1=value1
#   PARAM__2=value2
#   etc...
# * "vector" if the parameter is defined with the form
#   PARAM=value1[,] value2[,] etc...
# * "value" if the parameter is a single value
#   PARAM=value
# * "undefined" if the parameter is not definied in the parameters file
#   This should give an error in the script.

# For run.def :
function ORCHIDEE_grep
{
    IGCM_debug_PushStack "ORCHIDEE_grep"
    typeset Output

    if [ -f ${SUBMIT_DIR}/PARAM/run.def ] ; then
	set +A Output -- $( grep "^${1}__[0-9]*\ *=" ${SUBMIT_DIR}/PARAM/run.def )
	if [ X${Output[0]} != X ] ; then
	    SearchParam=line
	else
	    unset Output
	    set +A Output -- $( grep "^${1}\ *=\ *[^[:blank:]][^[:blank:]]*\ [^[:blank:]][^[:blank:]]*.*" ${SUBMIT_DIR}/PARAM/run.def )
	    if [ X${Output[0]} != X ] ; then
		SearchParam=vector
	    else
		unset Output
		set +A Output -- $( grep "^${1}\ *=" ${SUBMIT_DIR}/PARAM/run.def )
		if [ X${Output[0]} != X ] ; then
		    SearchParam=value
		else
		    unset Output
                    if [ -f ${SUBMIT_DIR}/PARAM/orchidee.default ] ; then
            		set +A Output -- $( grep "^${1}__[0-9]*\ *=" ${SUBMIT_DIR}/PARAM/orchidee.default )
            		if [ X${Output[0]} != X ] ; then
            		    SearchParam=line
            		else
            		    unset Output
            		    set +A Output -- $( grep "^${1}\ *=\ *[^[:blank:]][^[:blank:]]*\ [^[:blank:]][^[:blank:]]*.*" ${SUBMIT_DIR}/PARAM/orchidee.default )
            		    if [ X${Output[0]} != X ] ; then
            			SearchParam=vector
            		    else
            			unset Output
            			set +A Output -- $( grep "^${1}\ *=" ${SUBMIT_DIR}/PARAM/orchidee.default )
            			if [ X${Output[0]} != X ] ; then
            			    SearchParam=value
            			else
				    SearchParam=undefined
				fi
			    fi
			fi
		    else
			SearchParam=undefined
			echo "$@ : file not found !!"
		    fi
		fi
	    fi
	fi
    fi
    IGCM_debug_PopStack "ORCHIDEE_grep"
    return $RET
}



