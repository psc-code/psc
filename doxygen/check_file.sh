#! /bin/bash

#########################################################
# check_file.sh - A shell script to aid in the doxygen  #
# adaption of the PSC source. Specify a path to one or  #
# more source files. Will check the doxygen compatible  #
# status of the documentation in the files and output   #
# any warnings.                                         #
#                                                       #
# Command line options other than -f should not be      #
# necesary unless the doxygen conf file psc_dox.conf is # 
# modified.                                             #
#########################################################


WARNCONF_D=""
BUILDCONF_D="`dirname "$0"`/psc_dox.conf"
CHFILES=""

SEP=":::::::::::"
SEP1="----------"
SEP2="~~~~~~~~~~"


usage()
{
    cat << EOF
usage: $0 <options>  -f [file1]  -f [file2] ... 
 Options are.. 
  -h
    Print this message and exit
  
 -f <path_to_file>
    Specify a file to check. 

 -b <doxygen_conf>
     Specify doxygen build file which would be used to
     build the actual documentation
     DEFAULT: ${BUILDCONF_D}
EOF
    exit 1
}

BUILDCONF=${BUILDCONF_D}

while getopts ":f:e:b:h" Option
do
    case $Option in
	f) CHFILES="${CHFILES} ${OPTARG}"
	    ;;
	e) ERRLOG=${OPTARG}
	    ;;
	b) BUILDCONF=${OPTARG}
	    ;;
	h) 
	    usage
	    ;;
	*) 
	    echo >&2 "Unrecognized option ${OPTION}" 
	    usage
	    ;;
	
    esac
done

if [ -z "${CHFILES}" ]
then
    echo "Please specify atleast one file to check!"
    usage
fi

echo ${SEP1}

for file in ${CHFILES}
do
    if [ ! -e "${file}" ]
    then
	echo "File ${file} not found!"
	echo $SEP1
	read line
	continue
    fi
    
    echo "File : ${file}"
    
    echo $SEP2

# Check for file tags
    check=`grep -c '[\]file' ${file}`
    if [ "${check}" -eq "0" ]
    then 
	echo "Missing '\file' tag."
	echo "Please add a special comment block to the end of the file containing"
	echo "'\file ${file} <description>'"
	echo "Without this tag, this script will be unable to find some undocumented elements."
	echo "Rerun $0 after adding the \file tag"
        echo $SEP2
    elif [ "${check}" -gt "1" ]
    then
	echo "More than one '\file' tag found."
	grep '[\]file' ${file}
	echo "Please check that you want all these tags."
        echo $SEP2
    fi
    

# Check for untagged FIXME's

    check=`grep -c "[ /c!]FIXME" ${file}`
    if [ ! "${check}" -eq "0" ]
    then
	echo "Detected ${check} FIXMEs which will not be converted to todos"
	echo "Please consider elaborating on the following FIXMEs, enclosing them in a doxygen comment block"
	echo "and adding the tag prefix (FIXME->\FIXME) so they will show up in the generated documentation"
	echo
	grep -n "[ /c!]FIXME" ${file}
	echo
    fi 
    
    
    echo ${SEP2}
# Check for warnings    
    
    WARNS=$(mktemp)

    cat ${BUILDCONF} | awk -F"=" '{if( /^EXTRACT_ALL|^GENERATE_/ ) { sub(/YES/,"NO",$0); print} else if( /^INPUT/ ) { sub(/..\//,"'"${file}"'",$0); print } else if( /^WARN_LOGFILE/ ) { print $1 "=" } else { print } }' | doxygen - >> /dev/null 2>${WARNS}

    NWARN=`grep -c "Warning" ${WARNS}`
    NNDOC=`grep -c "not documented" ${WARNS}`
    
    echo "Number of warnings: ${NWARN}"
    echo "Number of undocumented: ${NNDOC}"
    
    NCOMP=`grep -c "Compound" ${WARNS}`
    NMEM=`grep -c "Member" ${WARNS}`
    NFUNC=`grep "Member" ${WARNS} |grep -c "(function)"`
    NVAR=`grep "Member" ${WARNS} |grep -c "(variable)"`
    
    echo -e "\t Compounds: ${NCOMP}"
    echo -e "\t Members: ${NMEM}"
    echo -e "\t \t functions: ${NFUNC}"
    echo -e "\t \t variables: ${NVAR}" 

    echo $SEP1
    
    cat ${WARNS} | sort -nt ":" -k2 

    echo $SEP1

    echo "Press <enter> to continue"
    read line
    rm ${WARNS}
done

echo "Finished checking files: ${CHFILES}"