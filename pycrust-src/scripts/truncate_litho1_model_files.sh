#! /bin/bash

## This script will process the litho1 model files
## and leave behind only the lithosphere layers 

options='r:t:h'

rset=false
tset=false

while getopts $options option
do
    case $option in
        r  )    rawpath=$OPTARG ;
                rset=true ;
                echo "raw model path: " $OPTARG;;
        t  )    truncatedpath=$OPTARG ;
                tset=true ;
                echo "truncated model path: " $OPTARG;;
        h  )    echo "truncate... -r PATH_TO_RAW -t PATH_TO_TRUNCATED" ;
                exit ;;
        \? )    echo "Unknown option: -$OPTARG" ;;
        * )     exit "Missing option argument: -$OPTARG" ;
                exit -1 ;;
    esac
done

set -e


if  ! ( $rset && $tset)
then
    echo "Provide path to raw model (-r) and a place to store the truncated model (-t)"
    exit 1
fi

mkdir -p $truncatedpath

cd $rawpath
for file in *.model
do
    sed -n -e '/ASTHENO-TOP/,$p' $file > $truncatedpath/${file}_tr
done
