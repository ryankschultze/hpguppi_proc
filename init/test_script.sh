#!/bin/bash
# Test shell script

if [ "$1" = 'cbf' ]
then
    if [ $# -eq 2 ]
    then
        echo '2 arguments.'
    elif [ $# -eq 3 ]
    then
        echo '3 arguments.'
    elif [ $# -eq 4 ]
    then
        echo '4 arguments.'
    elif [ $# -eq 5 ]
    then
        echo '5 arguments.'
    else
        echo 'No arguments.'
    fi
else
    echo 'No mode chosen.'
    exit
fi


