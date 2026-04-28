#!/bin/bash

baseDir=/data/bff-test
#check if $1 is given
if [ -n "$1" ]; then
    baseDir=$1
    echo "using base directory: $baseDir"
else
    echo "using default base directory: $baseDir"
fi


#check if the base directory exists
if [ ! -d "$baseDir" ]; then
    echo "Base directory does not exist"
    exit 0
fi
#remove old bam files
rm -f "$baseDir"/bam/*
rm -f "$baseDir"/bam.done/*

#remove old breakpoint files
rm -f "$baseDir"/breakpoint_files/*
rm -f "$baseDir"/final_breakpoint_files/*

#remove old summary files
rm -f "$baseDir"/old_summary*




