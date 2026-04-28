#!/bin/bash
shopt -s extglob
inputArgs=("$@")
printenv
if [ -n "$filterends" ]; then
  echo " filterEnds.pl  ${inputArgs[@]} > $filterends"
  eval " filterEnds.pl  ${inputArgs[@]} > $filterends"
else
  echo "filterSam.pl ${inputArgs[@]}"
  eval "filterSam.pl ${inputArgs[@]}"
fi
