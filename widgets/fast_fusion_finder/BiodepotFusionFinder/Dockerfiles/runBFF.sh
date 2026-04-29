#!/bin/bash
shopt -s extglob
inputArgs=("$@")
printenv
if [ -n "$filterends" ]; then
  echo " filterEnds.pl  ${inputArgs[@]} > $filterends"
  eval " filterEnds.pl  ${inputArgs[@]} > $filterends"

  if [ -n "$BAM" ] && [ -n "$REFERENCE" ]; then
    bp_file=""
    for ((i=0; i<${#inputArgs[@]}; i++)); do
      case "${inputArgs[i]}" in
        -b|--breakpointFile)        bp_file="${inputArgs[i+1]}";;
        --breakpointFile=*)         bp_file="${inputArgs[i]#--breakpointFile=}";;
        -b=*)                       bp_file="${inputArgs[i]#-b=}";;
      esac
    done
    : "${bp_file:?could not find -b/--breakpointFile in args for breakpointPileup.pl}"

    pileup_out="${pileup:-${filterends}.fa}"
    extra_args=()
    [ -n "$WINDOW" ]       && extra_args+=(--window "$WINDOW")
    [ -n "$pileupTable" ]  && extra_args+=(--pileupOut "$pileupTable")
    echo "breakpointPileup.pl --reads $filterends --bam $BAM --ref $REFERENCE -b $bp_file ${extra_args[@]} > $pileup_out"
    breakpointPileup.pl --reads "$filterends" --bam "$BAM" --ref "$REFERENCE" -b "$bp_file" "${extra_args[@]}" > "$pileup_out"
  else
    echo "BAM and/or REFERENCE not set; skipping breakpointPileup.pl"
  fi
else
  echo "/data/debug/filterSam.pl ${inputArgs[@]}"
  eval "/data/debug/filterSam.pl ${inputArgs[@]}"
fi
