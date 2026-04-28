#!/bin/bash
working_done="/data/bff-test/data/working.done"
echo "Cleaning up..."
[ -z "$basecallext" ] && basecallext="fast5"
[ -z "$fastqext" ] && fastqext="fastq"
function mkdirp() {
    local dir=$1
    if [ -d "$dir" ]; then
        echo "directory exists"
    else
        echo "directory does not exist - creating it"
        mkdir -p "$dir" && sync 
        rm -rf "$dir"/* && sync 
    fi
}
function rmrf() {
    local dir=$1
    echo "removing $dir"
    while [ -d "$dir" ]; do
        echo "waiting for $dir to be removed"
        #run in background and wait for the pid to finish   
        rm -rf "$dir" & 
        wait $!
        sleep 2
    done
}
if [ -n "$workingDir" ]; then
  #mark that the working directory is done
  #loop through workingDir and touch a file in the working_done directory
  echo "Marking files in $workingDir as done..."
  # Create working_done directory if it doesn't exist
  if [ ! -d "$working_done" ]; then
    echo "Creating working_done directory: $working_done"
    mkdir -p "$working_done"
  fi
  # Find all files in workingDir and touch them in working_done
  if [ -d "$workingDir" ] && [ "$(ls -A "$workingDir" 2>/dev/null)" ]; then
    while IFS= read -r file; do
      filename=$(basename "$file")
      echo "Marking $filename as done"
      touch "$working_done/$filename"
      sync
    done < <(find "$workingDir" -type f)
    echo "All files in $workingDir marked as done"
  else
    echo "No files found in $workingDir"
  fi

fi
if [ -n "$outputDir" ]; then 
  echo "moving $outputDir/*.$fastqext"
  if compgen -G "$outputDir/*.$fastqext" > /dev/null; then
     mkdir -p "$outputDir.done/" && eval "mv -f $outputDir/*.$fastqext $outputDir.done/."
  else
    echo "no matches"
  fi
fi
[ -n "$modelDir" ] && echo "removing $modelDir/*" && rm -rf "$modelDir/*"
[ -n "$referenceDir" ] && echo "removing $referenceDir/*" && rm -rf "$referencedir/*"
exit 0;