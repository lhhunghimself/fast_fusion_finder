#!/bin/bash
echo "runDorado.sh"
shopt -s extglob

#Get the filter file if it exists
#!/bin/bash

# Check if the $filterlist variable is set, is a file and is empty - then just write an empty fastq file
printenv
if [[ -n "$filterlist" && -f "$filterlist" && ! -s "$filterlist" ]]; then
    # Check if $fastqfilename variable is set
    echo "filterlist is empty"
    if [[ -n "$fastqfilename" ]]; then
        echo "writing empty file to $outputDir/$fastqfilename"
        # Write/replace with an empty file
        : > "$outputDir/$fastqfilename"
        exit 0
    else
        echo "Error: fastqfilename environment variable is not set."
        exit 1
    fi
fi
inputArgs=("$@")
#remove the first and last quotes in doradocmd
doradocmd=$(echo "$doradocmd" | sed 's/^"//; s/"$//')

# Set maximum number of retries
MAX_RETRIES=20
attempt=1
success=false
sleep_time=10

echo "$doradocmd $modeldir $inputDir ${inputArgs[@]} > $outputDir/$fastqfilename"

# Try running the command, retry if it fails up to MAX_RETRIES
while [ $attempt -le $MAX_RETRIES ] && [ "$success" = false ]; do
    echo "Attempt $attempt of $MAX_RETRIES..."
    
    # Run the command and capture exit status
    eval "$doradocmd $modeldir $inputDir ${inputArgs[@]}" > "$outputDir/$fastqfilename"
    exit_status=$?
    
    if [ $exit_status -eq 0 ]; then
        echo "Command succeeded on attempt $attempt"
        success=true
        exit $exit_status
    else
        echo "Command failed with exit status $exit_status on attempt $attempt"
        if [ $attempt -lt $MAX_RETRIES ]; then
            echo "Retrying..."
            sleep $sleep_time  # Wait 5 seconds before retrying
        else
            echo "Maximum retry attempts reached. Exiting."
            exit $exit_status
        fi
    fi
    #add sleep 5 seconds between attempts
    sleep $sleep_time
    ((attempt++))
done
