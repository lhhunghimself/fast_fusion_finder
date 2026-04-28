#!/bin/bash

# remote_server="128.208.252.232"
# remote_dir="/srv/lhhung/output"
# local_dir="$PWD/input"
# working_dir="$PWD/working"
# sshUser=lhhung
nthreads=1
working_done="/data/bff-test/data/working.done"
working_copy_done="/data/bff-test/data/pod5.done"
running_dir="/data/bff-test/data/working.running"

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
function add_to_running_dir() {
    # Create running_dir if it doesn't exist
    if [ ! -d "$running_dir" ]; then
        echo "Creating running directory: $running_dir"
        mkdir -p "$running_dir"
        sync
    fi
    
    # Loop through all files in working_dir and touch them in running_dir
    if [ -d "$working_dir" ] && [ "$(ls -A "$working_dir" 2>/dev/null)" ]; then
        while IFS= read -r file; do
            filename=$(basename "$file")
            echo "Marking $filename as running"
            touch "$running_dir/$filename"
            sync
        done < <(find "$working_dir" -type f)
        echo "All files in working_dir marked as running"
    else
        echo "No files found in working_dir to mark as running"
    fi
}
function update_running_dir() {
    #remove any files in the running_dir that are not in the working_dir or the done_files directory    
    if [ -d "$running_dir" ]; then
        while IFS= read -r file; do
            filename=$(basename "$file")
            if [ ! -f "$working_dir/$filename" ] || [ -f "$working_done/$filename" ]; then
                echo "Removing $filename from running_dir"
                rm -f "$running_dir/$filename"
            fi
        done < <(find "$running_dir" -type f)
    fi
}

function update_working_dir() {
    if [ -d "$working_done" ]; then
        echo "working_done directory exists"
    else
        echo "working_done directory does not exist - creating it"
        mkdir -p "$working_done"
    fi
    # Create a dictionary of files in working_done directory
    declare -A done_files
    
    # Find all files in working_done directory
    if [ -d "$working_done" ] && [ "$(ls -A "$working_done" 2>/dev/null)" ]; then
        while IFS= read -r file; do
            filename=$(basename "$file")
            done_files["$filename"]=1
        done < <(find "$working_done" -type f)
    fi
    
    # Loop through files in working_dir and remove those that are in done_files
    if [ -d "$working_dir" ] && [ "$(ls -A "$working_dir" 2>/dev/null)" ]; then
        while IFS= read -r file; do
            filename=$(basename "$file")
            if [ "${done_files[$filename]}" = "1" ]; then
                echo "File $filename already processed, removing from working directory"
                #check if the file is in the working_copy_done directory - if not copy it to the working_copy_done directory
                if [ ! -f "${working_copy_done}/${filename}" ]; then
                    echo "File $filename not in working_copy_done directory, copying to it"
                    cp -r "$file" "$working_copy_done"
                    sync
                fi
                #remove the file from the working directory
                rm -f "$file"
                sync
                #remove the file from the running_dir directory
                rm -f "$running_dir/$filename"
                sync    
            fi
        done < <(find "$working_dir" -type f)
    fi
}
update_local_dir() {
    #make a dictionary of files from the done_files directory and the working_dir directory
    declare -A done_files
    declare -A working_files
    local local_files=()
    #check if there are local files
    if [ -d "$local_dir" ] && [ "$(ls -A "$local_dir" 2>/dev/null)" ]; then
        echo "local directory is not empty"
        local_files=$(ls "$local_dir")
    else
        echo "local directory is empty"
        return
    fi
    #get the list of files in the done_files directory
    local files=$(ls "$working_done")
    echo "files in working_done: $files"    
    #get the list of files in the working_dir directory
    local files=$(ls "$working_dir")
    echo "files in working_dir: $files"
    #combine the two lists
    local files=$(echo "$files" | sort -u)
    echo "combined files: $files"
    #create a dictionary of the files
    for file in $files; do
        done_files["$file"]=1
        working_files["$file"]=1
    done
    echo "done_files: $done_files"
    echo "working_files: $working_files"
    #loop through the files in the $local_dir directory and remove those that are in the done_files dictionary
    for file in $local_files; do
        if [ "${done_files[$file]}" = "1" ]; then
            echo "File $file already processed, removing from local directory"
            rm -f "$local_dir/$file"
            sync
        fi
    done
}
update_working_done() {
    #remove all the files in the working_done directory
    rm -rf "$working_done"/*
    sync    
}

# Check if ssh directory is present, unless localcopy is set
if [ -z "$localcopy" ]; then
    if [ -d "$HOME/.ssh" ]; then
        echo "ssh directory present"
    else
        echo "ssh directory not present - copying from $sshDir"
        cp -r "$sshDir" "$HOME/.ssh"
    fi
fi
 
# Create working directory if it doesn't exist
local_temp_dir="$local_dir.temp"
mkdirp "$working_dir"
mkdirp "$local_dir"
mkdirp "$running_dir"
mkdirp "$local_temp_dir"
#create a copying directory so that we know that someting is being copied
copying_dir="$local_dir.copying"
mkdirp "$copying_dir"




# Function to check if the input directory is empty of files ending with inputExt
dir_full() {
    local dir=$1
    local ext=$2
    if [ -z "$(ls -A "${dir}" | grep "${ext}$")" ]; then
        echo "Directory $dir is empty" >&2
        return 0 # directory is empty
    else
        echo "Directory $dir is not empty" >&2
        return 1 # directory is not empty
    fi
}

# Function to check if the file is being copied or exists locally or is being worked on
file_status_check() {
    local file=$1
    if [ -f "${working_dir}/${file}" ] || [ -f "${local_dir}/${file}" ] || [ -f "${copying_dir}/${file}" ] || [ -f "${working_done}/${file}" ] || [ -f "${working_copy_done}/${file}" ] ; then
            return 0 # file is being copied or exists locally or is being worked on
    else
        return 1 # file is not being copied or exists locally or is being worked on
    fi
}

# Function to perform the file copy operation
copy_file() {
    local file=$1
    local file_name=$(basename "$file")
    local temp_file="${local_temp_dir}/${file_name}.tmp"
    local temp_copying_dir="${copying_dir}/${file_name}"
    # Mark the file as being copied by creating a directory
    mkdirp "$temp_copying_dir"
    echo "will copy $file to $temp_file"
    # Copy the file to a temporary file in local directory
    if [ -n "$localcopy" ]; then
        if cp "${file}" "${temp_file}"; then
            echo "File ${file_name} copied to local directory."
            # Move the temporary file to the final file name
            cp "${temp_file}" "${local_dir}/${file_name}"
            rm -f "${temp_file}"
            # Remove the placeholder directory on success
            rmdir "$temp_copying_dir"
            return 0
        else
            echo "Failed to copy ${file_name}."
            # Remove the temporary file and placeholder directory if cp failed
            rm -f "${temp_file}"
            rmdir "$temp_copying_dir"
            return 1
        fi
    else
        if scp "${sshUser}@${remote_server}:${file}" "${temp_file}"; then
            echo "File ${file_name} copied to local directory."
            # Move the temporary file to the final file name
            echo "will copy ${temp_file} to ${local_dir}/${file_name}"
            cp "${temp_file}" "${local_dir}/${file_name}" && rm -f "${temp_file}"
            # Remove the placeholder directory on success
            rmdir "$temp_copying_dir"
            return 0
        else
            echo "Failed to copy ${file_name}."
            # Remove the temporary file and placeholder directory if scp failed
            rm -f "${temp_file}"
            rmdir "$temp_copying_dir"
            return 1
        fi
    fi
}
check_for_files_to_download() {
    echo "Checking for files to download"
    filesToBeDownloaded=()
    # Check each file
    count=0
    if [ -z "$localcopy" ]; then
        file_list=$(ssh "${sshUser}@${remote_server}" "find '${remote_dir}' -name '*.$inputExt'")
        echo "Remote files found: $file_list"
    else
        cmd="find '${local_scan_dir}' -name '*.$inputExt'"
        echo "cmd: $cmd"
        file_list=$(eval "$cmd")
        echo "Local files found: $file_list"
    fi
    for remote_file in $file_list; do
        file_name=$(basename "$remote_file")
        # Check if file is being copied or exists locally or has been processed
        if file_status_check "$file_name"; then
            continue
        else
            echo "file $file_name is not being copied or exists locally or has been processed"
            # Add file to the array and increment counter
            filesToBeDownloaded+=("$remote_file")
            ((count++))
        fi
        # Break the loop if the maximum number of files has been reached
        [ -n "$maxFiles" ] && [ $count -ge "$maxFiles" ] && break
    done

    # Download the files using up to nthreads if the list is not empty
    if [ ${#filesToBeDownloaded[@]} -gt 0 ]; then
        echo "will copy ${#filesToBeDownloaded[@]} files"
        for file in "${filesToBeDownloaded[@]}"; do
            copy_file "$file"
            echo "copied $file"
        done
        echo "All downloads completed."
    else
        if [ $((nwaits % 6)) -eq 0 ]; then
            # Find total time waited
            total_time=$((nwaits * 10))
            echo "No files to download. Slept for ${total_time} seconds..."
        fi
        nwaits=$((nwaits + 1))
        sleep 10
    fi
}
nwaits=0
# Main loop
echo "starting main loop"

while true; do
    # First check if working_done directory has processed files
    dir_full "$working_done" "$inputExt"
    working_done_is_full=$?
    
    if [ $working_done_is_full -eq 1 ]; then
        echo "Found processed files in working_done, updating directories..."
        update_working_dir
        update_local_dir
        update_running_dir
        update_working_done
        sleep 10
        continue
    fi

    # Then check if working directory has files to process
    dir_full "$working_dir" "$inputExt"
    working_dir_is_full=$?
    if [ $working_dir_is_full -eq 1 ]; then
        echo "Working directory has files to process"
        
        # Check if running_dir is empty (no job currently running)
        if [ ! -d "$running_dir" ] || [ -z "$(ls -A "$running_dir")" ]; then
            echo "No job running - sending signal to start a new job"
            add_to_running_dir
            sleep 10
            exit 0
        else
            echo "Job already running - waiting for it to finish"
            
            # Wait for the job to complete
            wait_count=0
            while true; do
                dir_full "$working_done" "$inputExt"
                if [ $? -eq 0 ]; then
                    # working_done is empty, check running directory
                    if [ -z "$(ls -A "$running_dir")" ]; then
                        echo "Job appears to have completed"
                        break
                    fi
                else
                    echo "Job completed with results"
                    break
                fi
                
                wait_count=$((wait_count + 1))
                if [ $((wait_count % 6)) -eq 0 ]; then
                    echo "Waited for $(($wait_count/6)) minutes for job to complete..."
                fi
                sleep 60
            done
        fi
        
        continue
    fi
    
    # At this point, we know working_dir is empty. Proceed with file download and transfer
    # without redundant empty checks until we try to copy files from local_dir to working_dir
    #check if the local directory is empty
    while [ -z "$(ls -A "${local_dir}" 2>/dev/null)" ]; do
        echo "Local directory is empty - sleeping for 10 seconds"
        check_for_files_to_download
        sleep 10
    done
    echo "Local directory has files - moving to working directory"
    cp -r "$local_dir"/* "$working_dir" && sync
    rmrf "$local_dir"
    echo "Removed local directory and created new empty local directory"
    sleep 10 && sync
    mkdirp "$local_dir"
    add_to_running_dir
    exit 0
done
