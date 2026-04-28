#!/bin/bash

compare_two_files() {
    local file1=$1
    local file2=$2
    unset email_file
    #if file1 does not exist or is empty then return
    if [ ! -s "$file1" ]; then
        echo "Error: File not found or empty - $file1"
        return 1
    fi
    #if file2 does not exist or is empty then copy file1 to file2 and return as if different
    if [ ! -s "$file2" ]; then
        echo "Warning: File not found or empty - $file2"
        cp "$file1" "$file2"
        email_file="$file2"
        return 0
    fi 
    # Check if the files are the same
    if cmp -s "$file1" "$file2"; then
        echo "Files are the same"
        return 1
    else
        echo "Files are different"
        cp "$file1" "$file2"
        email_file="$file2"
        return 0
    fi
}


# read the email recipients from the environment variable 

#remove any \" from the email recipients using sed
EMAIL_RECIPIENTS=$(echo "$EMAIL_RECIPIENTS" | sed 's/\\//g')
EMAIL_RECIPIENTS_CONTENT=$(echo "$EMAIL_RECIPIENTS" | sed 's/[][]//g; s/\"//g' | tr '\n' ' ')

# Split the recipients into an array
IFS="," read -ra ADDR_ARRAY <<< "$EMAIL_RECIPIENTS_CONTENT"
IFS=" "
ADDR_ARRAY=("${ADDR_ARRAY[@]//#/}")

# copy the configuration file 
cp $CONFIG_FILE /etc/ssmtp/ssmtp.conf

compare_two_files $SUMMARY_FILE $OLD_SUMMARY_FILE
if [ $? -eq 0 ]; then
  # Send email to each recipient using ssmtp
  for RECIPIENT in "${ADDR_ARRAY[@]}"; do
    # Skip empty strings (which may result from parsing)
    if [ -n "$RECIPIENT" ]; then
      echo "Sending email to $RECIPIENT"    
      # Use a here-document to include the From, To, Subject headers, and then concatenate the body
      {
        echo -e "From: BFF-MAIL <lhhunghimself@gmail.com>"
        echo -e "To: $RECIPIENT"
        echo -e "Subject: $EMAIL_SUBJECT"
        echo -e ""
        cat "$SUMMARY_FILE"
      } | ssmtp $RECIPIENT
    fi
  done
fi
exit 0
