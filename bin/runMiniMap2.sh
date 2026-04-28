#!/bin/bash
makeArrayString() {
	echo $1 | sed 's/[][]//g' | sed 's/\,/ /g'
}
unquotedFile() {
	#remove any quotes that are escaped \"
	echo $* | sed 's/\\\"//g' | sed 's/\"//g'
}
align() {
	local file="$1"
	local unquoted filename bambase outputfile exit_code
	unquoted=$(unquotedFile "$file")
	filename=$(basename -- "$unquoted")
	bambase="$outputdir/${filename%%.*}"
	# check for case of empty file
	if [ ! -s "$unquoted" ]; then
		echo "file $unquoted is empty - skipping"
		return 0
	fi
	if [ -n "$saveassam" ]; then
		outputfile="$bambase.sam"
		# try to align the file nTries times
		exit_code=1
		for i in $(seq 1 "$nTries"); do
			echo "minimap2 -ax map-ont $optionalflags ${inputArgs[*]} $indexfile $unquoted > $outputfile"
			minimap2 -ax map-ont $optionalflags "${inputArgs[@]}" "$indexfile" "$unquoted" > "$outputfile"
			exit_code=$?
			if [ $exit_code -eq 0 ]; then
				break
			fi
			sleep "$nTriesSleep"
		done
		if [ $exit_code -ne 0 ]; then
			echo "failed to align $file after $nTries tries"
			return $exit_code
		fi
		if [ -n "$mergesam" ]; then
			if [ -f "$outputdir/$mergesam" ]; then
				grep -v '^[#@]' "$outputfile" >> "$outputdir/$mergesam"
			elif [ -f "$outputfile" ] && [ "$outputfile" != "$outputdir/$mergesam" ]; then
				cp "$outputfile" "$outputdir/$mergesam"
			else
				touch "$outputdir/$mergesam"
			fi
		fi
	else
		#try to align the file nTries times
		exit_code=1
		for i in $(seq 1 "$nTries"); do
			echo "minimap2 -ax map-ont $optionalflags ${inputArgs[@]} $indexfile $unquoted | samtools view -bS | samtools sort -O BAM | tee $bambase.bam | samtools index - $bambase.bai"  	
			minimap2 -ax map-ont $optionalflags "${inputArgs[@]}" "$indexfile" "$unquoted" | samtools view -bS | samtools sort -O BAM | tee "$bambase.bam" | samtools index - "$bambase.bai"
			exit_code=$?
			if [ $exit_code -eq 0 ]; then
				break
			fi
			sleep "$nTriesSleep"
		done
		if [ $exit_code -ne 0 ]; then
			echo "failed to align $file after $nTries tries"
			return $exit_code
		fi
	fi
}
shopt -s extglob
inputArgs=("$@")
set -ef -o pipefail
mkdir -p $outputdir
nTries=3
nTriesSleep=10
files=($(makeArrayString $inputfiles))
for file in "${files[@]}"; do
	align "$file"
done
#check if the output file exists - otherwise write an empty file
[ -n "$saveassam" ] && [ -n "$mergesam" ] && [ ! -f "$outputdir/$mergesam" ] && echo "writing empty file to $outputdir/$mergesam" && touch "$outputdir/$mergesam"
exit 0;