outdir=$1
dataset_name=$2
shift 2
dirs=$@
printf "${dirs[@]}"
all_chrs=()
long_suffix=".long.intra.bedpe"
short_suffix=".shrt.vip.bed"
function extract_chrs (){
        dir_arg=$1
        #echo "looking for: " $dir_arg
        long_pattern=$dir_arg"/*"$long_suffix
        short_pattern=$dir_arg"/*"$short_suffix
        #echo $files_pattern1
        long_files_list=$(ls $long_pattern)
	#printf "long files: ${long_files_list[@]}"
        short_files_list=$(ls $short_pattern)
        #echo "\nbye\n"
        #echo $files1
        long_chrs=()
        for f in ${long_files_list[@]}
        do
                #printf "$f\n"
                if [[ $f =~ chr[a-zA-Z0-9_]+ ]]; then
                        strresult=${BASH_REMATCH[0]}
                        long_chrs+=($strresult)
			all_chrs+=($strresult)
                        #printf "result: $strresult\n"
                fi
        done
        IFS=$'\n'
        sorted_long_chrs=($(sort <<<"${long_chrs[*]}"))
        unset IFS

        short_chrs=()
        for f in ${short_files_list[@]}
        do
                #printf "$f\n"
                if [[ $f =~ chr[a-zA-Z0-9_]+ ]]; then
                        strresult=${BASH_REMATCH[0]}
                        short_chrs+=($strresult)
			all_chrs+=($strresult)
                        #printf "result: $strresult\n"
                fi
        done
        IFS=$'\n'
        sorted_short_chrs=($(sort <<<"${short_chrs[*]}"))
        unset IFS
	#echo "${all_chrs[@]}"
}
for f in ${dirs[@]}
do
	extract_chrs $f
done
#printf "allchrsfound "${#all_chrs[@]}"\n"
#echo ${all_chrs[@]}
#printf "see: %s\n" "${all_chrs[@]}"
eval unique_chrs=($(printf "%q\n" "${all_chrs[@]}" | sort -u))
#unique_chrs=($(echo "${all_chrs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
#printf "hey: %s\n" "${a[@]}"
#printf "thisisoutdir: $outdir\n"
unset all_chrs
#printf "${all_chrs[@]}"

for c in ${unique_chrs[@]}
do
	chr_long_files=()
	chr_short_files=()
	for f in ${dirs[@]}
	do 
		chr_long_files+=($f"/*"$c$long_suffix)
		chr_short_files+=($f"/*"$c$short_suffix)
	done
	outfile_long=$outdir"/"$dataset_name"."$c$long_suffix
	outfile_short=$outdir"/"$dataset_name"."$c$short_suffix
	#echo "calling: ${chr_long_files[@]} > $outfile_long 2>/dev/null \n"
	cat ${chr_long_files[@]} > $outfile_long 2>/dev/null
	cat ${chr_short_files[@]} > $outfile_short 2>/dev/null
done
