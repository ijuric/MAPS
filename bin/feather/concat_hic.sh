outdir=$1
dataset_name=$2
hic_dir=$3
shift 3
dirs=$@
printf "${dirs[@]}"
all_chrs=()
hic_affix="hic."
function extract_chrs (){
        dir_arg=$1
        #echo "looking for: " $dir_arg
        long_pattern=$dir_arg"/*"$long_suffix
        short_pattern=$dir_arg"/*"$short_suffix
	hic_pattern=$dir_arg"/*"$hic_affix"*"
        #echo $files_pattern1
        hic_files_list=$(ls $hic_pattern)
	#printf "hic files: ${hic_files_list[@]}"
        #short_files_list=$(ls $short_pattern)
        #echo "\nbye\n"
        #echo $files1
        chrs_pairs=()
        for f in ${hic_files_list[@]}
        do
                #printf "$f\n"
                if [[ $f =~ chr[a-zA-Z0-9_]+.chr[a-zA-Z0-9_]+ ]]; then
                        strresult=${BASH_REMATCH[0]}
                        chrs_pairs+=($strresult)
			all_chrs+=($strresult)
                        #printf "result: $strresult\n"
                fi
        done
        IFS=$'\n'
        sorted_long_chrs=($(sort <<<"${long_chrs[*]}"))
        unset IFS
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
	chr_hic_files=()
	for f in ${dirs[@]}
	do 
		chr_hic_files+=($f"/*"$hic_affix$c)
	done
	outfile_hic=$outdir"/"$hic_dir"/"$dataset_name"."$hic_affix$c
	echo "calling: ${chr_hic_files[@]} > $outfile_hic 2>/dev/null \n"
	cat ${chr_hic_files[@]} > $outfile_hic 2>/dev/null
done
hic_files_list=$(ls $outdir"/"$hic_dir"/"*hic.*)
cat ${hic_files_list[@]} > $outdir"/"$dataset_name".hic.input"
