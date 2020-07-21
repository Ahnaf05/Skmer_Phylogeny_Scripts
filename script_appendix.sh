#!/bin/bash
#declare -a gene_array
declare -a file_name
declare -A tmp_GTR
declare -A tmp
declare -A x_AC
declare -A x_AG
declare -A x_AT
declare -A x_CG
declare -A x_CT
declare -A x_GT
count=0
gtr_index_count=0
GTRfile="GTRMatrix.txt"
mkdir -p "save"
rm $GTRfile
read_matrix() {
    local i=0
    local line
    local j
    # Ignore the first 2 lines containing size of the matrix
    while read -r line; do
        j=0
        # split on spaces
        for v in `echo $line`; do
            tmp[$i,$j]="$v"
            j=$((j+1))
        done
        i=$((i+1))
    done
}
##for getting original files##
rewrite()
{
	local j=0
	for i in `ls -1 ref_dir`
	do
		cp "save/$i" "ref_dir/$i"
		j=$((j+1))
	done
}
copy_matrix()
{
	for((i=0;i<count;i++))
	do
		for((j=0;j<count;j++))
		do
			$1[$i,$j]=$tmp[$i,$j]
		done
	done
}
#read_matrix<ref-dist-mat.txt
#echo ${tmp[1,2]}
for filename in `ls -1 ref_dir`
do
	#gene_array[count]=$(<"ref_dir/$filename")
	cp "ref_dir/$filename" "save/$filename" 
	file_name[count]=$filename
	#echo "${gene_array[count]}"
	count=$((count+1))
done
### calculate x_ac ##########

## Replace G and T with C ####
for i in `ls -1 ref_dir`
do
	sed -i 's/G/C/g' "ref_dir/$i"
	sed -i 's/T/C/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=0;i<=count;i++))
do
	for((j=0;j<=count;j++))
	do
		x_AC[$i,$j]=${tmp[$i,$j]}
	done
done
rewrite

## Replace G and T with A ####
for i in `ls -1 ref_dir`
do
	sed -i 's/G/A/g' "ref_dir/$i"
	sed -i 's/T/A/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=1;i<=count;i++))
do
	for((j=1;j<=count;j++))
	do
		tmp1=${x_AC[$i,$j]}
		tmp2=${tmp[$i,$j]}
		x_AC[$i,$j]=`awk "BEGIN {print $tmp1+$tmp2; exit}"`
	done
done
rewrite

### Replace G with T###
for i in `ls -1 ref_dir`
do
	sed -i 's/G/T/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=1;i<=count;i++))
do
	for((j=1;j<=count;j++))
	do
		tmp1=${x_AC[$i,$j]}
		tmp2=${tmp[$i,$j]}
		x_AC[$i,$j]=`awk "BEGIN {print $tmp1-$tmp2; exit}"`
	done
done
rewrite

#######calculate x_ag##########

## Replace C and T with G ####
for i in `ls -1 ref_dir`
do
	sed -i 's/C/G/g' "ref_dir/$i"
	sed -i 's/T/G/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=0;i<=count;i++))
do
	for((j=0;j<=count;j++))
	do
		x_AG[$i,$j]=${tmp[$i,$j]}
	done
done
rewrite

## Replace C and T with A ####
for i in `ls -1 ref_dir`
do
	sed -i 's/C/A/g' "ref_dir/$i"
	sed -i 's/T/A/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=1;i<=count;i++))
do
	for((j=1;j<=count;j++))
	do
		tmp1=${x_AG[$i,$j]}
		tmp2=${tmp[$i,$j]}
		x_AG[$i,$j]=`awk "BEGIN {print $tmp1+$tmp2; exit}"`
	done
done
rewrite

### Replace C with T###
for i in `ls -1 ref_dir`
do
	sed -i 's/C/T/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=1;i<=count;i++))
do
	for((j=1;j<=count;j++))
	do
		tmp1=${x_AG[$i,$j]}
		tmp2=${tmp[$i,$j]}
		x_AG[$i,$j]=`awk "BEGIN {print $tmp1-$tmp2; exit}"`
	done
done
rewrite

#######calculate x_at##########

## Replace C and G with T ####
for i in `ls -1 ref_dir`
do
	sed -i 's/C/T/g' "ref_dir/$i"
	sed -i 's/G/T/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=0;i<=count;i++))
do
	for((j=0;j<=count;j++))
	do
		x_AT[$i,$j]=${tmp[$i,$j]}
	done
done
rewrite

## Replace C and G with A ####
for i in `ls -1 ref_dir`
do
	sed -i 's/C/A/g' "ref_dir/$i"
	sed -i 's/G/A/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=1;i<=count;i++))
do
	for((j=1;j<=count;j++))
	do
		tmp1=${x_AT[$i,$j]}
		tmp2=${tmp[$i,$j]}
		x_AT[$i,$j]=`awk "BEGIN {print $tmp1+$tmp2; exit}"`
	done
done
rewrite

### Replace C with G###
for i in `ls -1 ref_dir`
do
	sed -i 's/C/G/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=1;i<=count;i++))
do
	for((j=1;j<=count;j++))
	do
		tmp1=${x_AT[$i,$j]}
		tmp2=${tmp[$i,$j]}
		x_AT[$i,$j]=`awk "BEGIN {print $tmp1-$tmp2; exit}"`
	done
done
rewrite

### calculate x_cg ##########

## Replace A and T with G ####
for i in `ls -1 ref_dir`
do
	sed -i 's/A/G/g' "ref_dir/$i"
	sed -i 's/T/G/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=0;i<=count;i++))
do
	for((j=0;j<=count;j++))
	do
		x_CG[$i,$j]=${tmp[$i,$j]}
	done
done
rewrite

## Replace A and T with C ####
for i in `ls -1 ref_dir`
do
	sed -i 's/A/C/g' "ref_dir/$i"
	sed -i 's/T/C/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=1;i<=count;i++))
do
	for((j=1;j<=count;j++))
	do
		tmp1=${x_CG[$i,$j]}
		tmp2=${tmp[$i,$j]}
		x_CG[$i,$j]=`awk "BEGIN {print $tmp1+$tmp2; exit}"`
	done
done
rewrite

### Replace A with T###
for i in `ls -1 ref_dir`
do
	sed -i 's/A/T/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=1;i<=count;i++))
do
	for((j=1;j<=count;j++))
	do
		tmp1=${x_CG[$i,$j]}
		tmp2=${tmp[$i,$j]}
		x_CG[$i,$j]=`awk "BEGIN {print $tmp1-$tmp2; exit}"`
	done
done
rewrite

### calculate x_ct ##########

## Replace A and G with T ####
for i in `ls -1 ref_dir`
do
	sed -i 's/A/T/g' "ref_dir/$i"
	sed -i 's/G/T/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=0;i<=count;i++))
do
	for((j=0;j<=count;j++))
	do
		x_CT[$i,$j]=${tmp[$i,$j]}
	done
done
rewrite

## Replace A and G with C ####
for i in `ls -1 ref_dir`
do
	sed -i 's/A/C/g' "ref_dir/$i"
	sed -i 's/G/C/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=1;i<=count;i++))
do
	for((j=1;j<=count;j++))
	do
		tmp1=${x_CT[$i,$j]}
		tmp2=${tmp[$i,$j]}
		x_CT[$i,$j]=`awk "BEGIN {print $tmp1+$tmp2; exit}"`
	done
done
rewrite

### Replace A with G###
for i in `ls -1 ref_dir`
do
	sed -i 's/A/G/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=1;i<=count;i++))
do
	for((j=1;j<=count;j++))
	do
		tmp1=${x_CT[$i,$j]}
		tmp2=${tmp[$i,$j]}
		x_CT[$i,$j]=`awk "BEGIN {print $tmp1-$tmp2; exit}"`
	done
done
rewrite

### calculate x_gt ##########

## Replace A and C with T ####
for i in `ls -1 ref_dir`
do
	sed -i 's/A/T/g' "ref_dir/$i"
	sed -i 's/C/T/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
echo ${tmp[1,2]}
for((i=0;i<=count;i++))
do
	for((j=0;j<=count;j++))
	do
		x_GT[$i,$j]=${tmp[$i,$j]}
	done
done

rewrite

## Replace A and C with G ####
for i in `ls -1 ref_dir`
do
	sed -i 's/A/G/g' "ref_dir/$i"
	sed -i 's/C/G/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
for((i=1;i<=count;i++))
do
	for((j=1;j<=count;j++))
	do
		tmp1=${x_GT[$i,$j]}
		tmp2=${tmp[$i,$j]}
		x_GT[$i,$j]=`awk "BEGIN {print $tmp1+$tmp2; exit}"`
	done
done
rewrite

### Replace A with C###
for i in `ls -1 ref_dir`
do
	sed -i 's/A/C/g' "ref_dir/$i"
	#cat "ref_dir/$i"
done
skmer reference ref_dir
read_matrix<ref-dist-mat.txt
echo ${tmp[1,2]}
for((i=1;i<=count;i++))
do
	for((j=1;j<=count;j++))
	do
		tmp1=${x_GT[$i,$j]}
		tmp2=${tmp[$i,$j]}
		x_GT[$i,$j]=`awk "BEGIN {print $tmp1-$tmp2; exit}"`
	done
done
rewrite
for((i=1;i<=count-1;i++))
do
	for((j=i+1;j<=count;j++))
	do
		idx1=$((i-1))
		idx2=$((j-1))
		echo "GTR ARRAY $i (${file_name[idx1]}),$j (${file_name[idx2]})">>$GTRfile
		echo -e "\n">>$GTRfile
		echo -e "0  \c">>$GTRfile
		echo -e "${x_AC[$i,$j]} \c">>$GTRfile
		echo -e "${x_AG[$i,$j]} \c">>$GTRfile
		echo -e "${x_AT[$i,$j]} \c">>$GTRfile
		tmp_GTR[0,0]=0
		tmp_GTR[0,1]=${x_AC[$i,$j]}
		tmp_GTR[0,2]=${x_AG[$i,$j]}
		tmp_GTR[0,3]=${x_AT[$i,$j]}
		echo -e "\n">>$GTRfile
		echo -e "${x_AC[$i,$j]} \c">>$GTRfile
		echo -e "0  \c">>$GTRfile
		echo -e "${x_CG[$i,$j]} \c">>$GTRfile
		echo -e "${x_CT[$i,$j]} \c">>$GTRfile
		tmp_GTR[1,0]=${x_AC[$i,$j]}
		tmp_GTR[1,1]=0
		tmp_GTR[1,2]=${x_CG[$i,$j]}
		tmp_GTR[1,3]=${x_CT[$i,$j]}
		echo -e "\n">>$GTRfile
		echo -e "${x_AG[$i,$j]} \c">>$GTRfile
		echo -e "${x_CG[$i,$j]} \c">>$GTRfile
		echo -e "0  \c">>$GTRfile
		echo -e "${x_GT[$i,$j]} \c">>$GTRfile
		tmp_GTR[2,0]=${x_AG[$i,$j]}
		tmp_GTR[2,1]=${x_CG[$i,$j]}
		tmp_GTR[2,2]=0
		tmp_GTR[2,3]=${x_GT[$i,$j]}
		echo -e "\n">>$GTRfile
		echo -e "${x_AT[$i,$j]} \c">>$GTRfile
		echo -e "${x_CT[$i,$j]} \c">>$GTRfile
		echo -e "${x_GT[$i,$j]} \c">>$GTRfile
		echo -e "0  \c">>$GTRfile
		tmp_GTR[3,0]=${x_AT[$i,$j]}
		tmp_GTR[3,1]=${x_CT[$i,$j]}
		tmp_GTR[3,2]=${x_GT[$i,$j]}
		tmp_GTR[3,3]=0
		echo -e "\n">>$GTRfile
		gtr_index_count=$((gtr_index_count+1))
	done
done

#echo $count
