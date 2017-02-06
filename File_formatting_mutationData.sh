for i in $(ls *.txt) ; do awk -F'\t' '{print ($1,$9,$16) }' $i | grep -v -E "Silent|RNA" | sed 's/ /\t/g'  >> "Filtered_"$i; done 
for i in $(ls Filtered*) ; do  cut -d"-" -f1,2,3,4 $i >> "filtered_"$i; done
rm Filtered*
