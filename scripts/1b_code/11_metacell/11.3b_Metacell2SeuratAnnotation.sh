#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 1st July 2020
# Title: 11.3b_Metacell2SeuratAnnotation.sh
# Goal: Prepare annotations from Metacell output
# Usage: 11.3b_Metacell2SeuratAnnotation.sh metacell_dir
# Where: metacell_dir is a directory containing out_fig, which contains "metacell.list.txt"
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Initialisation
# ------------------------------------------------------------------
display_usage() {
	echo "ERROR: 1 argument expected for 11.3b_Metacell2SeuratAnnotation.sh - exiting!"

	echo -e "\nUsage:\ns11.3b_Metacell2SeuratAnnotation.sh /path/to/metacell_dir\n"
	}

# Check if correct number of arguments (n = 1) provided
if [ $# != 1 ]; then
	display_usage
    exit 1
fi

metacell_dir=$1
cd "$metacell_dir"/out_fig

# delete the header/rowname at top
tail -n+2 metacell.list.txt | awk '{print $2, $1}' | sed 's/ /	/g' | sort | sort -nk1,1 > metacell.list.txt_v2
tail -n+2 metacell.groupings.txt > metacell.groupings.txt_v2
join -a1 -e"unassigned" -o '1.2,0,2.2' metacell.list.txt_v2 metacell.groupings.txt_v2 | sed 's/ /	/g' > metacell_metadata.txt
join -a1 -1 3 -2 2 -e"unassigned" -o '1.1,1.2,0,2.1' <(sort -k3,3 metacell_metadata.txt ) <(cut -f 2,4 metacell.colourchart.txt | sort -k2,2 ) | sed 's/ /	/g' > metacell_metadataGroups.txt
sort -u metacell_metadataGroups.txt -o metacell_metadataGroups.txt

rm metacell.list.txt_v2 metacell.groupings.txt_v2
rm metacell_metadata.txt
