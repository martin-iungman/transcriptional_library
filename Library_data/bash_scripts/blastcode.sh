#!/bin/sh

while IFS= read -r line; do
  if [[ $line =~ ^\> ]]; then     seq_id=${line#>}
    sequence=$(awk -v id="$seq_id" -v RS='>' '$0 ~ id {print $0}' "$query_file")    
    query_tmp=$(mktemp);     echo ">query" > "$query_tmp";     echo "$sequence" >> "$query_tmp"
    blastn -query "$query_tmp" -db "$db_dir/reference" -outfmt '6 sstart send' | awk -v id="$seq_id" -v OFS="\t" '{print id, $1, $2}' >> "$output_file"
fi
done
