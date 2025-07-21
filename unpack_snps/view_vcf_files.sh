ls ~/gtrd/allele_counts/ | awk -F'_' '{print $1}' | sort | uniq -c | sort -nr
