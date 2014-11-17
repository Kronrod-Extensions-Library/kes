awk -e '/RULE/ {S = ""; for (i = 3; i <= 3+$2; i++) {S = S $(i) " ";}; print S;}' recs*.dat | sort | uniq > all_rules.dat
