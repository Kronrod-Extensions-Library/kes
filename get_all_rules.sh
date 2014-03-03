awk -e '/RULE/ {S = ""; for (i = 5; i < 5+$2; i++) {S = S $(i) " ";}; print S;}' recs*.dat | sort | uniq > all_rules.dat
