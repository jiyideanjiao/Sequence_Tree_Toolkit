grep " 0" OG*.phy > 0.txt
grep " 1" OG*.phy > 1.txt
grep " 2" OG*.phy > 2.txt
grep " 3" OG*.phy > 3.txt
grep " 4" OG*.phy > 4.txt
grep " 5" OG*.phy > 5.txt
grep " 6" OG*.phy > 6.txt
grep " 7" OG*.phy > 7.txt
grep " 8" OG*.phy > 8.txt
grep " 9" OG*.phy > 9.txt
cat *.txt > all.txt
sed -i 's/:/ /g' all.txt
awk '!a[$1]++' all.txt > all_uniq1.txt
awk '{print $3}' all_uniq1.txt > all_uniq2.txt
sort all_uniq2.txt | uniq -c > count_freq1.txt
awk '{print $2,$1}' count_freq1.txt > count_freq2.txt
sed -i 's/ /,/g' count_freq2.txt
mv count_freq2.txt count_freq2.csv
