grep "Relaxation/intensification parameter" 3Significance/OG*.out > kvalue
grep "Relaxation/intensification parameter" 2no_significance/OG*.out > kvalue1
grep "Likelihood ratio test " 3Significance/OG*.out > p2
grep "Likelihood ratio test " 2no_significance/OG*.out > p1
grep "non-synonymous/synonymous rate ratio for " 3Significance/OG*.out > dnds
grep "Evidence for" 3Significance/OG*.out > selection
