
after running some code blocks from `plot_assess_vary_windows.Rmd`

```
subset_df <- g1000_s30000_c0p01_out$dfoi[which(g1000_s30000_c0p01_out$dfoi$se==0.001 & g1000_s30000_c0p01_out$dfoi$avgr==0.6),]
max(subset_df$phasing_acc) #100
which.max(subset_df$phasing_acc) #returns the first index only sadly
```


for every combination of number of gametes, number of snps, coverage <br />
-- subset to a specific combo of specific average rate of recombination, and specific sequencing error rate <br />
-- -- find the window lengths that led to the maximum phasing accuracy (which.max == "argmax") <br />
 -- -- -- what's the window length that had (a) the most number of independent replicates that were an argmax (b) the largest such window length (shorter window lengths take longer to run) <br /> <br />
 -- -- -- -- then it's a linear model given all the input info to find the coefficients that best describe the ideal window length given the input parameters (# of gametes, # of snps, coverage, rate of recomb, sequencing error rate)