* `run_assess_varying_window.sh` -- loop over wanted iterations of window size, seq error, avg recomb, random seed given a number of gametes, a number of snps, and a coverage
* `submit_assess_varying_window.sh` -- actually submits the job to slurm
* `assess_varying_window.R` -- run's rhapsodi on a genModel output given a window size and run's the assessing of phasing, etc from the rhapsodi output
* `check_assess_varying_window.R` -- used to look for iterations that didn't successfully complete from `run_assess_varying_window.sh` & `assess_varying_window.R`; also used to build a csv file that lists number of gametes, number of snps, coverage, avg recomb, phasing acc, phasing ser, etc.
* `count_assess_varying_window.sh` -- gives a count of the number of iterations that `run_assess_varying_window.sh` & `assess_varying_window.R` will produce, given as an argument to `check_assess_varying_window.R`

running on MARCC for
* g2500_s5000_c0.693
* g1000_s30000_c0.01
* g50_s100000_c0.511

Downloading the `.csv` files from `check_assess_varying_window.R` locally. Added the centimorgan/window readout to the table there with `plot_assess_varying_window.Rmd`

We expected that we should always see a parabola when observing ln(centimorgans/window) on x axis and phasing accuracy on the y-axis such that for too large of window sizes, phasing isn't accurate, and for too small of window sizes, there's not enough shared information and therefore inaccurate phasing accuracy.

However, after expanding the window lengths we used in running rhapsodi on the generated data, I observed that for higher coverages, the exponential decay on phasing accuracy still exists. This is likely because the parabolic behavior for smaller window sizes isn't a problem given a high enough coverage. (Q: what is a high enough coverage?)

ran more combos on MARCC
* g2500_s5000_c0.001
* g2500_s5000_c0.01
* g2500_s5000_c0.1
* g2500_s5000_c0.223
* g1000_s30000_c0.001
* g1000_s30000_c0.1
* g1000_s30000_c0.223

Next steps are to find the ideal window size. I've elaborated on this in `set_window_meeting03012022.md`
