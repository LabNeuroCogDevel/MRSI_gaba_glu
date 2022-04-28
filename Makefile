out/gaba_glu.csv: shareGabaGluSheet.R readMRS.R 13MP20200207_LCMv2fixidx.csv lcm.xlsx
	./shareGabaGluSheet.R

out/slope_and_bar.pdf:
	./plot_GabaGlu_slope_and_bar.R
