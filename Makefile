# inputs from LCModel and visual inspection 
RAW := 13MP20200207_LCMv2fixidx.csv lcm.xlsx

.PHONY: all
all: out/gaba_glu.csv out/slope_and_bar.pdf

out/gaba_glu.csv: shareGabaGluSheet.R readMRS.R ${RAW}
	./shareGabaGluSheet.R

out/slope_and_bar.pdf: ${RAW}
	./plot_GabaGlu_slope_and_bar.R
