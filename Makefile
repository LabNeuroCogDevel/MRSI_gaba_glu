# inputs from LCModel and visual inspection 
RAW := 13MP20200207_LCMv2fixidx.csv lcm.xlsx

.PHONY: all
all: out/gamadj_wide.csv # imgs/thresholding_cnt_met_region.png

# 20230531 - previously used
#all: out/gaba_glu.csv out/slope_and_bar.pdf

# adj_wide is last made
out/long_thres.csv out/gamadj_long.csv: out/gamadj_wide.csv
out/gamadj_wide.csv: 13MP20200207_LCMv2fixidx.csv ./adjust_all.R
	./adjust_all.R

imgs/thresholding_cnt_met_region.png: imgs/gam_adjusted_Vs_Cr.png
imgs/gam_adjusted_Vs_Cr.png: out/long_thres.csv out/gamadj_wide.csv ./visualize.R
	./visualize.R


out/gaba_glu.csv: shareGabaGluSheet.R readMRS.R ${RAW}
	./shareGabaGluSheet.R

out/slope_and_bar.pdf: ${RAW}
	./plot_GabaGlu_slope_and_bar.R

# input from hand placed PFC coordinates run through lcmodel
13MP20200207_LCMv2fixidx.csv: ../txt/13MP20200207_LCMv2fixidx.csv
	cp ../txt/13MP20200207_LCMv2fixidx.csv $@
