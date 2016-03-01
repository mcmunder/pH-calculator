# This script generates a pH calibration curve 

# The ollowing section is for testig the calibrate script seperately.
# Comment out if calling as source() from calc-pH.R 

# rm(list=ls()) # clean up
# library("ggplot2")
# library("plyr")
# library("tcltk2")
# library("readr")
# input_dir_calibration = tk_choose.dir(caption=("Select input directory containing calibration data."))


# List of files in input directory

setwd(input_dir_calibration)
list_files = list.files(pattern = '.txt')
list_dataframes = llply(list_files, read_tsv)

list_norm_ratios = NULL
for(i in seq(2, length(list_dataframes), 2)){
  # Area should always normalize to 1! Good control! Just skip the [3] indexing
  # normalize to tritc/tritc channel
  #norm_ratio = (list_dataframes[[i-2]][3] / list_dataframes[[i]][3]) / (list_dataframes[[i-1]][3] / list_dataframes[[i]][3])
  norm_ratio = list_dataframes[[i-1]][3]  / list_dataframes[[i]][3] 
  list_norm_ratios = append(list_norm_ratios, norm_ratio)
  rm(i, norm_ratio)
}


# 8 pH conditions as numeric values

pH_values = c(4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0)


# compute mean and sd of ratio over 6 field of views

df_calibrate = NULL
for(j in seq(6, length(list_norm_ratios), 6)){
  pH = pH_values[j/6]
  mean_ratio = mean(unlist(list_norm_ratios[(j-5):j]))
  sd_ratio = sd(unlist(list_norm_ratios[(j-5):j]))
  sem_ratio = sd_ratio / sqrt(length(unlist(list_norm_ratios[(j-5):j])))
  df_calibrate_temp = data.frame(pH, mean_ratio, sd_ratio, sem_ratio)
  df_calibrate = rbind(df_calibrate, df_calibrate_temp)
  rm(j, mean_ratio, sd_ratio, sem_ratio, pH)
}


# plot

p = ggplot(df_calibrate, aes(x=pH, y=mean_ratio))
p = p + geom_point(size=4)
p = p + geom_errorbar(aes(ymin=mean_ratio-sd_ratio, ymax=mean_ratio+sd_ratio), width=.1, size=1)
p = p + scale_x_continuous(breaks=seq(5, 8, .5))
#p = p + scale_y_continuous(limits=c(0.5, 1.3), breaks=seq(0.5, 1.3, .1))
#p = p + geom_smooth(method="lm", formula = y~poly(x,2), colour='green') 
#p = p + geom_smooth(method="lm", formula = y~poly(x,3), colour='red')
p = p + geom_smooth(method="lm", formula = y~poly(x,4), colour='blue') # "best" fit 
p = p + labs(x="pH", y="Emission ratio (AU)")
p = p + theme_bw(base_size=24)
print(p) 


# Safe best for pH calculation (calc_pH.R and calc-pH_time.R)
x = df_calibrate$pH
y = df_calibrate$mean_ratio

# Swop x and y for computation of oH in the main script
best_fit = lm(x~poly(y,4,raw=TRUE))


