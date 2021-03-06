# Measure pH with pHluorin2 
#
# This script analyses a time series of pH measurements. For a closer 
# description see <add bitbucket repo>.


# Clean up environment

rm(list=ls()) # clean up
  
# Load required packeges

library("ggplot2")
library("plyr")
library("readr")
library("tcltk")
library("tcltk2")
library("gdata")

# Please adjust all of the following parameter(s)
time_resolution =  10 # in minutes
field_of_views = 6
timepoints = 3
name = 'Heatshock_shaker' # Used as plot title and plot filename


# Set directories ---------------------------------------------------------

working_dir = tk_choose.dir(caption=("Choose working directory!"))
input_dir_measurements = tk_choose.dir(caption=("Select input directory containing measurement results."))
input_dir_calibration = tk_choose.dir(caption=("Select directory containing calibration data."))
output_dir = tk_choose.dir(caption=("Select output directory. Dataframes and plots will be saved here."))


# Define functions ------------------------------------------------------------
calc_pH = function(emission_ratio, model){
  coefs = coef(model)
  #y = a + bx + cx^2 + dx^3 + ex^4
  calculated_pH = coefs[1] + 
                  (coefs[2] * emission_ratio) + 
                  (coefs[3] * emission_ratio^2) + 
                  (coefs[4] * emission_ratio^3) + 
                  (coefs[5] * emission_ratio^4)
  return(calculated_pH)
}


# Calibrate ---------------------------------------------------------------

setwd(working_dir)
source("calibrate.R")


# Measure -----------------------------------------------------------------

setwd(input_dir_measurements)

# Loop throgh timepoints, calculate dapi/fitc to fifc/fitc ratios, 
# concatenate ratios for each timepoint and compute mean and sd of
# ratios. If segmentation in Fiji would be perfect (each ROI=1 cell), 
# the mean would be the mean ratio per cell, sd would show how much 
# the pH varies between cells etc. Right now segmentation is NOT perfect!

df_all = NULL
df_mean = NULL
for (i in 1:timepoints) {
  # Generate file list for fitc/fitc and dapi/fitc channels for 
  # one timepoint:
  if (i < 10) {
    list_files_df = list.files(pattern = paste('_df_', '0', i, sep=''))
    list_files_ff = list.files(pattern = paste('_ff_', '0', i, sep=''))
  } else {
    list_files_df = list.files(pattern = paste('_df_', i, sep=''))
    list_files_ff = list.files(pattern = paste('_ff_', i, sep=''))
  }
  # Read in files with Hadley's readr function read_tsv
  #browser()
  list_data_df = llply(list_files_df, read_tsv)
  list_data_ff = llply(list_files_ff, read_tsv)
  # Combine all measurements for one timepoint... So stupid! Ugliest way ever to
  # avoid a nested loop. How the hell can I slice a list in R? Or do an
  # apply(list$Mean, mean) on list_data_df etc...
  values_df = unlist(list_data_df)
  values_ff = unlist(list_data_ff)
  intensities_df = values_df[c(grep("Mean", names(values_df)))]
  intensities_ff = values_ff[c(grep("Mean", names(values_ff)))]
  ratios = intensities_df / intensities_ff
  # calculate pH values;  mean and sd
  pH_all = calc_pH(ratios, best_fit)
  mean_pH = mean(calc_pH(ratios, best_fit))
  sd_pH = sd(calc_pH(ratios, best_fit))
  df_all_temp = data.frame(timepoint=(i-1)*time_resolution, ratios=ratios, pH=pH_all)
  df_mean_tmp = data.frame(timepoint=(i-1)*time_resolution, mean_df = mean(intensities_df),
                         mean_ff = mean(intensities_ff), mean_pH=mean_pH, sd_pH=sd_pH)
  df_mean = rbind(df_mean, df_mean_tmp)
  df_all = rbind(df_all, df_all_temp)
}


# Plot mean

p = ggplot(df_mean, aes(x=timepoint, y=mean_pH))
p = p + geom_line(size=1, colour="red")
p = p + geom_point(size=3)
p = p + geom_errorbar(aes(ymin=mean_pH-sd_pH, ymax=mean_pH+sd_pH), size=1, width=4)
#p = p + scale_x_continuous(limits=c(-2, 125), breaks=seq(0, 125, 10))
#p = p + scale_y_continuous(breaks=seq(5.2, 7.9, 0.2))
p = p + labs(title=name, x="time [min]", y="pH")
p = p + theme_bw(base_size=24)
p 

# Save plot and dataframe
setwd(output_dir)
ggsave(p, file=paste(name, ".pdf", sep=""), width=6, height=5)
write_csv(df_mean, paste('df_mean_', name, '.csv', sep=''))
setwd(working_dir)


# Plot all

p = ggplot(df_all, aes(x=factor(timepoint), y=pH))
p = p + geom_boxplot()
p = p + labs(title=name, x="time [min]", y="pH")
p = p + theme_bw(base_size=24)
p 

# Save plot and dataframe
setwd(output_dir)
ggsave(p, file=paste(name, ".pdf", sep=""), width=6, height=5)
write_csv(df_all, paste('df_all_', name, '.csv', sep=''))
setwd(working_dir)


# # Plot stuff together
# name1 = "Energy-depletion_pH55"
# name2 = "Energy-depletion_pH70"
# 
# setwd(output_dir)
# pH55 = read_csv(paste('df_pH_', name1, '.csv', sep=''))
# pH70 = read_csv(paste('df_pH_', name2, '.csv', sep=''))
# 
# final = combine(pH55, pH70)
# 
# p = ggplot(final, aes(x=timepoint, y=mean_pH, color=factor(source)))
# p = p + geom_line(size=1)
# p = p + geom_point(size=3)
# p = p + geom_errorbar(aes(ymin=mean_pH-sd_pH, ymax=mean_pH+sd_pH), size=1, width=4)
# #p = p + scale_x_continuous(limits=c(-2, 125), breaks=seq(0, 125, 10))
# #p = p + scale_y_continuous(breaks=seq(5.2, 7.9, 0.2))
# p = p + labs(title="Energy-depletion_combined", x="time [min]", y="pH")
# p = p + theme_bw(base_size=24)
# p 
# 
# ggsave(p, file=paste("Eenrgy-depletion_combined", ".pdf", sep=""), width=8, height=5)
# 
# 
# 
# 
