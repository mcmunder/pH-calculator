# Measure pH with pHluorin2 
#
# This script analyses a time series of pH measurements. For a closer 
# description see <add bitbucket repo>.


# Clean up environment

rm(list=ls()) # clean up
  
# Load required packeges

library("ggplot2")
library("plyr")
library("hash")
library("readr")
library("tcltk")
library("tcltk2")


# Please adjust all of the following parameter(s)
time_resolution =  5 # in minutes
field_of_views = 6
timepoints = 19
name = 'Energy depletion' # Used as plot title and plot filename


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

stats = function(df){
  mean<- mean(df$pH)
  sd <- sd(df$pH)
  sem <- sd / sqrt(length(df))
  data.frame(mean, sd, sem)
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

df_pH = NULL
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
  list_data_df = llply(list_files_df, read_tsv)
  list_data_ff = llply(list_files_ff, read_tsv)
  # Combine all measurements for one timepoint... So stupid! Ugliest way ever to
  # avoid a nested loop. How the hell can I slice a list in R? Or do an
  # apply(list$Mean, mean) on list_data_df etc...
  values_df = unlist(list_data_df)
  values_ff = unlist(list_data_ff)
  ratios = values_df[c(grep("Mean", names(values_df)))] / 
           values_ff[c(grep("Mean", names(values_ff)))]
  # calculate pH values;  mean and sd
  mean_pH = mean(calc_pH(ratios, best_fit))
  sd_pH = sd(calc_pH(ratios, best_fit))
  df_pH_tmp = data.frame(timepoint=(i-1)*time_resolution, mean_pH=mean_pH, sd_pH=sd_pH)
  df_pH = rbind(df_pH, df_pH_tmp)
}


# Plot

p = ggplot(df_pH, aes(x=timepoint, y=mean_pH))
p = p + geom_line(size=1, colour="red")
p = p + geom_point(size=3)
p = p + geom_errorbar(aes(ymin=mean_pH-sd_pH, ymax=mean_pH+sd_pH), size=1, width=4)
p = p + scale_x_continuous(limits=c(-2, 95), breaks=seq(0, 95, 10))
#p = p + scale_y_continuous(breaks=seq(6.6, 7.8, 0.1))
p = p + labs(title=name, x="time [min]", y="pH")
p = p + theme_bw(base_size=24)
p 

# Save plot
ggsave(p, file=paste(output_dir, "/", name, ".pdf", sep=""), width=6, height=5)



