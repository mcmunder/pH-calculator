#### Measure pH with pHluorin2 ###
##
## Timeseries in CellASICS chamber##
##
####


# Clean up --------------------------------------------------------------------
rm(list=ls()) # clean up
  

# Load required packages ------------------------------------------------------

library("ggplot2")
library("plyr")
library("tcltk2")
library("gdata")


# Set parameters --------------------------------------------------------------
# Please adjust all of the following parameters

name = "pH56" # Must match the name of result files
timeRes = 10 # in minutes


# Set input and output directory ------------------------------------------

input_dir_measurements = tk_choose.dir(caption=("Select input directory containing measurement results."))
input_dir_calibration = tk_choose.dir(caption=("Select input directory containing calibration data."))
output_dir = tk_choose.dir(caption=("Select output directory. Dataframes and plots will be saved here."))


# Define functions ------------------------------------------------------------
calc_pH = function(emission_ratio, model){
  coefs = coef(model)
  #y = a + bx + cx^2 + bx^3 + ax^4
  calculated_pH = coefs[1] + (coefs[2] * emission_ratio) + (coefs[3] * emission_ratio^2) + (coefs[4] * emission_ratio^3) + (coefs[5] * emission_ratio^4)
  return(calculated_pH)
}

stats = function(df){
  mean<- mean(df$pH)
  sd <- sd(df$pH)
  sem <- sd / sqrt(length(df))
  data.frame(mean, sd, sem)
}

# Calibrate
source(paste(input_dir_calibration, "/calibrate.R", sep=""))

# Read files in a for loop
setwd(input_dir_measurements)

# FITC/FITC
df_ff = NULL
for(i in 1:6){
  #FITC/FITC
  df_ff_temp = read.delim(paste("Results_", name, "_0", i, "_ff", ".txt", sep=""))
  time = df_ff_temp$X * timeRes - 10
  replicate = rep(i, nrow(df_ff_temp))
  channel = rep("ff", nrow(df_ff_temp))
  df_ff_temp = cbind(df_ff_temp, time, replicate, channel)
  df_ff = rbind(df_ff, df_ff_temp)
} 

# DAPI/FITC
df_df = NULL
for(i in 1:6){
  #FITC/FITC
  df_df_temp = read.delim(paste("Results_", name, "_0", i, "_df", ".txt", sep=""))
  time = df_ff_temp$X * timeRes - 10
  replicate = rep(i, nrow(df_df_temp))
  channel = rep("df", nrow(df_df_temp))
  df_df_temp = cbind(df_df_temp, time, replicate, channel)
  df_df = rbind(df_df, df_df_temp)
} 

# Ratios
df_ratios = data.frame(cbind(df_ff$time, df_ff$replicate, df_df$Mean1/df_ff$Mean1))
colnames(df_ratios) = c("time", "replicate", "ratio")

# Normalize
#norm_ratios = df_ratios$ratio/mean_ratio_70

# Calc pH
pH = calc_pH(df_ratios$ratio, best_fit)

# Final dataframe
df_pH = cbind(df_ratios, pH)

# Get stats over time
df_final = ddply(df_pH, .(time), stats)

# Plot

p = ggplot(df_final, aes(x=time, y=mean))
p = p + geom_line(size=1, colour="red")
p = p + geom_point(size=3)
p = p + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size=1, width=4)
#p = p + scale_x_continuous(limits=c(0, 420), breaks=seq(0, 420, 60))
#p = p + scale_y_continuous(breaks=seq(6.6, 7.8, 0.1))
p = p + labs(title=name, x="time [min]", y="pH")
p = p + theme_bw(base_size=24)
p 

# Save plot
ggsave(p, file=paste(output_dir, "/", name, ".pdf", sep=""), width=9, height=6)



