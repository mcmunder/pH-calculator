#### Measure pH with pHluorin2 ###
##
## Add description herre ##
##
####


# Set working directory -------------------------------------------------------

rm(list=ls()) # clean up
directory = "/Users/munder/Desktop/pH/MN/"
setwd(directory)     


# Load required packages ------------------------------------------------------

library("ggplot2")
library("plyr")
library("tcltk")


# Set parameters --------------------------------------------------------------
# Please adjust all of the following parameters

#options(digits=16) # change numer of digits globally

sampleName = "log"
#dirOutput = "/Users/munder/PhD/pH-measurements/"
#folderOutput = "2014_04_17_pHluorin_DNP/" # give name of output folder here


# Define functions ------------------------------------------------------------

calc_pH = function(emission_ratio, model){
    coefs = coef(model)
    #y = a + bx + cx^2 + bx^3 + ax^4
    calculated_pH = coefs[1] + (coefs[2] * emission_ratio) + (coefs[3] * emission_ratio^2) + (coefs[4] * emission_ratio^3) + (coefs[5] * emission_ratio^4)
    return(calculated_pH)
}

stats = function(df){
  mean<- mean(df)
  sd <- sd(df)
  sem <- sd / sqrt(length(df))
  data.frame(mean, sd, sem)
}


# Calibrate -------------------------------------------------------------------

source(paste(directory, "calibrate/calibrate.R", sep=""))


# Read in measurements --------------------------------------------------------

temp = tk_choose.files(caption="Please select FITC/FITC measurements")
list_ff = sapply(temp, read.delim)
list_ff = t(list_ff)
data_ff = (as.vector(list_ff[ , "Mean1"]))[[1]]

temp = tk_choose.files(caption="Please select DAPI/FITC measurements")
list_df = sapply(temp, read.delim)
list_df = t(list_df)
data_df = (as.vector(list_df[ , "Mean1"]))[[1]]

rm(list_ff, list_df)

# Calculate and normalize ratios ----------------------------------------------

ratios = data_df/data_ff

# normalize ratios to pH 7.0 from calibration experiment
norm_ratios = ratios/mean_ratio_70

# stats ratios
stats_ratios = stats(norm_ratios)


# calculate pH ----------------------------------------------------------------

pH = calc_pH(norm_ratios, best_fit)

# stats pH
stats_pH = stats(pH)

# Plot calculated pH onto calibration curve -----------------------------------

p = p + geom_point(aes(x=stats_pH$mean, y=stats_ratios$mean), size=3, colour="red")
#p = p + geom_errorbar(aes(x=stats_pH$mean, y=stats_ratios$mean,
#                          ymin=stats_ratios$mean-stats_ratios$sem,
#                          ymax=stats_ratios$mean+stats_ratios$sem), 
#                          width=.1, size=1, colour="red")
p = p + geom_line(aes(x=stats_pH$mean), size=.5, colour="red")
p = p + geom_line(aes(y=stats_ratios$mean), size=.5, colour="red")

p

ggsave(p, file=paste(sampleName, ".pdf", sep=""), width=7, height=5)


# print("---------------------------------------------------------")
# print(paste("---------------------------------------------------------",
#             "The calculated pH is:", stats_pH$mean, sep=" "))
