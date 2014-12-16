#### Measure pH with pHluorin2 ###
##
## Timeseries in CellASICS chamber##
##
####


# Set working directory -------------------------------------------------------

rm(list=ls()) # clean up
directory = "~/Git/pH-calculator/"
setwd(directory)     

# Load required packages ------------------------------------------------------

library("ggplot2")
library("plyr")
#library("tcltk")


# Set parameters --------------------------------------------------------------
# Please adjust all of the following parameters

condition = "UnitB"
dir_output = "/Users/munder/Git/pH-calculator/"
folder_output = "output/" # give name of output folder here
timeRes = 10 # in minutes
field_of_views = 6


# Define functions ------------------------------------------------------------
calc.pH = function(emission_ratio, model){
  #browser()
  coefs = coef(model)
  #y = a + bx + cx^2 + bx^3 + ax^4
  calculated_pH = coefs[1] + (coefs[2] * emission_ratio) + (coefs[3] * emission_ratio^2) + (coefs[4] * emission_ratio^3) + (coefs[5] * emission_ratio^4)
  return(calculated_pH)
}

stats = function(df){
  #browser()
  mean<- mean(df$pH)
  sd <- sd(df$pH)
  sem <- sd / sqrt(length(df))
  data.frame(mean, sd, sem)
}

# remove background
rm.bg = function(df_bg, ff_bg){
  background = seq(0, 100, length=1000)
  for(i in seq_along(background)){
    #browser()
    bg = background[i]
    cor_df_bg = df_bg - bg
    cor_ff_bg = ff_bg - bg
    ratio = cor_df_bg / cor_ff_bg
    mean_ratio = mean(ratio)
    pH = as.numeric(calc.pH(mean_ratio, best_fit))
    if(round(pH, digits=3)==7.400){
      #print(bg)
      #print(pH)
      return(bg)
    }
  }
}

write.reload = function(filename){
  # write
  write.table(get(filename), 
              paste(dir_output, folder_output, filename, ".txt", sep=""), sep="\t")
  # reload
  assign(filename, 
         read.delim(paste(dir_output, folder_output, filename, ".txt", sep="")))  
}

save.plot = function(plot, file_name, width, height){
  setwd(paste(dir_output, folder_output, sep=""))
  ggsave(p, file=file_name, width=width, height=height)
  setwd(directory)
}

# Calibrate
source(paste(directory, "calibrate/calibrate.R", sep=""))

# Read files in a for loop
setwd(paste(directory, "input_CELLASIC", sep=""))

# FITC/FITC
df_ff = NULL
for(i in 1:field_of_views){
  #FITC/FITC
  df_ff_temp = read.delim(paste("Results_", condition, "_0", i, "_ff", ".txt", sep=""))
  time = df_ff_temp$X*timeRes-10
  replicate = rep(i, nrow(df_ff_temp))
  channel = rep("ff", nrow(df_ff_temp))
  df_ff_temp = cbind(df_ff_temp, time, replicate, channel)
  df_ff = rbind(df_ff, df_ff_temp)
} 

# DAPI/FITC
df_df = NULL
for(i in 1:field_of_views){
  #FITC/FITC
  df_df_temp = read.delim(paste("Results_", condition, "_0", i, "_df", ".txt", sep=""))
  time = df_ff_temp$X*timeRes-10
  replicate = rep(i, nrow(df_df_temp))
  channel = rep("df", nrow(df_df_temp))
  df_df_temp = cbind(df_df_temp, time, replicate, channel)
  df_df = rbind(df_df, df_df_temp)
} 

# Correct for autofluorescence background (this assumes that intracellular pH is at pH 7.4 
# in log phase cells and uses this asumption to substract background from both channels

# df_bg = ((ddply(df_df, .(time), summarize, mean=mean(Mean1)))$mean)[1:6]
# ff_bg = ((ddply(df_ff, .(time), summarize, mean=mean(Mean1)))$mean)[1:6]
# 
# bg = rm.bg(df_bg, ff_bg)

# bg = 50
#   
# df_df$Mean1 = df_df$Mean1 - bg
# df_ff$Mean1 = df_ff$Mean1 - bg


# Ratios
df_ratios = data.frame(cbind(df_ff$time, df_ff$replicate, df_df$Mean1/df_ff$Mean1))
colnames(df_ratios) = c("time", "replicate", "ratio")

# Normalize
#norm_ratios = df_ratios$ratio/mean_ratio_70

# Calc pH
pH = calc.pH(df_ratios$ratio, best_fit)

# Final dataframe
df_pH = cbind(df_ratios, pH, rep(condition, length(pH)))
colnames(df_pH) = c("time", "replicate", "ratio", "pH", "condition")

# Get stats over time
df_final = ddply(df_pH, .(time, condition), stats)
assign(paste("df_", condition, sep=""), df_final)
write.reload(paste("df_", condition, sep=""))

# Plot


p = ggplot(df_final, aes(x=time, y=mean))
p = p + geom_line(size=1.8, colour="blue")
p = p + geom_point(size=4)
p = p + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size=1.2, width=5)
#p = p + scale_x_continuous(limits=c(-5, 300), breaks=seq(0, 300, 60))
#p = p + scale_y_continuous(breaks=seq(6.6, 7.8, 0.1))
p = p + labs(x="time [min]", y="pH")
p = p + theme_bw(base_size=28)
p 

# Save plot
setwd(paste(directory, "output", sep=""))
ggsave(p, file=paste(condition, ".pdf", sep=""), width=9, height=6)

# # Plot multiple pH together
# 
# setwd(directory)
# 
# df_UnitA = read.delim("output/df_UnitA.txt")
# df_UnitB = read.delim("output/df_UnitB.txt")
# df_UnitC = read.delim("output/df_UnitC.txt")
# df_UnitD = read.delim("output/df_UnitD.txt")
# 
# df = rbind(df_UnitC, df_UnitD)
# 
# p = ggplot(df, aes(x=time, y=mean, colour=condition))
# p = p + geom_line(size=2)
# p = p + geom_point(size=2)
# p = p + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size=1, width=5)
# p = p + scale_x_continuous(limits=c(-5, 300), breaks=seq(0, 300, 60))
# #p = p + scale_y_continuous(breaks=seq(6.6, 7.8, 0.1))
# p = p + labs(x="time [min]", y="pH")
# p = p + theme_bw(base_size=28)
# p 
# 
# # save.plot(p, "pH_2DG_AntiA.pdf", 9, 6)







