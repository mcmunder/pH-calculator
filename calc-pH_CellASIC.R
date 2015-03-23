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

unit = "B"
dir_output = "/Users/munder/Git/pH-calculator/"
folder_output = "output/" # give name of output folder here
time_res = 10 # in minutes
timepoints = 67
field_of_views = 6

# Colour 
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


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
source(paste(directory, "/calibrate/calibrate.R", sep=""))  #/old

# loop for plotting differnt measurements together

# units = c("A", "B", "C", "D")
# 
# for(i in 1 :4){
# 
#   unit = units[i]


# Read files
setwd(paste(directory, "input_CellASIC", sep=""))

temp <- list.files(pattern=unit)
list <- sapply(temp, read.delim)
list <- t(list)
mean = sapply(list[,"Mean"], mean) 
raw_data = ldply(mean)

# tags

# field of view
fov = NULL
for(i in 1:field_of_views){
  fov_temp = rep(i, timepoints*2)
  fov = c(fov, fov_temp)
}

#channel
channel = c(rep("df", timepoints), rep("ff", timepoints))

# timestamp
time = (seq(0, timepoints-1, 1)) * time_res

dataframe = cbind(unit, channel, time, fov, raw_data)
colnames(dataframe) = c("unit", "channel", "time", "fov", "id", "counts")

# Test more substraction
#bg = 10
#dataframe$counts = dataframe$counts - bg

ratios = (subset(dataframe, channel=="df"))$counts / (subset(dataframe, channel=="ff"))$counts

# Calc pH
pH = calc.pH(ratios, best_fit)

# Final dataframe
df_pH = as.data.frame(cbind(unit, time, ratios, pH))
df_pH$pH = as.numeric(as.character(df_pH$pH))
df_pH$time = as.numeric(as.character(df_pH$time))


# Get stats over time
df_final = ddply(df_pH, .(time, unit), stats)
#assign("df_2DG_pH55", subset(df_final, time<=200))

#df_sorb_6mM$time = df_sorb_6mM$time - 300

#write.reload("df_2DG_pH70")

# 
#test = subset(df_final, time<=200)
test = subset(df_final, time>=300 & time<=500)
# sorb_6mM$time = sorb_6mM$time - 300

# Plot
p = ggplot(test, aes(x=time, y=mean))
p = p + geom_line(size=1.8, colour="blue")
p = p + geom_point(size=4)
p = p + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size=1.2, width=5)
#p = p + scale_x_continuous(limits=c(-5, 200), breaks=seq(0, 200, 60))
#p = p + scale_y_continuous(breaks=seq(6.6, 7.8, 0.1))
p = p + labs(x="time [min]", y="pH")
p = p + theme_bw(base_size=28)
p 

# Save plot
# setwd(paste(directory, "output", sep=""))
# ggsave(p, file=paste(sorb_6mM, ".pdf", sep=""), width=9, height=6)

#Plot multiple pH together

# setwd(directory)
# # # # 
#df_2DG_pH55 = read.delim("output/df_2DG_pH55.txt")
#df_2DG_pH70 = read.delim("output/df_2DG_pH70.txt")
# df_sorb1mM$unit = "C"
# df_sorb2mM = read.delim("output/sorb_2mM.txt")
# df_sorb4mM = read.delim("output/df_sorb4mM.txt")
# df_sorb4mM$unit = "B"
# df_sorb6mM = read.delim("output/df_sorb6mM.txt")
# df_sorb6mM$unit = "A"
# 
# # 
# # 
# # df_A = read.delim("output/df_A.txt")
# # df_B = read.delim("output/df_B.txt")
# # df_C = read.delim("output/df_C.txt")
# # df_D = read.delim("output/df_D.txt")
# # 
# # 
# df = rbind(df_sorb1mM, df_sorb4mM, df_sorb6mM)
# #write.reload("df")
# # df = rbind(df_A, df_B, df_D)
# # #df$mean = df$mean - 0.2
# # 
# # 
# p = ggplot(df, aes(x=time, y=mean, colour=unit))
# p = p + geom_line(size=2)
# p = p + geom_point(size=2)
# p = p + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size=1, width=5)
# p = p + scale_colour_manual(values=cbPalette)
# p = p + scale_x_continuous(limits=c(-5, 200), breaks=seq(0, 300, 60))
# p = p + scale_y_continuous(breaks=seq(5, 8, 0.5))
# p = p + labs(x="time [min]", y="pH")
# p = p + theme_bw(base_size=28)
# p 
# # 
# #save.plot(p, "pH_sorbic_acid_all.pdf", 7, 6)
