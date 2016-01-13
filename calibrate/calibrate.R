# This script generates a pH calibration curve based on measurements
# taken on September 3, 2013. For these measurements cells were treated 
# exactly as descibed in Brett et al. 2005, (DOI: 10.1091/mbc.E04-11-0999),
# Figure 2B.

setwd(input_dir_calibration)

# 8 pH conditions (values)
pH_values = c(4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0)





# List of files in input directory
list_files = list.files(pattern = '.txt')
list_dataframes = llply(list_files, read.delim)

list_norm_ratios = NULL
for(i in seq(2, length(list_dataframes), 2)){
  # Area should always normalize to 1! Good control! Just skip the [3] indexing
  # normalize to tritc/tritc channel
  #norm_ratio = (list_dataframes[[i-2]][3] / list_dataframes[[i]][3]) / (list_dataframes[[i-1]][3] / list_dataframes[[i]][3])
  norm_ratio = list_dataframes[[i-1]][3]  / list_dataframes[[i]][3] 
  list_norm_ratios = append(list_norm_ratios, norm_ratio)
  rm(i, norm_ratio)
}

df_stats = NULL
for(j in seq(6, length(list_norm_ratios), 6)){
  pH = pH_values[j/6]
  mean_ratio = mean(unlist(list_norm_ratios[(j-5):j]))
  sd_ratio = sd(unlist(list_norm_ratios[(j-5):j]))
  sem_ratio = sd_ratio / sqrt(length(unlist(list_norm_ratios[(j-5):j])))
  df_stats_temp = data.frame(pH, mean_ratio, sd_ratio, sem_ratio)
  df_stats = rbind(df_stats, df_stats_temp)
  rm(j, mean_ratio, sd_ratio, sem_ratio, pH)
}

# curve fitting 
p = ggplot(df_stats, aes(x=pH, y=mean_ratio))
p = p + geom_point()
p


# what's the "best" fit?
x = df_stats$pH
y = df_stats$mean_ratio

plot(x,y,pch=19)

#fit first degree polynomial equation:
fit  <- lm(y~x)
#second degree
fit2 <- lm(y~poly(x,2,raw=TRUE))
#third degree
fit3 <- lm(y~poly(x,3,raw=TRUE))
#fourth degree
fit4 <- lm(y~poly(x,4,raw=TRUE))
#generate range of 1000 numbers between highers and lowest dataframe$mean
xx <- seq(min(x), max(x), length=1000)
plot(x,y,pch=19,ylim=c(4.5, 8.0))
lines(xx, predict(fit, data.frame(x=xx)), col="red")
lines(xx, predict(fit2, data.frame(x=xx)), col="green")
lines(xx, predict(fit3, data.frame(x=xx)), col="blue")
lines(xx, predict(fit4, data.frame(x=xx)), col="purple") # seems to be "best" here

best_fit = lm(y~poly(x,4,raw=TRUE))

# plot
p = ggplot(df_stats, aes(x=pH, y=mean_ratio))
p = p + geom_point(size=4)
p = p + geom_errorbar(aes(ymin=mean_ratio-sem_ratio, ymax=mean_ratio+sem_ratio), width=.1, size=1)
p = p + scale_x_continuous(breaks=seq(5, 8, .5))
#p = p + scale_y_continuous(limits=c(0.5, 1.3), breaks=seq(0.5, 1.3, .1))
#p = p + geom_smooth(method="lm", formula = y~poly(x,2), colour='green') 
p = p + geom_smooth(method="lm", formula = y~poly(x,3), colour='red')
p = p + geom_smooth(method="lm", formula = y~poly(x,4), colour='blue') # "best" fit 
p = p + labs(x="pH", y="Emission ratio (AU)")
p = p + theme_bw(base_size=24)
print(p) 

