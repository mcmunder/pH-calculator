# This script generates a pH calibration curve based on measurements
# taken on September 3, 2013. For these measurements cells were treated 
# exactly as descibed in Brett et al. 2005, (DOI: 10.1091/mbc.E04-11-0999),
# Figure 2B.

setwd(paste(directory, "calibrate", sep=""))

temp <- list.files(pattern="df")
list_df <- sapply(temp, read.delim)
list_df <- t(list_df )
data_df <- data.frame(list_df[,"Mean1"])

temp <- list.files(pattern="ff")
list_ff <- sapply(temp, read.delim)
list_ff <- t(list_ff )
data_ff <- data.frame(list_ff[,"Mean1"])

# ratios
ratios <- data_df/data_ff
mean_ratios <- sapply(ratios,mean)

# normalize ratios to pH7
mean_ratio_70 = mean_ratios[[5]]
norm_mean_ratios = mean_ratios/mean_ratio_70

# sd and sem of normalized ratios
sd_ratio<- sapply(ratios,sd)
sem_ratio <- sd_ratio/sqrt(length(data_df))

# create final dataframe
pH = c(5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0)
dataframe = data.frame(cbind(pH, mean_ratios, sd_ratio, sem_ratio))
colnames(dataframe) = c('pH', 'Mean', 'SD','SEM')

# curve fitting 

# what's the "best" fit?
x = dataframe$Mean
y = dataframe$pH

plot(x,y,pch=19)

#fit first degree polynomial equation:
fit  <- lm(y~x)
#second degree
fit2 <- lm(y~poly(x,2,raw=TRUE))
#third degree
fit3 <- lm(y~poly(x,3,raw=TRUE))
#fourth degree
fit4 <- lm(y~poly(x,4,raw=TRUE))
#generate range of 100 numbers between highers and lowest dataframe$mean
xx <- seq(min(x), max(x), length=1000)
plot(x,y,pch=19,ylim=c(4.8, 8.2))
lines(xx, predict(fit, data.frame(x=xx)), col="red")
lines(xx, predict(fit2, data.frame(x=xx)), col="green")
lines(xx, predict(fit3, data.frame(x=xx)), col="blue")
lines(xx, predict(fit4, data.frame(x=xx)), col="purple") # seems to be "best" here

best_fit = lm(y~poly(x,4,raw=TRUE))


# plot
p = ggplot(dataframe, aes(x=pH, y=Mean))
p = p + geom_point(size=4)
p = p + geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.1, size=1)
p = p + scale_x_continuous(breaks=seq(5, 8, .5))
#p = p + scale_y_continuous(limits=c(0.5, 1.3), breaks=seq(0.5, 1.3, .1))
p = p + geom_smooth(method="lm", formula = y~poly(x,4), colour='blue') # "best" fit
p = p + labs(x="pH", y="normalized emission ratio")
p = p + theme_bw(base_size=24)
print(p) 

setwd(directory)

rm(list=setdiff(ls(), c("directory", "folderOutput", "calc_pH", "stats", "best_fit", "name", "timeRes", "mean_ratio_70", "p")))

