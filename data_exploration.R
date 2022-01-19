# Set Working Directory (where data and R script are located)
wd <- "~/Documents/Artificial Intelligence/Master/2122 Sem. 1/Bayesian Networks/Assignment1/Assignment_1"
setwd(wd)

# Set Seed
set.seed(123)

# Import Packages
library('car')

# Import Data
d <- read.csv("forestfires.csv", colClasses=c("integer","integer","factor","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

# Show head of the Data
head(d)

# Our data has a few categorical variables.
# These are: month and day.
# First let's fix the month variable, and change it into a numeric variable
d$month <- ordered( d$month, levels = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 
                                        'jul', 'aug', 'sep', 'oct', 'nov', 'dec'))
d_new <- unclass(d$month)
d = subset(d, select = -c(month))
d$month <- as.numeric(d_new)

# Check if day is informative.
days <- c('mon', 'tue', 'wed', 'thu', 'fri', 'sat', 'sun')
d$day <- ordered(d$day, levels = days)

png(file="~/Documents/Artificial Intelligence/Master/2122 Sem. 1/Bayesian Networks/Assignment1/Assignment_1/plts/days.png",
    width=600, height=350)
par(las=2)
barplot(table(d$day), main = "Amount of Datapoints per Day", xlab = "Day", ylab = "Count")
dev.off()

# We don't want to use day, so drop it with X and Y.
d <- subset(d, select = -c(X, Y, day))

# And we now have a dataset with the variables we want to use.
head(d)


## Exploring the variables in more detail

# Density of the Area variable:
plot(density(d$area), main = 'Density of Area')
head(sort(d$area, decreasing = TRUE))
qqPlot(d$area)

# We can tell this is heavily skewed to the right, but that this is not necessarily caused by outliers.
# This skewness can be fixed by doing a log transform of this variable
d$area <- log(d$area + 1)
plot(density(d$area), main = 'Density of log-transformed Area')
outliers = qqPlot(d$area, main = "QQ-plot of Area", ylab = "Area (in ha)")
d <- d[-outliers,]
plot(density(d$area), main = 'Density of log-transformed Area after outlier removal')
qqPlot(d$area)

# Density of the rain variable:
plot(density(d$rain), main = 'Density of Rain')
head(sort(d$rain, decreasing = TRUE))
outliers = qqPlot(d$rain, main = "QQ-plot of Rain", ylab = "Rain (in mm/m2)")

# We can tell this is heavily skewed to the right again, but after inspecting the values, this definitaly is an outlier:
# So let's remove this datapoint from the dataset.
d <- d[-outliers[1],]

plot(density(d$rain), main = 'Density of Rain after Outlier Removal')
qqPlot(d$rain, main = "QQ-plot of Rain", ylab = "Rain (in mm/m2)")

# Density of the wind variable:
plot(density(d$wind), main = 'Density of Wind')
outliers = qqPlot(d$wind, main = "QQ-plot of Wind", ylab = "Wind (in km/h)")
# This looks reasonable well to fit a Gaussian to.

# Density of the RH variable:
plot(density(d$RH), main = 'Density of Relative Humidity')
# This also looks reasonable well to fit a Gaussian to.
qqPlot(d$RH)

# Density of the temperature variable:
plot(density(d$temp), main = 'Density of Temperature')
qqPlot(d$temp, main = "QQ-plot of Temperature")
# This also looks reasonable well to fit a Gaussian to.


# Density of the Initial Spread Index variable:
plot(density(d$ISI), main = 'Density of Initial Spread Index')
head(sort(d$ISI, decreasing = TRUE))
outliers = qqPlot(d$ISI, main = 'QQ-plot of Initial Spread Index')
# We see again a heavily right-skewed distribution, but after inspecting the ISI values, there is one that is extremely large compared to others.
# So let's remove than one as well.
d <- d[-outliers[1],]
plot(density(d$ISI), main = 'Density of Initial Spread Index after Outlier Removal')
outliers = qqPlot(d$ISI, main = 'QQ-plot of Initial Spread Index after Outlier Removal')

# Density of the Drought Code variable:
plot(density(d$DC), main = 'Density of Drought Code')
outliers = qqPlot(d$DC, main = 'QQ-plot of Drought Code')

# This is a problematic distribution, as it looks like it's bimodal.
# However, there is nothing we can do about this.


# Density of the Duff Moisture Code variable:
plot(density(d$DMC), main = 'Density of Duff Moisture Code')
outliers = qqPlot(d$DMC, main = 'QQ-plot of Duff Moisture Code')
# Again, looks bimodal.


# Density of the Fine Fuel Moisture Code variable:
png(file="~/Documents/Artificial Intelligence/Master/2122 Sem. 1/Bayesian Networks/Assignment1/Assignment_1/plts/density_ffmc.png", width=600, height=350)
plot(density(d$FFMC), main = 'Density of FFMC')
dev.off()

head(sort(d$FFMC, decreasing = FALSE))
outliers = qqPlot(d$FFMC, main = 'QQ-plot of Fine Fuel Moisture Code')

# Now we have a left-skewed distribution, with an outlier (18.7).
# Let's remove this as well.
d <- d[-outliers,]
plot(density(d$FFMC), main = 'Density of Fine Fuel Moisture Code after Outlier Removal')
outliers = qqPlot(d$FFMC)
d <- d[-outliers[1],]
plot(density(d$FFMC), main = 'Density of Fine Fuel Moisture Code after Outlier Removal')
qqPlot(d$FFMC)

d$FFMC <- log(max(d$FFMC)-d$FFMC+1)

png(file="~/Documents/Artificial Intelligence/Master/2122 Sem. 1/Bayesian Networks/Assignment1/Assignment_1/plts/density_processed_ffmc.png", width=600, height=350)
plot(density(d$FFMC), main = 'Density of FFMC after Preprocessing')
dev.off()
qqPlot(d$FFMC)

## Now we have removed the outliers using QQ-plots, for the following variables:

# Area - 2 outliers removed, those with the highest area were substantially larger than other values.
# Rain - 1 outlier removed, the one with the highest value (6.4) which was again substantially larger than other values
# ISI - 1 outlier removed, the one with the highest value (around 60), idem as before
# FFMC - 3 outliers removed, the three lowest values. 

# let's save this dataset to a new .csv file
write.csv(d, paste(wd,"/explored_forestfires.csv", sep=""), row.names = FALSE)

