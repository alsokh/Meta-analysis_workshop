# to set the working directory
setwd("~/Documents/Main/Meta_analysis")


# to load needed packages
library(tidyverse)
library(meta)
library(metafor)
library(dmetar)
library(dmetatools)

# loading sample data for meta-analyses
data("asthma")
data("woodyplants")
data("Olkin1995")
data("Fleiss93")
data("dat.bcg")
data("dat.konstantopoulos2011")
data("dat.begg1989")
data("dat.assink2016")
data("dat.baker2009")
data("dat.viechtbauer2021")
data("BdiScores")
data("DepressionMortality")
data("HealthWellbeing")
data("ThirdWave")

# to assign data frames to related variable
df_mean <- BdiScores
df_percentage <- dat.baker2009
df_SMD <- SuicidePrevention
df_OR <- DepressionMortality
df_Cor <- HealthWellbeing
df_three_level <- dat.assink2016 # this model is not being used usually
df_bivariate <- asthma # this model is not being used usually
df_gen <- ThirdWave


# save the smd data into csv
write_csv(x = df_SMD, file = "smd.csv")

# metaprop --> for single proportion
# metabin --> for more than one group of proportions
# metamean --> for single mean
# metacont --> for SMD effect in meta-analysis
# metacor --> for correlation coefficient meta-analysis
# metagen --> a general package for conducting meta-analysis using caculated effect



# to look at the mean data
head(df_mean)
glimpse(df_mean)

# to run the model for single mean
model_mean <- metamean(n = n, mean = mean, sd = sd, studlab = author, 
                       data = df_mean)
# check the model
summary(model_mean)
model_mean


# to run the model for single proportion
model_percent <- metaprop(event = exac, n = total, studlab = study, 
                          data = df_percentage)
# check the model
model_percent
summary(model_percent)


# to run the model for continuous data
model_smd <- metacont(n.e = n.e, mean.e = mean.e, sd.e = sd.e, 
                      n.c = n.c, mean.c = mean.c, sd.c = sd.c,
                      studlab = author, data = df_SMD)

# check the model
model_smd
summary(model_smd)

head(df_OR)
# to run the model for OR
model_OR <- metabin(event.e = event.e, n.e = n.e,
                    event.c = event.c, n.c = n.c, 
                    studlab = author, data = df_OR)

# check the model
model_OR
summary(model_OR)


head(df_Cor)
# to run the model for correlation data
model_cor <- metacor(cor = cor, n = n, studlab = author, 
                             data = df_Cor)

# check the model
model_cor
summary(model_cor)


head(ThirdWave)
# to run the model with the effect size calculated
model_gen <- metagen(TE = TE, seTE = seTE, studlab = Author, data = ThirdWave)

# check the model
model_gen
summary(model_gen)



# loading needed Packages needed 
library(meta)
library(metafor)
library(rmeta)
library(netmeta)
library(mada)
library(dmetar)
library(metasens)


effect <- escalc(measure = "SMD", m1i = mean.e, m2i = mean.c, 
                 sd1i = sd.e, sd2i = sd.c, 
                 n1i = n.e, n2i = n.c, data = df_SMD)

head(effect)



#let's do meta-analysis with meta
df_meta <- readr::read_csv("smd.csv")
head(df_meta)

meta_model_1 <- metacont(n.e = n.e, mean.e = mean.e, sd.e = sd.e, 
                         n.c = n.c, mean.c = mean.c, sd.c = sd.c,
                         studlab = author, data = df_meta)

summary(meta_model_1)

# funnel plot
meta::forest(meta_model_1)
meta::funnel(meta_model_1)

# Publication bias
meta::metabias(meta_model_1)
meta::metabias(meta_model_1, k.min = 3)

# sensitivity analysis, outliers, leave-one-out, etc
outliers <- dmetar::find.outliers(meta_model_1)
sensitivity <- metainf(meta_model_1)
summary(sensitivity)
# leave_one_out <- leave1out(meta_model_2)
influence <- InfluenceAnalysis(meta_model_1)
summary(influence)
plot(influence)

trim <- trimfill(meta_model_1)
funnel(trim)

# gosh <- gosh.rma(meta_model_2)
# a <- dmetar::gosh.diagnostics(gosh)
# plot(a)

plot(meta_model_1)
# let's do meta-analysis with metafor
meta_model_2 <- rma(effect)

gosh <- gosh.rma(meta_model_2)
gosh <- gosh.rma(meta_model_2)
gosh <- gosh.rma(meta_model_2)
