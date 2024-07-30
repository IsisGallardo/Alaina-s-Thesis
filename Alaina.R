rm(list = ls())

library(ggplot2)
library(GGally)
library(MASS)
library(caret)
library(rpart)
library(tidyverse)



data = read.csv("C:/Users/isisa/Desktop/Data/R/Alaina's data/flower_data.csv")

#delete date and comments columns
data = subset(data, select = -c(date, comments))

#Create lists of species, sites
species = unique(data$species_name)
sites = unique(data$site_id)
sites_high = c('H_H1', 'H_H2', 'H_H3', 'B_H1', 'B_H2', 'B_H3')
sites_low = c('H_L1', 'H_L2', 'H_L3', 'B_L1', 'B_L2', 'B_L3')
sites_unburned = c('H_U1', 'H_U2', 'H_U3', 'B_U1', 'B_U2', 'B_U3')

#Change DAFOR Ratings into numbers
r = 2.5
o = 10.0
f = 22.5
a = 40.0
d = 75.0
for (x in 1:length(data$DAFOR)) {
  if (is.na(data$DAFOR[x])){
    data$DAFOR[x] = 0
  } else if (data$DAFOR[x] == 'R') {
    data$DAFOR[x] = r
  } else if (data$DAFOR[x] == 'O') {
    data$DAFOR[x] = o
  } else if (data$DAFOR[x] == 'F') {
    data$DAFOR[x] = f
  } else if (data$DAFOR[x] == 'A') {
    data$DAFOR[x] = a
  } else if (data$DAFOR[x] == 'D') {
    data$DAFOR[x] = d
  }
}

#Split the data into sample periods
data_s1 = subset(data, subset = sample_period == 1)
data_s2 = subset(data, subset = sample_period == 2)
data_s3 = subset(data, subset = sample_period == 3)
data_s4 = subset(data, subset = sample_period == 4)


#Function to count species per burn site
count_species = function(df) {
  count = data.frame(species)
  for (x in 1: length(sites)){
    hold = c()
    for (i in 1: length(species)) {
      num = 0
      for (j in 1:length(df$DAFOR)) {
        if (identical(sites[x], df$site_id[j]) && identical(df$species_name[j], species[i])) {
          num = num + as.numeric(df$DAFOR[j])
        }
      }
      hold = append(hold, num)
    }
    count = cbind(count, hold)
  }
  colnames(count) = c('species', sites)
  count
}


#Output data. 
#Abundance, relative abundance and entropy.

datas1 = count_species(data_s1)


abundance <- function(site, data) {
 sum(data[[site]])
}




relative_abundance <- function(site, data) {
  nonzero <- data %>%
    select(all_of(site)) 
  total <- sum(nonzero)
  result <- mutate(nonzero, probabilities = nonzero / total)
  return(result)
}




entropy <- function(site, data) {
  P <- relative_abundance(site, data)[,2]
  result <- -sum(P * log(P),na.rm=TRUE)
  return(result)
}




#Table of estimates for abundance and Shannon-Wiener index per season.
table_of_indicies <- function(data) {
  tempAbundance <- c()
  tempShannonWiener <- c()
  for (site in sites) {
    tempAbundance = append(tempAbundance, abundance(site, data))
    tempShannonWiener = append(tempShannonWiener, entropy(site, data))
  }
  table=cbind("Sites"=sites, "Abundance"=tempAbundance, "Shannon-Wiener"=tempShannonWiener)
  return(table)
}




y <- c('data_s1', 'data_s2', 'data_s3', 'data_s4')
#Produce tables for each season
for (season in y) {
  datas <- get(season)  # Accessing the data frame using the season variable
  z <- count_species(datas)
  z <- table_of_indicies(z)
  print(z)  # Print the result for the current season
}

#species <- N(s)
#abundance <- sum(N(s),s)
#relative_abundance <- table(N(s) / abundance)
#ShannonWienerInxe <- entropy(relative_abundance)
