rm(list = ls())

#### Data import  and checking ####
df <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZyFamilyReadsPerKbp.csv")
rft <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZymeResponse_p_FPR.csv")

df[1:6]

library(roperators)
library(vegan)
library(data.table)
library(ggplot2)
library(Hmisc)
library(doBy)

#### Data Cleaning and checking ####

# Normalise per Gbp
df[, 4:length(df)] <- df[, 4:length(df)] / (df$Reads/1e9) 

# Check the data
head(df[,1:10])
hist(unlist(df[, 4:length(df)]), breaks = 100) # Range of read counts
which(colSums(df[4:length(df)]) == 0) # Any fully zero CAZy columns?
head(rft)

#### Were the CAZymes thoroughly sampled?

sp2 <- specaccum(df[,4:length(df)], "random")
sp2
summary(sp2)
plot(sp2, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sp2, col="yellow", add=TRUE, pch="+")

raremax <- min(rowSums(df[,4:length(df)]))
raremax
rarecurve(round(df[,4:length(df)]), step = 20, sample = raremax)

# Select the significantly less rich CAZymes
SLR <- rft[rft$p_adj < 0.1 & rft$log2FC_raw < -1,]

# SSelect the significantly more rich CAZymes
SMR <- rft[rft$p_adj < 0.1 & rft$log2FC_raw > 1,]

head(SLR)
head(SMR)

#### How abundant were different CAZy types? ####

ctf <- data.table(Sample = df$Sample
           , Treatment = df$Treatment
           , GH = rowSums(df[, grep("GH", names(df))])
           , AA = rowSums(df[, grep("AA", names(df))])
           , CBM = rowSums(df[, grep("CBM", names(df))])
           )

ctf[, list(Mean = mean(GH), SD = sd(GH)), by = Treatment]
ctf[, list(Mean = mean(AA), SD = sd(AA)), by = Treatment]
ctf[, list(Mean = mean(CBM), SD = sd(CBM)), by = Treatment]

mean_ci <- function(x){
  mean_se(x, mult = 1.96)
}
ggplot(melt(ctf, id = c("Sample", "Treatment"))
       , aes(x = Treatment, y = value)) +
  geom_point() +
  facet_wrap(~variable, scales = "free") +
  stat_summary(fun.data = "mean_ci")

ctf
sum(ctf$GH)
sum(ctf$AA)
sum(ctf$CBM)

#### Was total CAZyme abundance different between treatments? ####

head(df[1:5])

# Create a total CAZyme column
tf <- df
tf$CAZymes <- rowSums(df[,4:length(df)])

# Plot the data
plot(CAZymes ~ Treatment, data = tf)
points(CAZymes ~ Treatment, data = tf, pch = 16, col = 'grey30')

# Create a model to test if there are differences
mod1 <- lm((CAZymes) ~ Treatment, data = tf)
par(mfrow=c(2,2));plot(mod1);par(mfrow=c(1,1))
car::qqPlot(residuals(mod1, type = 'pearson'))

# Was treatment important?
drop1(mod1, test = 'F')
anova(mod1)
summary(mod1)
# Perform multiple comparisons test
library(multcomp)
summary(glht(mod1, linfct = mcp(Treatment = "Tukey")))

# Create a sensible 
mean_ci <- function(x){
  mean_se(x, mult = 1.96)
}

# Save the multiple comparison output for plotting
TukeyLabels <- cld(glht(mod1, linfct = mcp(Treatment = "Tukey")), level = 0.05)


TukeyLabels <- data.frame(Treatment = names(TukeyLabels$mcletters$Letters)
           , Letters = TukeyLabels$mcletters$Letters
           , CAZymes = rbindlist(tapply(tf$CAZymes, tf$Treatment, mean_ci))$ymax
           )

#tf <- merge(tf, TukeyLabels, by = "Treatment")
# Find the mean and confidence intervals for each group



BlackoutPalette <- c('#c2a5cf'
                     ,'#7b3294'
                     ,'#a6dba0'
                     ,'#008837')

p <- ggplot(tf, aes(x = Treatment, y = CAZymes
               , colour = Treatment)) + 
  geom_point(colour = 'grey50', position = position_jitter(0.15)) +
  stat_summary(fun.data = mean_ci, size = 1.1) +
  theme_classic() + 
  scale_colour_manual(breaks = levels(tf$Treatment)
                      , values = BlackoutPalette) +
  scale_x_discrete(breaks=levels(tf$Treatment)
                   , labels=c("1-Year\nFallow", "10-Year\nFallow"
                              , "1-Year\nGrassland", "10-Year\nGrassland")) +
  guides(colour = 'none') +
  ylab("Normalized reads mapping to CAZymes") +
  geom_text(data = TukeyLabels, aes(label = Letters
                                    , x = Treatment
                                    , y = CAZymes + 30), colour = 'black')
  
p  

# Save the figure
ggsave("Figures/CAZy/CAZyme_Abundance_Treatment_Metagenome.pdf", plot = p
         , device = "pdf", width = 4, height = 4, units = "in")


######

# Which are the CAZymes which we are interested in
cr <- unique(c(as.character(SLR$CAZyme), as.character(SMR$CAZyme)))
cr

ResponsiveCazymes <- df[, names(df) %in% cr]
UnresponsiveCazymes <- df[, names(df) %ni% cr]
Meta <- UnresponsiveCazymes[1:3]
UnresponsiveCazymes <- UnresponsiveCazymes[4:length(UnresponsiveCazymes)]

#########################

# Check Distribution of reads 
# between responsive and unresponsive CAZymes

par(mfrow=c(2,1))
hist(unlist(UnresponsiveCazymes), breaks = 120, xlim = c(0, 40)
     )
hist(unlist(ResponsiveCazymes)#, breaks = 12
     , xlim = c(0, 40)
     )
par(mfrow=c(1,1))

# Many low abundance responsive cazymes

# Compare the two groups:

UnresponsiveCazymes$Treatment <- Meta$Treatment
ResponsiveCazymes$Treatment <- Meta$Treatment


# Add a grouping factor to them all and make them plotable
uc <- melt(UnresponsiveCazymes, id = "Treatment")
head(uc)
uc$Pop <- "Other CAZymes"


rc <- melt(ResponsiveCazymes, id = "Treatment")
head(rc)
rc$Pop <- rc$variable


rf <- rbind(rc,uc) # Combine them
rf
head(rf)

# Plot the data, adding a line for the 10% quantile as
# a cutoff for read abundance

BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

ggplot(rf, aes(x = Pop, y = log10(value+1))) + 
  geom_point(position = 'jitter'
             , alpha = 0.3
             , aes(colour = Treatment)
             , shape = 16) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = log10(quantile(uc$value[uc$value!=0], 0.1)+1)
             , linetype = "dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90
                                   , hjust = 1)) + 
  scale_colour_manual(breaks = rf$Treatment, values = BlackoutPalette) + 
  ylab(expression("Normalized Read Abundance")) +
  xlab("CAZyme") + 
  ylim(0, 0.5) 


head(rf)

rf_s <- summaryBy(value ~ Pop, data = rf, FUN = median)

rf_sb <- rf_s[rf_s$value.median > quantile(uc$value[uc$value!=0], 0.1),]

# Keep a list of CAZymes where the median value > 10% of all other CAZymes
write.csv(rf_sb$Pop[-21]
          , "CAZy/Cleaned_Data/ResponsiveCAZymes_FurtherExploration.csv"
          , row.names = FALSE)

head(df[1:5])

#### Describe the composition of CAZymes ####

# Calculate The treatment level means for each CAZyme
cf <- melt(summaryBy(. ~ Treatment, data = df[-c(2,3)], FUN = mean, keep.names = TRUE)
           , id = "Treatment")

# Calculate total number of CAZyme reads per sample (for calculating percentages)
SampleCAZymeReads <- summaryBy(value ~ Sample, FUN = sum
                               , data = melt(summaryBy(. ~ Sample, data = df[c(2,4:length(df))]
                                                       , FUN = sum, keep.names = TRUE)
                                             , id = "Sample")
                               )

# How many normalized reads mapped to CAZymes?
head(df[1:5])

CAZymeReads <- sum(unlist(df[4:length(df)]))
CAZymeReads

# Calculate the percentage of reads mapping to each type of CAZyme (of reads mapping to CAZymes)
dfc <- df
dfc[4:length(dfc)] <- dfc[4:length(dfc)] / SampleCAZymeReads$value.sum * 100
head(dfc)

head(dfc[1:5])

source("Functions/ci.R")

##### What were the abundant CAZy families? ####
# Mean percentage of CAZymes overall
CAZyPercentage <- apply(dfc[4:length(dfc)], 2, mean)
CAZyPercentageCI <- apply(dfc[4:length(dfc)], 2, sd)

AbundantNames <- names(CAZyPercentageCI[names(CAZyPercentageCI) %in% 
                                          names(CAZyPercentage[CAZyPercentage > 4])
                                        ])
# Mean
round(rev(sort(CAZyPercentage[CAZyPercentage > 4])), 2)
# CI
round(rev(sort(CAZyPercentageCI[names(CAZyPercentageCI) %in% AbundantNames])), 2)

# What percentage of CAZymes belongs to other CAZy families?
sum(CAZyPercentage[CAZyPercentage < 4])
length(CAZyPercentage[CAZyPercentage < 4])

#### What percentage of reads mapped to GH, AA, or CBMs? ####

# GH
sum(CAZyPercentage[grep("GH", names(CAZyPercentage))])
length(grep("GH", names(CAZyPercentage)))

# AA
sum(CAZyPercentage[grep("AA", names(CAZyPercentage))])
length(grep("AA", names(CAZyPercentage)))

# CBM
sum(CAZyPercentage[grep("CBM", names(CAZyPercentage))])
length(grep("CBM", names(CAZyPercentage)))

sum(CAZyPercentage[grep("GH", names(CAZyPercentage))]) + sum(CAZyPercentage[grep("AA", names(CAZyPercentage))]) +sum(CAZyPercentage[grep("CBM", names(CAZyPercentage))])

##### Percentage of CAZymes per treatment #####

#### Data Sorting #####
# Find the means and 95% CIs
cf <- summaryBy(. ~ Treatment, data = dfc[c(1,4:length(dfc))]
          , FUN = c(mean, ci), keep.names = TRUE)

# Put them in plottable format
cf_m <- melt(cf[, grep("Treatment|mean", names(cf))], id = "Treatment")
cf_i <- melt(cf[, grep("Treatment|ci", names(cf))], id = "Treatment")
# Make names pretty
names(cf_m)[3] <- "Abundance"
names(cf_i)[3] <- "Abundance.ci"
head(cf_m)
head(cf_i)

length(unique(cf_m$variable))

# sort out the CAZyme names
cf_m$variable <- sub("\\.mean", "", as.character(cf_m$variable))
cf_i$variable <- sub("\\.ci", "", as.character(cf_i$variable))

# Whack 'em together and name stuff sensibly
cf_l <- cbind(cf_m, cf_i$Abundance.ci)
names(cf_l)[4] <- "Abundance.ci"
head(cf_l)

# Find the cazymes where the mean abundance across all samples > 1% of CAZymes
AbundantCAZymes <- unique(cf_l$variable)[tapply(cf_l$Abundance
                                                , INDEX = cf_l$variable
                                                , FUN = mean) > 1]

# keep the abundant CAZymes
cf_a <- cf_l[cf_l$variable %in% AbundantCAZymes, ]

table(cf_a$variable)

nrow(cf_m)
nrow(cf_l)

tail(cf_l[order(cf_l$Abundance),])


#### Plotting abundance of abundant CAZymes between treatments ####
ggplot(cf_a
       , aes(x = Treatment, y = Abundance, fill = Treatment)) +
  geom_bar(stat = 'identity') +
  guides(fill = 'none') +
  facet_wrap(. ~ variable) +
  geom_errorbar(aes(ymax = Abundance + Abundance.ci
                    , ymin = Abundance - Abundance.ci)
                , width = 0.3)

#### Test for differences in the abundance of abundant CAZymes ####
head(df[1:5])

AbundantCAZymeResponses <- apply(df[4:length(df)]#[names(df) %in% AbundantCAZymes]
                                 , 2
      , FUN = function(x){
  m <- lm(x ~ df$Treatment)
  Results <- drop1(m, test = 'Chisq')
  Results[5][2,1]
})

AbundantCAZymeResponses <- AbundantCAZymeResponses[names(AbundantCAZymeResponses) %in% AbundantCAZymes]

####
par(mfrow=c(2,2))
plot(lm(df[,names(df) %in% AbundantCAZymes[1]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[2]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[3]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[4]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[5]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[6]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[7]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[8]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[9]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[10]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[11]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[12]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[13]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[14]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[15]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[16]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[17]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[18]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[19]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[20]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[21]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[22]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[23]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[24]] ~ df$Treatment))
plot(lm(df[,names(df) %in% AbundantCAZymes[25]] ~ df$Treatment))
par(mfrow=c(1,1))

####

hist(p.adjust(AbundantCAZymeResponses, method = 'bonferroni'), breaks = 50)
P_Vals <- p.adjust(AbundantCAZymeResponses, method = 'bonferroni')
head(cf_l)
P_Vals <- data.frame(variable = names(P_Vals), p_adj = P_Vals, Treatment = "b-New")

# Correct order for CAZymes
OrderIWant <- rev(names(sort(tapply(cf_a$Abundance, cf_a$variable, FUN = mean))))
cf_a$variable <- factor(cf_a$variable, levels = OrderIWant)

RawData <- melt(dfc[names(dfc) %in% c("Treatment", unique(as.character(cf_a$variable)))])

# Whack the p-values into the data frame
cf_a2 <- merge(cf_a, P_Vals, by = c("variable", "Treatment"), all.x = TRUE)

TotalCazyme <- sum(RawData$value)

ggplot(cf_a2
       , aes(x = Treatment, y = Abundance, fill = Treatment)) +
  geom_bar(stat = 'identity') +
  guides(fill = 'none') +
  facet_wrap(. ~ variable, scales = 'free_y') +
  geom_errorbar(aes(ymax = Abundance + Abundance.ci
                    , ymin = Abundance - Abundance.ci)
                , width = 0.3) +
  geom_text(aes(x = "b-Old", y = 1.5*(Abundance + Abundance.ci)
                , label = ifelse(is.na(p_adj)
                                 , yes = ""
                                 , no = paste0("p = "
                                               , round(p_adj, 3)
                                               )
                                 )
                )) +
  geom_point(data = RawData
             , aes(x = Treatment, y = value#/TotalCazyme * 100
                   ))



  #ylim(0,9) +
  #scale_y_continuous(trans = 'log10')







 