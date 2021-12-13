# Take the means of each CAZyme in each treatment, and compare to the mean of 
# all other treatments combined. This shows whether the group was divergent
# or not. Then run Gaussian GLMs for this and correct the p-values. Plot 
# corrected p-values against log2-fold change (volcano plot), then identify 
# CAZymes which were responsive to the treatment. Then plot absolute 
# differences in richness for each of these CAZymes and find a suitable model 

rm(list = ls())

library(roperators)

#### Data import  and checking ####
df <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZyFamilyRichness.csv")

df[1:6]

#### Data Cleaning and checking ####

df[, 4:length(df)] <- df[, 4:length(df)] / (df$Reads/1e9) # Normalise richness per Gbp of library for each sequence

head(df)[1:6]

#### Take the means of each CAZyme in each treatment, #### 
####     and mean of all other treatments combined. ####

BOld <- ifelse(df$Treatment == "b-Old", yes = "b-Old", no = "Pop")
BNew <- ifelse(df$Treatment == "b-New", yes = "b-New", no = "Pop")
COld <- ifelse(df$Treatment == "c-Old", yes = "c-Old", no = "Pop")
CNew <- ifelse(df$Treatment == "c-New", yes = "c-New", no = "Pop")
BOld <- relevel(factor(BOld), ref = 'Pop')
BNew <- relevel(factor(BNew), ref = 'Pop')
COld <- relevel(factor(COld), ref = 'Pop')
CNew <- relevel(factor(CNew), ref = 'Pop')

# Calculate fold change data for each treatment + for each gene

ci <- function(x){
  ci <- 1.96*(sd(x)/sqrt(length(x)))
  return(ci)
}

# Blackout Old
BO <- apply(df[,4:length(df)], MARGIN = 2
            , FUN = function(x){
              BO_mean = mean(x[BOld != 'Pop'])
              Pop_mean = mean(x[BOld == 'Pop'])
              BO_ci = ci(x[BOld != 'Pop'])
              Pop_ci = ci(x[BOld == 'Pop'])
              
              j <- data.frame(Treatment = "b-Old"
                              , Group = c(rep(c("b-Old", "Pop")
                                              , each = length(BO_mean)))
                              , M = c(BO_mean, Pop_mean)
                              , CI = c(BO_ci, Pop_ci)
              )
              
              return(j)
            })

# Blackout New
BN <- apply(df[,4:length(df)], MARGIN = 2
            , FUN = function(x){
              BO_mean = mean(x[BNew != 'Pop'])
              Pop_mean = mean(x[BNew == 'Pop'])
              BO_ci = ci(x[BNew != 'Pop'])
              Pop_ci = ci(x[BNew == 'Pop'])
              
              j <- data.frame(Treatment = "b-New"
                              , Group = c(rep(c("b-New", "Pop")
                                              , each = length(BO_mean)))
                              , M = c(BO_mean, Pop_mean)
                              , CI = c(BO_ci, Pop_ci)
              )
              
              return(j)
            })

# Blackout Old
CO <- apply(df[,4:length(df)], MARGIN = 2
            , FUN = function(x){
              BO_mean = mean(x[COld != 'Pop'])
              Pop_mean = mean(x[COld == 'Pop'])
              BO_ci = ci(x[COld != 'Pop'])
              Pop_ci = ci(x[COld == 'Pop'])
              # RichnessCI = c(ci(x[BOld != 'Pop'])
              #                  , ci(x[BOld == 'Pop']))
              
              
              j <- data.frame(Treatment = "c-Old"
                              , Group = c(rep(c("c-Old", "Pop")
                                              , each = length(BO_mean)))
                              , M = c(BO_mean, Pop_mean)
                              , CI = c(BO_ci, Pop_ci)
              )
              
              return(j)
            })

# Control New
CN <- apply(df[,4:length(df)], MARGIN = 2
            , FUN = function(x){
              BO_mean = mean(x[CNew != 'Pop'])
              Pop_mean = mean(x[CNew == 'Pop'])
              BO_ci = ci(x[CNew != 'Pop'])
              Pop_ci = ci(x[CNew == 'Pop'])
              # RichnessCI = c(ci(x[BOld != 'Pop'])
              #                  , ci(x[BOld == 'Pop']))
              
              
              j <- data.frame(Treatment = "c-New"
                              , Group = c(rep(c("c-New", "Pop")
                                              , each = length(BO_mean)))
                              , M = c(BO_mean, Pop_mean)
                              , CI = c(BO_ci, Pop_ci)
              )
              
              return(j)
            })

# Now make these into a single dataframe
BO <- do.call("rbind", BO)
BN <- do.call("rbind", BN)
CO <- do.call("rbind", CO)
CN <- do.call("rbind", CN)

# Add a CAZyme column--easier to do this now than later
BO$CAZyme <- sub("?\\.[1-2]$", "", row.names(BO))
BN$CAZyme <- sub("?\\.[1-2]$", "", row.names(BN))
CO$CAZyme <- sub("?\\.[1-2]$", "", row.names(CO))
CN$CAZyme <- sub("?\\.[1-2]$", "", row.names(CN))


rawMeans <- do.call("rbind", list(BO,BN,CO,CN))

# Check all CAZyme names are accurate
any(rawMeans$CAZyme %ni% names(df[4:length(df)]))
head(rawMeans)

# Write this to a file for later analysis. Useful for plotting 
# absolute differences between treatment and the rest of the population
write.csv(rawMeans, "CAZy/Cleaned_Data/Henfaes_MeanCAZymeRichness_Treatment_v_Pop.csv"
          , row.names = FALSE)

# Now calculate log2 fold change +1
# Subset to the treatment vs the control. Thankfully these stay in order
# So dividing the first ro of one by the other is from the same data output
Trt <- rawMeans[rawMeans$Group != 'Pop',]
Pop <- rawMeans[rawMeans$Group == 'Pop',]

# Check that this has worked
unique(Trt$Group) # should return "b-Old b-New c-Old c-New"
unique(Pop$Group) # Should return "Pop"

any(Trt$CAZyme %ni% names(df[4:length(df)])) # All CAZymes accounted for-- should ensure merging is possible

# Calculate log2 fold change +1
Trt$log2FC_raw <- log2(Trt$M+1) - log2(Pop$M+1)

# Look at the data
head(Trt)
Trt[Trt$CAZyme == "AA10",] # Check all factor levels of Treatment are present

# Check how responsive the CAZymes were
hist(Trt$log2FC_raw)

#### Run Gaussian GLMs on all CAZymes and correct the p-values. ####

# Fit all of the models and take the drop1 results, and the model 
# coefficients

Model_P_Values <- apply(X = df[,c(4:length(df))]
                        , MARGIN = 2
                        , FUN = function(x){
                          # m1 <- glm(df$AA10 ~ BOld)
                          
                          # Fit models for each of the treatments
                          m1 <- glm(x ~ BOld)
                          m2 <- glm(x ~ BNew)
                          m3 <- glm(x ~ COld)
                          m4 <- glm(x ~ CNew)
                          
                          # Calculate p-values
                          m1_p <- drop1(m1, test = "LR")#[5][,1][2]
                          m2_p <- drop1(m2, test = "LR")#[5][,1][2]
                          m3_p <- drop1(m3, test = "LR")#[5][,1][2]
                          m4_p <- drop1(m4, test = "LR")#[5][,1][2])
                          
                          # Results data frame
                          resf <- data.frame(
                            # Which subpopulation is being tested
                            
                            Treatment = c("b-Old", "b-New", "c-Old", "c-New")
                            
                            # p-values
                            , p = c(m1_p$`Pr(>Chi)`[2]
                              , m2_p$`Pr(>Chi)`[2]
                              , m3_p$`Pr(>Chi)`[2]
                              , m4_p$`Pr(>Chi)`[2]
                            )
                            # Estimate of the response (already on response scale)
                            , Est = 
                              c(m1$coefficients[2]
                                , m2$coefficients[2]
                                , m3$coefficients[2]
                                , m4$coefficients[2])
                            
                            # Estimate of the Intercept (already on response scale) 
                            , Int = 
                              c(m1$coefficients[1]
                                , m2$coefficients[1]
                                , m3$coefficients[1]
                                , m4$coefficients[1])
                            
                            , EstSE = 
                              c(summary(m1)$coefficients[,2][2]
                                , summary(m2)$coefficients[,2][2]
                                , summary(m3)$coefficients[,2][2]
                                , summary(m4)$coefficients[,2][2])
                            
                            , IntSE = 
                              c(summary(m1)$coefficients[,2][1]
                                , summary(m2)$coefficients[,2][1]
                                , summary(m3)$coefficients[,2][1]
                                , summary(m4)$coefficients[,2][1])
                            
                          )   
                          
                          # adjusted p-value for multiple testing - should be done afterwards
                          # resf$p_adj <- p.adjust(p = resf$p, method = "holm")
                          
                          resf
                        })

Model_P_Values[1:3] # Check output

# Put data from all CAZymes into a single data frame
rft <- do.call("rbind", Model_P_Values)
head(rft)

# Add a CAZyme column
CAZyme <- rep(names(Model_P_Values), sapply(Model_P_Values, nrow))
any(CAZyme %ni% names(df[4:length(df)])) # Check all CAZyme names are valid
rft <- cbind(CAZyme, rft) # Append the column
rft$Treatment
nrow(rft)/4

# Check the output
head(rft)

# Which correction method and alpha to choose?

pholm <- p.adjust(rft$p, method = 'holm')
pBH <- p.adjust(rft$p, method = 'BH')



# Plot all possible values
plot(sort(pholm) ~ sort(rft$p)
     , bty = 'l'
     , xaxs = 'i'
     , yaxs = 'i'
     , col = 'red'
     , lwd = 3
     , type = 'n'
     , xlab = "p value"
     , ylab = "Adjusted p value")
# Add significance levels
polygon(x = c(0,0,1,1), y = c(-1,0.1, 0.1,-1), border = NA
        , col = scales::alpha("black",0.3))
polygon(x = c(0,0,1,1), y = c(-1,0.05, 0.05,-1), border = NA
        , col = scales::alpha("black",0.3))
polygon(x = c(0,0,1,1), y = c(-1,0.01, 0.01,-1), border = NA
        , col = scales::alpha("black",0.3))
# Add p-values and adjusted p-values
abline(a = 0, b = 1, lwd = 3, col = "black") # p v p
lines(sort(pholm) ~ sort(rft$p), lwd = 3, col = 'red3')
lines(sort(pBH) ~ sort(rft$p), lwd = 3, col = 'darkgreen')

# Add a legend
legend("bottomright", legend = c("Unadjusted"
                                 , "False Positive Rate"
                                 , "Holm")
       , lty = 1
       , lwd = 3
       , col = c("black", "darkgreen", "red3"))

# How many values for further exploration do we get from each cutoff?
# Unadjuasted
length(rft$p[rft$p < 0.01]);length(rft$p[rft$p < 0.05]);length(rft$p[rft$p < 0.1])

# holm adjuasted
length(pholm[pholm < 0.01]);length(pholm[pholm < 0.05]);length(pholm[pholm < 0.1])

# FPR adjuasted
length(pBH[pBH < 0.01]);length(pBH[pBH < 0.05]);length(pBH[pBH < 0.1])

# This is a screening excercise so it's not super important. 
# The aim is to identify which CAZymes have responded significantly
# to the treatment, and so a non-conservative approach here is 
# fine.
# Could even set alpha to 0.1 as replication is so low and the 
# system is so variable.

# Go with FPR adjustment, and alpha of 0.1. Variable system, 
# exploratory analysis

rft$p_adj <- p.adjust(rft$p, method = 'BH')

## Calculate fold change. Look up how to log REAAALY big negative numbers 
hist((rft$Int))
min((rft$Int)) # 0 bounded
hist(rft$Est) # Huge tails in bioth directions. Needs transforming so that 
              # log2 can be applied.

hist(rft$Int + rft$Est, breaks = 100)
min(rft$Int + rft$Est)  # hmmmmm. Artefact of using a gaussian model...?

# +1, calculate fold change
rft$log2FC_model <- log2(rft$Est+rft$Int+1) - log2(rft$Int+1)
hist(rft$log2FC_model)


#### Plot corrected p-values against log2-fold change (volcano plot). #### 

# First merge the data frames by CAZyme and treatment
head(Trt)
head(rft)

# Check again that all treatmernt levels are available
Trt[Trt$CAZyme == 'AA10',]

# Perform the merge
rft <- merge(Trt[,c(1,5,6)], rft
              , by = c("CAZyme", "Treatment"))
# Check the output
head(rft)
write.csv(rft
          , "CAZy/Cleaned_Data/Henfaes_CAZymeResponse_p_FPR.csv"
          , row.names = FALSE)

# Check agreement between model and actual log2 fold change
par(mfrow = c(2,1))
hist(rft$log2FC_raw, xlim = c(-6,2))
hist(rft$log2FC_model, xlim = c(-6,2))
par(mfrow = c(1,1))

# First volcano plot--model predictions and unadjusted p values
plot(-log10(p) ~ log2FC_model, data = rft)
abline(v = -1, lty = 2)
abline(v = 1, lty = 2)
abline(h= -log10(0.05), lty = 2)
abline(h= -log10(0.01), lty = 2)

# Second volcano plot--data-based fold change and unadjusted p values
plot(-log10(p) ~ log2FC_raw, data = rft
     , pch = 16
     , type = 'n'
     , xlab = expression("Log"[2]*" Fold Change in Richness")
     , ylab = expression("-Log"[10]*"(p value)"))

points(-log10(p) ~ log2FC_raw
       , data = rft[rft$p_adj < 0.05
                    & rft$log2FC_raw < -1,]
       , pch = 16, col = scales::alpha(1:4, 0.6)[rft$Treatment]
       )
points(-log10(p) ~ log2FC_raw
       , data = rft[rft$p_adj >= 0.05
                    #,
                    ,]
       , pch = 16, col = scales::alpha(1:4, 0.2)[rft$Treatment]
)

abline(v = -1, lty = 2)
abline(v = 1, lty = 2)
abline(h = -log10(0.01), lty = 2)
abline(h = -log10(0.05), lty = 2)

# Good agreement from model and raw results (as it should be)
# Try this with results from Gamma(log) model.

# Third volcano plot. The one we will use--data-based fold change and FPR adjusted p values

# Create a palette 
BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

# Now make an empty plot with appropriate limits
plot(log10(-log10(p_adj)) ~ log2FC_raw
     , data = rft
     , pch = 16
     , type = 'n'
     , xlab = expression("Log"[2]*" Fold Change in Richness")
     , ylab = expression("Log"[10]*"( -Log"[10]*"(p"[adj]*" value) )"))

# Add the significantly less rich CAZymes
points(log10(-log10(p_adj)) ~ log2FC_raw
       , data = rft[rft$p_adj < 0.05
                    & rft$log2FC_raw < -1,]
       , pch = 16
       , col = scales::alpha(BlackoutPalette, 0.7)[rft$Treatment]
)

# Add the significantly more rich CAZymes
points(log10(-log10(p_adj)) ~ log2FC_raw
       , data = rft[rft$p_adj < 0.05
                    & rft$log2FC_raw > 1,]
       , pch = 16
       , col = scales::alpha(BlackoutPalette, 0.7)[rft$Treatment]
)

# Add CAZymes which were significantly less than 1 fold 
# different from the rest of the population
points(log10(-log10(p_adj)) ~ log2FC_raw
       , data = rft[rft$p_adj < 0.05
                    & rft$log2FC_raw < 1
                    & rft$log2FC_raw > -1,]
       , pch = 16
       , col = scales::alpha(BlackoutPalette, 0.25)[rft$Treatment]
)

# Add points which were nonsignificant
points(log10(-log10(p_adj)) ~ log2FC_raw
       , data = rft[rft$p_adj >= 0.05,]
       , pch = 16
       , col = scales::alpha(BlackoutPalette, 0.25)[rft$Treatment]
)

# Add cutoff lines
abline(v = -1, lty = 2)
abline(v = 1, lty = 2)
#abline(h = -log10(0.01), lty = 2)
#abline(h = log10(-log10(0.1)), lty = 2)
abline(h = log10(-log10(0.05)), lty = 2)
#abline(h = log10(-log10(0.01)), lty = 2)

# Add a legend
legend("bottomright"
       , legend = levels(rft$Treatment)
       , pch = 16
       #,  = 3
       , col = BlackoutPalette)


#### Identify CAZymes which were responsive to the treatment (p_adj<0.01). ####

head(rft)

# Significantly less rich CAZymes
SLR <- rft[rft$p_adj < 0.1 & rft$log2FC_raw < -1,]

# Significantly more rich CAZymes
SMR <- rft[rft$p_adj < 0.1 & rft$log2FC_raw > 1,]

length(SLR[,1])
length(SMR[,1])

table(SLR$Treatment)
table(SMR$Treatment)

SLR$CAZyme
SMR$CAZyme

length(grep("CBM", c(SMR$CAZyme, SLR$CAZyme)))/length(c(SMR$CAZyme, SLR$CAZyme))*100

# ResponsiveGroups <-

#### Plot absolute differences in richness for each responsive CAZymes ####
head(df)[1:6]
ggplot(df, aes(x = Treatment, y = AA12)) +
  geom_boxplot() + 
  stat_summary(fun.data = 'mean_se', colour = 'red')



####

plotframe <- df[, names(df) %in% c("Treatment"#, "GH135"
                                   #, "GH35"
                                   , "CBM4.GH9" # Endoglucanase
                                   , "CBM22.CBM6.GH10" # Xylanase
                                   )]

plotframe <- reshape2::melt(plotframe, id = "Treatment")
names(plotframe) <- c("Treatment", "CAZyme", "Richness")

BlackoutPalette <- c('#c2a5cf'
                     ,'#7b3294'
                     ,'#a6dba0'
                     ,'#008837')

ggplot(plotframe, aes(x = Treatment, y = Richness)) +
  geom_point(alpha = 0.5, colour = 'grey50') + 
  stat_summary(fun.data = 'mean_se'
               , aes(colour = Treatment)
               , size = 1.2
               ) +
  facet_wrap( ~ CAZyme, scales = "free", nrow = 2) + 
  scale_colour_manual(values = BlackoutPalette) + 
  theme_bw() +
  theme(panel.grid = element_blank()
        , text = element_text(size = 18))






#####################
ggplot(df, aes(x = Treatment, y = GH135)) +
  geom_point() + 
  stat_summary(fun.data = 'mean_se', colour = 'red')

ggplot(df, aes(x = Treatment, y = GH35)) +
  geom_point() + 
  stat_summary(fun.data = 'mean_se', colour = 'red')

ggplot(df, aes(x = Treatment, y = CBM4.GH9)) +
  geom_point() + 
  stat_summary(fun.data = 'mean_se', colour = 'red')

ggplot(df, aes(x = Treatment, y = CBM22.CBM6.GH10)) +
  geom_point() + 
  stat_summary(fun.data = 'mean_se', colour = 'red')


#### Find a suitable model for each of these (and correct p values) #### 

