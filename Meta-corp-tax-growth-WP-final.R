rm(list = ls()) #clear list

#automatic installation of required packages
packages <- c("xlsx","calibrate","stargazer","sandwich","lmtest","getopt","CausalGAM","ggplot2","reshape2","xts",
              "lattice","gridExtra","gtable","plm","lfe","lmtest","car","tis","foreign","MASS","quantreg","ggrepel",
              "dplyr","stringr","datasets","rio","psych","systemfit","MatchIt","CRTgeeDR","eurostat","plyr","zoo","ggthemes",
              "robumeta","metafor","dplyr","clubSandwich","Hmisc","metafor","pracma","pkgs","broom","sjPlot", "here", "data.table", "pscore")
ipak(packages)

#load packages
library(xlsx) #Excel-Paket laden
library(calibrate) #Laden des Pakets, das f??r Datenbeschriftung n??tig ist
library (stargazer) #Laden des Pakets, mit dem R-Regressionsoutput in Latex-Tabellen ??bergef??hrt werden kann
library(sandwich)
library(lmtest)
library(getopt)
library(CausalGAM)
library(ggplot2)
library(reshape2)
library(xts)
library(lattice)
library(gridExtra)
library(gtable)
library(plm)
library(lfe)
library(lmtest)
library(car)
library(tis)
library(foreign)
library(MASS)
library(quantreg)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggplot2)
library(datasets)
library(rio)
library(psych)
library(systemfit)
library(foreign)
library(MatchIt)
library(CRTgeeDR)
library(eurostat)
library(plyr)
library(zoo)
library(ggthemes)
library("robumeta")
library("metafor")
library("dplyr")
library(clubSandwich)
library(Hmisc)
library(metafor)
library(pracma)
library(broom)
library(sjPlot)
library(here)
library(data.table)
library(pscore)

#Load data
dat <- fread(here("data/Coding-taxes-and-growth-WP-final.csv"))

#additional variables and data transformations
#calculate the partial correlation coefficient
dat$PartialCorrelationCoefficient<- dat$TstatisticCorrected / (sqrt((dat$TstatisticCorrected^2)+dat$DegreesofFreedom))

#calculate standard error of the partial correlation coefficient
dat$StandardErrorPartialCorrelation <- sqrt((1-(dat$PartialCorrelationCoefficient)^2)/dat$DegreesofFreedom)

#Precision of partial correlation
dat$PrecSE <- 1 / dat$StandardErrorPartialCorrelation
#Inverse standard error of partial correlation
dat$InverseSE <- 1 / dat$StandardErrorCorrected
#Inverse standard error of corrected correlation coefficient
dat$InverseSECorrected <- 1 / dat$StandardErrorCorrected

#Variance partial correlation
dat$Variance <- dat$StandardErrorPartialCorrelation^2
#Variance of corrected correlation coefficient
dat$VarianceSECorrected <- dat$StandardErrorCorrected^2
#PrecVariance partial correlation
dat$PrecVariance <- 1 / dat$Variance
#PrecVariance corrected
dat$PrecVarianceSECorrected <- 1 / dat$VarianceSECorrected

#Proxy for standard error
dat$proxySE <- 1 / (sqrt(dat$Observations))

#taking logs
#dat$YearofPublication <- log(dat$YearofPublication)
dat$NumberofCountries <- log(dat$NumberofCountries)
dat$Citations <- log(dat$Citations)

dat$MeanYearData<- (dat$StartYearCorr+dat$EndYearCorr)/2
dat$InstrumentSE <- 1 / (sqrt(dat$DegreesofFreedom))
dat$InstrumentVariance <- dat$InstrumentSE^2
dat$PrecInstrumentVariance <- 1 / dat$InstrumentVariance
dat$PrecNumberRegressions <- 1 / dat$NumberRegressions
#normalize impact factor
dat$JournalImpactFactor <- as.numeric(dat$JournalImpactFactor)
dat$MaxImpactFactor <- max(dat$JournalImpactFactor)
dat$MaxImpactFactor <- as.numeric(dat$MaxImpactFactor)
dat$NormalizedImpactFactor <- dat$JournalImpactFactor / max(dat$JournalImpactFactor)

#alternative weight based on degrees of freedom
dat$weightDF <- 1/dat$DegreesofFreedom
dat$alternativeweights <- 1/dat$NumberRegressions

#transform data
dat_long <- melt(dat, id=1:69)

#sub-sets of the data
#corrected coefficients
dat_long_corrected <- subset(dat_long, YNCorrected %in% c('1'))
dat_long_corrected$PrecStandardErrorCorrected <- 1/dat_long_corrected$StandardErrorCorrected

#preferred estimates only
dat_long_preferred <- subset(dat_long_corrected, Preferred %in% c('1'))

#without inferior estimates
dat_long_without_inferior <- subset(dat_long_corrected, Preferred %in% c('1', '0'))

#Descriptive statistics
#all-set

#minimum, maximum, standard deviation
max(dat_long_corrected$CorrelationCoefficientCorrected)
min(dat_long_corrected$CorrelationCoefficientCorrected)
sd(dat_long_corrected$CorrelationCoefficientCorrected)
mean(dat_long_corrected$CorrelationCoefficientCorrected)
median(dat_long_corrected$CorrelationCoefficientCorrected)

#winsorising at the 2nd and 98th percentiles
dat_2p98p <- dat_long

dat_2p98p$StandardErrorCorrected <- winsorizor(dat_long$StandardErrorCorrected, c(0.02), na.rm=TRUE)
dat_2p98p$CorrelationCoefficientCorrected <- winsorizor(dat_long$CorrelationCoefficientCorrected, c(0.02), na.rm=TRUE)
dat_2p98p$InstrumentSE <- winsorizor(dat_long$InstrumentSE, c(0.02), na.rm=TRUE)
dat_2p98p$InstrumentVariance <- dat_2p98p$InstrumentSE^2
dat_2p98p$PrecInstrumentVariance <- 1 / dat_2p98p$InstrumentVariance

#calculate standard errors corrected after precision was winsorised
dat_2p98p$InverseSECorrected <- 1 / dat_2p98p$StandardErrorCorrected

#Variance of corrected correlation coefficient
dat_2p98p$VarianceSECorrected <- dat_2p98p$StandardErrorCorrected^2
#PrecVariance corrected
dat_2p98p$PrecVarianceSECorrected <- 1 / dat_2p98p$VarianceSECorrected

#Funnel plot
plot_funnel_CorrelationCoefficientCorrected_2p98p <- ggplot(data=dat_2p98p,
                                                      aes(x=CorrelationCoefficientCorrected, y=InverseSECorrected)) +
  geom_point(size=1.5, color="blue") +
  xlab("Standardised coefficient") +
  ylab("Inverse of standard error (precision)") +
  ggtitle("")+
  theme(title=element_text(size=10, face='bold'))+
  theme(panel.background = element_rect(fill = "white")) +
  theme(legend.position="bottom")+
  theme(legend.title=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.1),
        axis.line.y = element_line(color="black", size = 0.1)) +
  geom_vline(xintercept=0, colour="black", linetype=2)+
  geom_vline(xintercept=-0.01869439, colour="black", linetype=1)+
  theme(legend.text = element_text(colour="black", size = 4))+
  theme(axis.text.x=element_text(size=13))+
  theme(axis.title.x=element_text(size=13)) +
  theme(axis.text.y=element_text(size=13))+
  theme(axis.title.y=element_text(size=13))
plot_funnel_CorrelationCoefficientCorrected_2p98p

#save figure
ggsave(plot = plot_funnel_CorrelationCoefficientCorrected_2p98p, 
       filename = paste0(here("figures/funnel_w2_98"), 
                         ".pdf"),
       width = 6.5, height = 4)

#winsorise at the 5th and 95th percentiles
dat_5p95p <- dat_long

#winsorise StandardErrorCorrected (p2, p98)
quantile(dat$StandardErrorCorrected, c(0.05, 0.95), na.rm=TRUE)

dat_5p95p$StandardErrorCorrected <- ifelse(dat_5p95p$StandardErrorCorrected>0.182261307, 0.182261307, dat_5p95p$StandardErrorCorrected)
dat_5p95p$StandardErrorCorrected <- ifelse(dat_5p95p$StandardErrorCorrected<0.002584475, 0.002584475, dat_5p95p$StandardErrorCorrected)

#calculate standard errors corrected after precision was winsorised
dat_5p95p$InverseSECorrected <- 1 / dat_5p95p$StandardErrorCorrected

#Variance of corrected correlation coefficient
dat_5p95p$VarianceSECorrected <- dat_5p95p$StandardErrorCorrected^2
#PrecVariance corrected
dat_5p95p$PrecVarianceSECorrected <- 1 / dat_5p95p$VarianceSECorrected

#winsorise CorrelationCoefficientCorrected (p2, p98)
dat_5p95p$CorrelationCoefficientCorrected <- ifelse(dat_5p95p$CorrelationCoefficientCorrected>0.08133333, 0.08133333, dat_5p95p$CorrelationCoefficientCorrected)
dat_5p95p$CorrelationCoefficientCorrected <- ifelse(dat_5p95p$CorrelationCoefficientCorrected< -0.18200000, -0.18200000, dat_5p95p$CorrelationCoefficientCorrected)

#winsorise Instrument SE
dat_5p95p$InstrumentSE <- ifelse(dat_5p95p$InstrumentSE>0.25000000, 0.25000000, dat_5p95p$InstrumentSE)
dat_5p95p$InstrumentSE <- ifelse(dat_5p95p$InstrumentSE<0.02513388, 0.02513388, dat_5p95p$InstrumentSE)

#Funnel Asymmetry Precision Effect Tests
#Table 1
#column (1)
#(precision-weighted) average
regwa <- lm(CorrelationCoefficientCorrected~1, data=dat_2p98p)
summary(regwa)
coef_test(regwa, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")
confint(regwa, level=0.95)

#column(2)
#baseline WLS
pubbias_1_CorrelationCoefficientCorrected_2p98p <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected, weights=PrecVarianceSECorrected, data=dat_2p98p)
summary(pubbias_1_CorrelationCoefficientCorrected_2p98p)

coef_test(pubbias_1_CorrelationCoefficientCorrected_2p98p, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#column (3)
#median
pubbias_1_CorrelationCoefficientCorrected_median <- lm(coeffmed ~ StandardErrorCorrected, weights=PrecVarianceSECorrected, data=dat_2p98p)
summary(pubbias_1_CorrelationCoefficientCorrected_median)

coef_test(pubbias_1_CorrelationCoefficientCorrected_median, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#column (4)
#IV
library(AER)

#with ivreg package
pubbias_1_CorrelationCoefficientCorrected_IV <- AER::ivreg(CorrelationCoefficientCorrected ~ StandardErrorCorrected | InstrumentSE, weights=1/DegreesofFreedom, data=dat_2p98p)
summary(pubbias_1_CorrelationCoefficientCorrected_IV)

coef_test(pubbias_1_CorrelationCoefficientCorrected_IV, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#column (5)
#Partial correlation
pubbias_1_partialcorr <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation, weights=PrecVariance, data=dat)
summary(pubbias_1_partialcorr)

coef_test(pubbias_1_partialcorr, vcov = "CR0", 
          cluster = dat$paperid, test = "naive-t")

#stargazer standard errors
ses_pubbias_1_CorrelationCoefficientCorrected_2p98p <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_2p98p, vcov = "CR0", 
                                                                       cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1_CorrelationCoefficientCorrected_2p98p <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_2p98p, vcov = "CR0", 
                                                                         cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_pubbias_1_CorrelationCoefficientCorrected_2p98p <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_2p98p, vcov = "CR0", 
                                                                         cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


ses_pubbias_1_CorrelationCoefficientCorrected_median <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_median, vcov = "CR0", 
                                                                    cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1_CorrelationCoefficientCorrected_median <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_median, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_pubbias_1_CorrelationCoefficientCorrected_median <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_median, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


ses_pubbias_1_CorrelationCoefficientCorrected_IV <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_IV, vcov = "CR0", 
                                                                    cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1_CorrelationCoefficientCorrected_IV <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_IV, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_pubbias_1_CorrelationCoefficientCorrected_IV <- list(coef_test(pubbias_1_CorrelationCoefficientCorrected_IV, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_pubbias_1_CorrelationCoefficientCorrected_partialcorr <- list(coef_test(pubbias_1_partialcorr, vcov = "CR0", 
                                                                   cluster = dat$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1_CorrelationCoefficientCorrected_partialcorr <- list(coef_test(pubbias_1_partialcorr, vcov = "CR0", 
                                                                     cluster = dat$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_pubbias_1_CorrelationCoefficientCorrected_partialcorr <- list(coef_test(pubbias_1_partialcorr, vcov = "CR0", 
                                                                     cluster = dat$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#stargazer table funnel asymmetry tests (table 1; column (6) regarding the non-parametric test based on Andrews and Kasy (2019 is missing here))
stargazer(pubbias_1_CorrelationCoefficientCorrected_2p98p, pubbias_1_CorrelationCoefficientCorrected_median, pubbias_1_CorrelationCoefficientCorrected_IV, pubbias_1_partialcorr, t=list(unlist(tvals_pubbias_1_CorrelationCoefficientCorrected_2p98p), unlist(tvals_pubbias_1_CorrelationCoefficientCorrected_median), unlist(tvals_pubbias_1_CorrelationCoefficientCorrected_IV), unlist(tvals_pubbias_1_CorrelationCoefficientCorrected_partialcorr)), se=list(unlist(ses_pubbias_1_CorrelationCoefficientCorrected_2p98p), unlist(ses_pubbias_1_CorrelationCoefficientCorrected_median), unlist(ses_pubbias_1_CorrelationCoefficientCorrected_IV), unlist(ses_pubbias_1_CorrelationCoefficientCorrected_partialcorr)), p=list(unlist(pvals_pubbias_1_CorrelationCoefficientCorrected_2p98p), unlist(pvals_pubbias_1_CorrelationCoefficientCorrected_median), unlist(pvals_pubbias_1_CorrelationCoefficientCorrected_IV), unlist(pvals_pubbias_1_CorrelationCoefficientCorrected_partialcorr)))

#Table 3
#column (1)
MRA_1_CorrelationCoefficientCorrected_2p98p <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending, weights=PrecVarianceSECorrected, data=dat_2p98p)
summary(MRA_1_CorrelationCoefficientCorrected_2p98p)

coef_test(MRA_1_CorrelationCoefficientCorrected_2p98p, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#column (2)
#+estimation details + data characteristics
MRA_1_CorrelationCoefficientCorrected_est_data <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending + NoOtherTaxVariables + PMG + SURE + GMM + IV + OtherEstimator + TacklingEndogeneity + CountryFixedEffects + CrossSection + USAonly + IntraNational + NoGDPpercapita, weights=PrecVarianceSECorrected, data=dat_2p98p)
summary(MRA_1_CorrelationCoefficientCorrected_est_data)

coef_test(MRA_1_CorrelationCoefficientCorrected_est_data, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#column (3)
#+publication characteristics
dat_2p98p$YearofPublicationCorr <- dat_2p98p$YearofPublication - mean(dat_2p98p$YearofPublication, na.rm=TRUE)

MRA_1_CorrelationCoefficientCorrected_pubchar <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending + YearofPublicationCorr + Citations + NormalizedImpactFactor + AuthorOECDaffiliation, weights=PrecVarianceSECorrected, data=dat_2p98p)
summary(MRA_1_CorrelationCoefficientCorrected_pubchar)

coef_test(MRA_1_CorrelationCoefficientCorrected_pubchar, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#column (4)
#median
MRA_1_CorrelationCoefficientCorrected_median<- lm(coeffmed ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending, weights=PrecVarianceSECorrected, data=dat_2p98p)
summary(MRA_1_CorrelationCoefficientCorrected_median)

coef_test(MRA_1_CorrelationCoefficientCorrected_median, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#column (5)
#partial correlation
MRA_1_partialcorr <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending, weights=PrecVariance, data=dat)
summary(MRA_1_partialcorr)

coef_test(MRA_1_partialcorr, vcov = "CR0", 
          cluster = dat$paperid, test = "naive-t")

#Robustness checks (reported in the appendix)

#Table D1

#Column (2)
#no winsorising
MRA_1_CorrelationCoefficientCorrected_now <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending, weights=PrecVarianceSECorrected, data=dat_long)
summary(MRA_1_CorrelationCoefficientCorrected_now)

coef_test(MRA_1_CorrelationCoefficientCorrected_now, vcov = "CR0", 
          cluster = dat_5p95p$paperid, test = "naive-t")

#column (3)
#winsorising 5p95p
MRA_1_CorrelationCoefficientCorrected_5p95p <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending, weights=PrecVarianceSECorrected, data=dat_5p95p)
summary(MRA_1_CorrelationCoefficientCorrected_5p95p)

coef_test(MRA_1_CorrelationCoefficientCorrected_5p95p, vcov = "CR0", 
          cluster = dat_5p95p$paperid, test = "naive-t")

#column(4)
#including dummy preferred
MRA_1_CorrelationCoefficientCorrected_dummyPreferred <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending + Preferred, weights=PrecVarianceSECorrected, data=dat_2p98p)
summary(MRA_1_CorrelationCoefficientCorrected_dummyPreferred)

coef_test(MRA_1_CorrelationCoefficientCorrected_dummyPreferred, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#column (5)
#preferred estimates only
dat_long_preferred <- subset(dat_2p98p, Preferred %in% c('1'))

MRA_1_CorrelationCoefficientCorrected_onlypreferred <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending, weights=PrecVarianceSECorrected, data=dat_long_preferred)
summary(MRA_1_CorrelationCoefficientCorrected_onlypreferred)

coef_test(MRA_1_CorrelationCoefficientCorrected_onlypreferred, vcov = "CR0", 
          cluster = dat_long_preferred$paperid, test = "naive-t")

#column (6)
#without inferior values
dat_long_without_inferior <- subset(dat_2p98p, Preferred %in% c('1', '0'))

MRA_1_CorrelationCoefficientCorrected_noInferior <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending, weights=PrecVarianceSECorrected, data=dat_long_without_inferior)
summary(MRA_1_CorrelationCoefficientCorrected_noInferior )

coef_test(MRA_1_CorrelationCoefficientCorrected_noInferior , vcov = "CR0", 
          cluster = dat_long_without_inferior$paperid, test = "naive-t")

#Table D2

#column (2)
#exclude TotalTaxRevenues
MRA_1_CorrelationCoefficientCorrected_excludeTaxRev <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + GovernmentSpending, weights=PrecVarianceSECorrected, data=dat_2p98p)
summary(MRA_1_CorrelationCoefficientCorrected_excludeTaxRev)

coef_test(MRA_1_CorrelationCoefficientCorrected_excludeTaxRev, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#column (3)
#with study-fixed effects
MRA_1_CorrelationCoefficientCorrected_studyFE <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending + as.factor(paperid) - 1, weights=PrecVarianceSECorrected, data=dat_2p98p)
summary(MRA_1_CorrelationCoefficientCorrected_studyFE)

coef_test(MRA_1_CorrelationCoefficientCorrected_studyFE, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#column (4)
#Random Effects
MRA_1_CorrelationCoefficientCorrected_Random <- plm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending, index = c("paperid","pobsid"), model="random", data=dat_2p98p)
summary(MRA_1_CorrelationCoefficientCorrected_Random)

coef_test(MRA_1_CorrelationCoefficientCorrected_Random, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#column (5)
#OLS
MRA_1_CorrelationCoefficientCorrected_OLS <- lm(CorrelationCoefficientCorrected ~ StandardErrorCorrected + EATR + EMTR + ATR + CorpTaxStructure + DataNonOECDcountries + DataMixofCountries + LongRunExplicit + ShortRunExplicit + TotalTaxRevenues + GovernmentSpending, data=dat_2p98p)
summary(MRA_1_CorrelationCoefficientCorrected_OLS)

coef_test(MRA_1_CorrelationCoefficientCorrected_OLS, vcov = "CR0", 
          cluster = dat_2p98p$paperid, test = "naive-t")

#robust standard errors and t-values
ses_MRA_1_CorrelationCoefficientCorrected_now <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_now, vcov = "CR0", 
                                                                cluster = dat_long$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_now <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_now, vcov = "CR0", 
                                                                  cluster = dat_long$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_now <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_now, vcov = "CR0", 
                                                                  cluster = dat_long$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_1_CorrelationCoefficientCorrected_2p98p <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_2p98p, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_2p98p <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_2p98p, vcov = "CR0", 
                                                                        cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_2p98p <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_2p98p, vcov = "CR0", 
                                                                        cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_pubbias_1_CorrelationCoefficientCorrected_now <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_now, vcov = "CR0", 
                                                                    cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_pubbias_1_CorrelationCoefficientCorrected_now <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_now, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_pubbias_1_CorrelationCoefficientCorrected_now <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_now, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_1_CorrelationCoefficientCorrected_median <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_median, vcov = "CR0", 
                                                                       cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_median <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_median, vcov = "CR0", 
                                                                         cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_median <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_median, vcov = "CR0", 
                                                                         cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_1_CorrelationCoefficientCorrected_OLS <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_OLS, vcov = "CR0", 
                                                                    cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_OLS <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_OLS, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_OLS <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_OLS, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_1_CorrelationCoefficientCorrected_Random <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_Random, vcov = "CR0", 
                                                                   cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_Random <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_Random, vcov = "CR0", 
                                                                     cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_Random <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_Random, vcov = "CR0", 
                                                                     cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_1_CorrelationCoefficientCorrected_partialcorr <- list(coef_test(MRA_1_partialcorr, vcov = "CR0", 
                                                                            cluster = dat$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_partialcorr <- list(coef_test(MRA_1_partialcorr, vcov = "CR0", 
                                                                              cluster = dat$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_partialcorr <- list(coef_test(MRA_1_partialcorr, vcov = "CR0", 
                                                                              cluster = dat$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_1_CorrelationCoefficientCorrected_est_data <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_est_data, vcov = "CR0", 
                                                                   cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_est_data <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_est_data, vcov = "CR0", 
                                                                     cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_est_data <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_est_data, vcov = "CR0", 
                                                                     cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


ses_MRA_1_CorrelationCoefficientCorrected_pubchar <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_pubchar, vcov = "CR0", 
                                                                    cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_pubchar <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_pubchar, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_pubchar <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_pubchar, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


ses_MRA_1_CorrelationCoefficientCorrected_noInferior <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_noInferior, vcov = "CR0", 
                                                                            cluster = dat_long_without_inferior$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_noInferior  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_noInferior, vcov = "CR0", 
                                                                               cluster = dat_long_without_inferior$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_noInferior  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_noInferior, vcov = "CR0", 
                                                                               cluster = dat_long_without_inferior$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_1_CorrelationCoefficientCorrected_onlypreferred <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_onlypreferred, vcov = "CR0", 
                                                                       cluster = dat_long_preferred$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_onlypreferred  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_onlypreferred, vcov = "CR0", 
                                                                          cluster = dat_long_preferred$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_onlypreferred  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_onlypreferred, vcov = "CR0", 
                                                                          cluster = dat_long_preferred$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


ses_MRA_1_CorrelationCoefficientCorrected_studyFE <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_studyFE, vcov = "CR0", 
                                                                          cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_studyFE  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_studyFE, vcov = "CR0", 
                                                                             cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_studyFE  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_studyFE, vcov = "CR0", 
                                                                             cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_1_CorrelationCoefficientCorrected_excludeTaxRev <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_excludeTaxRev, vcov = "CR0", 
                                                                      cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_excludeTaxRev  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_excludeTaxRev, vcov = "CR0", 
                                                                         cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_excludeTaxRev  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_excludeTaxRev, vcov = "CR0", 
                                                                         cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_1_CorrelationCoefficientCorrected_dummyPreferred <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_dummyPreferred, vcov = "CR0", 
                                                                          cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_dummyPreferred  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_dummyPreferred, vcov = "CR0", 
                                                                             cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_dummyPreferred  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_dummyPreferred, vcov = "CR0", 
                                                                             cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_1_CorrelationCoefficientCorrected_2p98p <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_2p98p, vcov = "CR0", 
                                                                           cluster = dat_2p98p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_2p98p  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_2p98p, vcov = "CR0", 
                                                                              cluster = dat_2p98p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_2p98p  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_2p98p, vcov = "CR0", 
                                                                              cluster = dat_2p98p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

ses_MRA_1_CorrelationCoefficientCorrected_5p95p <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_5p95p, vcov = "CR0", 
                                                                  cluster = dat_5p95p$paperid, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_MRA_1_CorrelationCoefficientCorrected_5p95p  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_5p95p, vcov = "CR0", 
                                                                     cluster = dat_5p95p$paperid, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

pvals_MRA_1_CorrelationCoefficientCorrected_5p95p  <- list(coef_test(MRA_1_CorrelationCoefficientCorrected_5p95p, vcov = "CR0", 
                                                                     cluster = dat_5p95p$paperid, test = "naive-t")[,4]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


#Table 3 (main meta-regression table)
stargazer(MRA_1_CorrelationCoefficientCorrected_2p98p, MRA_1_CorrelationCoefficientCorrected_est_data, MRA_1_CorrelationCoefficientCorrected_pubchar, MRA_1_CorrelationCoefficientCorrected_median, MRA_1_partialcorr, t=list(unlist(tvals_MRA_1_CorrelationCoefficientCorrected_2p98p), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_est_data), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_pubchar), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_median), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_partialcorr)), se=list(unlist(ses_MRA_1_CorrelationCoefficientCorrected_2p98p), unlist(ses_MRA_1_CorrelationCoefficientCorrected_est_data), unlist(ses_MRA_1_CorrelationCoefficientCorrected_pubchar), unlist(ses_MRA_1_CorrelationCoefficientCorrected_median), unlist(ses_MRA_1_CorrelationCoefficientCorrected_partialcorr)), p=list(unlist(pvals_MRA_1_CorrelationCoefficientCorrected_2p98p), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_est_data), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_pubchar), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_median), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_partialcorr)))

#Table D1 (see appendix)
stargazer(MRA_1_CorrelationCoefficientCorrected_2p98p, MRA_1_CorrelationCoefficientCorrected_now, MRA_1_CorrelationCoefficientCorrected_5p95p, MRA_1_CorrelationCoefficientCorrected_dummyPreferred, MRA_1_CorrelationCoefficientCorrected_onlypreferred, MRA_1_CorrelationCoefficientCorrected_noInferior, t=list(unlist(tvals_MRA_1_CorrelationCoefficientCorrected_2p98p), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_now), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_5p95p), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_dummyPreferred), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_onlypreferred), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_noInferior)), se=list(unlist(ses_MRA_1_CorrelationCoefficientCorrected_2p98p), unlist(ses_MRA_1_CorrelationCoefficientCorrected_now), unlist(ses_MRA_1_CorrelationCoefficientCorrected_5p95p), unlist(ses_MRA_1_CorrelationCoefficientCorrected_dummyPreferred), unlist(ses_MRA_1_CorrelationCoefficientCorrected_onlypreferred), unlist(ses_MRA_1_CorrelationCoefficientCorrected_noInferior)), p=list(unlist(pvals_MRA_1_CorrelationCoefficientCorrected_2p98p), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_now), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_5p95p), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_dummyPreferred), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_onlypreferred), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_noInferior)))

#Table D2 (see appendix)
stargazer(MRA_1_CorrelationCoefficientCorrected_2p98p, MRA_1_CorrelationCoefficientCorrected_excludeTaxRev, MRA_1_CorrelationCoefficientCorrected_studyFE, MRA_1_CorrelationCoefficientCorrected_Random, MRA_1_CorrelationCoefficientCorrected_OLS, t=list(unlist(tvals_MRA_1_CorrelationCoefficientCorrected_2p98p), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_excludeTaxRev), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_studyFE), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_Random), unlist(tvals_MRA_1_CorrelationCoefficientCorrected_OLS)), se=list(unlist(ses_MRA_1_CorrelationCoefficientCorrected_2p98p), unlist(ses_MRA_1_CorrelationCoefficientCorrected_excludeTaxRev), unlist(ses_MRA_1_CorrelationCoefficientCorrected_studyFE), unlist(ses_MRA_1_CorrelationCoefficientCorrected_Random), unlist(ses_MRA_1_CorrelationCoefficientCorrected_OLS)), p=list(unlist(pvals_MRA_1_CorrelationCoefficientCorrected_2p98p), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_excludeTaxRev), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_studyFE), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_Random), unlist(pvals_MRA_1_CorrelationCoefficientCorrected_OLS)))
