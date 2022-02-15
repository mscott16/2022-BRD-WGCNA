#statistical testing
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(car)


stat17 = read.csv("2017_stats.csv")

shapiro.test(stat17$Neuper)
bartlett.test(Neuper~Class, data = stat17)

shapiro.test(stat17$Eosper)
bartlett.test(Eosper~Class, data = stat17)

shapiro.test(stat17$Basoper)
bartlett.test(Basoper~Class, data = stat17)

shapiro.test(stat17$Monoper)
bartlett.test(Monoper~Class, data = stat17)

shapiro.test(stat17$Lymphper)
bartlett.test(Lymphper~Class, data = stat17)

shapiro.test(stat17$NL_Ratio)
bartlett.test(NL_Ratio~Class, data = stat17)

shapiro.test(stat17$WBC)
bartlett.test(WBC~Class, data = stat17)

shapiro.test(stat17$RBC)
bartlett.test(RBC~Class, data = stat17)

shapiro.test(tat17$HGB)
bartlett.test(HGB~Class, data = stat17)

shapiro.test(stat17$HCT)
bartlett.test(HCT~Class, data = stat17)

shapiro.test(stat17$MCV)
bartlett.test(MCV~Class, data = stat17)

shapiro.test(stat17$MCH)
bartlett.test(MCH~Class, data = stat17)

shapiro.test(stat17$MCHC)
bartlett.test(MCHC~Class, data = stat17)

shapiro.test(stat17$RDW)
bartlett.test(RDW~Class, data = stat17)

shapiro.test(stat17$Platelets)
bartlett.test(Platelets~Class, data = stat17)

shapiro.test(stat17$FEC0)
bartlett.test(FEC0~Class, data = stat17)

shapiro.test(stat17$ADG12)
bartlett.test(ADG12~Class, data = stat17)

shapiro.test(stat17$ADG26)
bartlett.test(ADG26~Class, data = stat17)

shapiro.test(stat17$ADG82)
bartlett.test(ADG82~Class, data = stat17)

shapiro.test(stat17$Growth_Rate)
bartlett.test(Growth_Rate~Class, data = stat17)

shapiro.test(stat17$Sex)
bartlett.test(Sex~Class, data = stat17)

shapiro.test(stat17$Temp_D0)
bartlett.test(Temp_D0~Class, data = stat17)


#T-test (normal, continuous)
t.test(WBC~Class, data = stat17, var.equal = TRUE)

#Chisquare (categorical)
chisq.test(stat17$Class, stat17$Sex)

#T-test (non-normal, continuous)
t.test(RDW~Class, data = stat17, var.equal = FALSE)

probes=names(datExpr0)
modules=paste(probes, bwnet2$colors, moduleColorsBRD, sep=",")
write.csv(modules, file="colormods.csv", sep=",", quote = FALSE, row.names = FALSE)
