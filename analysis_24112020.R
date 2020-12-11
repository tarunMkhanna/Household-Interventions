library(metafor) #for meta analysis
library(dplyr)
library(ggplot2)
library(texreg) #for regression tables
library(tidyverse) #for read_csv
library(rmarkdown) #for rendering summary to word
library(devtools)
library(clubSandwich)
library(ggrepel) 

setwd("C:/Users/khat/Dropbox/01 Household Interventions/NCC Upload/")

# activate the functions used 
source("C:/Users/khat/Dropbox/01 Household Interventions/NCC Upload/functions_v1.R")

# SPECIFY INPUT FILE HERE
dfinc <- read_csv("C:/Users/khat/Dropbox/01 Household Interventions/NCC Upload/dataset_24112020.csv")

################################################################################################

# CALCULATING STANDARDIZED SIZES AND VARIANCES

# Make effect direction consistent
dfinc <- mutate(dfinc,
                  effect_direction = ifelse(statistical_technique %in% c("Difference of means", "ANOVA"), -1, effect_direction), # Difference of means are always calculated as positive when there is reduction in consumption
                  test_statistic = test_statistic * effect_direction * -1  # Effect direction is neagative, and therefore test statistic is positive, when there is reduction in consumption
            )
                  
# Calculating standardized effect sizes
dfinc <- mutate(dfinc,
          
          # calculating r coefficient for regression results using sample_size
          r_coefficient = ifelse(significance_test == "t-test", ((`test_statistic`)^2/((`test_statistic`)^2 + `total_sample_size`))^(1/2) * effect_direction * -1 ,  
                          ifelse((significance_test == "z-test" | significance_test == "X^2 Test"), (`test_statistic`^2/`total_sample_size`)^(1/2)* effect_direction * -1 , NA)),
                   
          # calculating smd for experiments with difference of means results
          pooled_sd = ifelse(is.na(pooled_sd) & study_design == "Pre-test/post-test", (((`treatment_sample_size`-1)*control_sd^2 + (`treatment_sample_size`-1)*treated_sd^2)/(`treatment_sample_size` + `treatment_sample_size` - 2))^(1/2), 
                      ifelse(is.na(pooled_sd) & study_design != "Pre-test/post-test", (((`control_sample_size`-1)*control_sd^2 + (`treatment_sample_size`-1)*treated_sd^2)/(`control_sample_size` + `treatment_sample_size` - 2))^(1/2),            
                             pooled_sd)),
          
          smd1= ifelse(`statistical_technique` == "Difference of means" | `statistical_technique` == "ANOVA" , `diff_mean`/`pooled_sd`, NA),
                
          smd2= ifelse(((`statistical_technique` == "Difference of means" | `statistical_technique` == "ANOVA") & (significance_test == "F-test" | significance_test == "H-test") & study_design == "Control-treatment"),         
                       (`test_statistic`/`control_sample_size`/`treatment_sample_size`*(`treatment_sample_size` + `control_sample_size`))^(1/2),
                ifelse(((`statistical_technique` == "Difference of means" | `statistical_technique` == "ANOVA") & (significance_test == "F-test" | significance_test == "H-test") & study_design == "Difference in Difference"),  
                       (`test_statistic`/`control_sample_size`/`treatment_sample_size`*(`treatment_sample_size` + `control_sample_size`))^(1/2),
                ifelse(((`statistical_technique` == "Difference of means" | `statistical_technique` == "ANOVA") & (significance_test == "F-test" | significance_test == "H-test") & study_design == "Pre-test/post-test"),        
                       (`test_statistic`/(`treatment_sample_size`)^2*(2*`treatment_sample_size`))^(1/2),   #using treatment sample = control sample     
                              
                ifelse(((`statistical_technique` == "Difference of means" | `statistical_technique` == "ANOVA") & (significance_test == "t-test" | significance_test == "z-test") & study_design == "Control-treatment"),          
                       `test_statistic` * (1/`control_sample_size`/`treatment_sample_size`*(`treatment_sample_size` + `control_sample_size`))^(1/2), 
                ifelse(((`statistical_technique` == "Difference of means" | `statistical_technique` == "ANOVA") & (significance_test == "t-test" | significance_test == "z-test") & study_design == "Difference in Difference"),   
                       `test_statistic` * (1/`control_sample_size`/`treatment_sample_size`*(`treatment_sample_size` + `control_sample_size`))^(1/2), 
                ifelse(((`statistical_technique` == "Difference of means" | `statistical_technique` == "ANOVA") & (significance_test == "t-test" | significance_test == "z-test") & study_design == "Pre-test/post-test"),         
                       `test_statistic` * (1/(`treatment_sample_size`)^2*(2*`treatment_sample_size`))^(1/2), #using treatment sample = control sample
                       
                       NA)))))),
          
          smd = smd1,
          smd = ifelse(is.na(smd1), smd2, smd1),
                             
          # converting d into r for difference of means regressons and ANNOVA
          r_coefficient = ifelse(`statistical_technique`  %in% c("Difference of means", "ANOVA"), smd/(smd^2+4)^(1/2), r_coefficient),

          zi = 0.5 * log((1+r_coefficient) / (1-r_coefficient)),
          
          )

# calculating size variance 
dfinc <- mutate(dfinc, 
          
          # variance of smd         
          var_smd = ifelse(control_sample_size>0, ((`control_sample_size`+`treatment_sample_size`)/(`control_sample_size`*`treatment_sample_size`)) + (smd^2/(2*(`control_sample_size`+`treatment_sample_size`))),
                                                  ((`treatment_sample_size`+`treatment_sample_size`)/(`treatment_sample_size`*`treatment_sample_size`)) + (smd^2/(2*(`treatment_sample_size`+`treatment_sample_size`)))),
          
          #variance of Zs                                                                                                                                           
          vi   = 1/(`total_sample_size` - 3)         
                
          )

################################################################################

dfinc$All <- 1 #to allow calculation for all models together

# CALCULATING AGGREGATE EFFECT SIZE AND COMPARISON ACROSS INTERVENTIONS

models <- c("All")

for (x in models) {
  
  assign(paste0("m.rma.DL.",x),

         rma(yi=zi, vi = vi, method="DL", data = dfinc[dfinc[x] ==1,],  test = "knha" 
         ),
  )

#knha is Knapp and Hartung Method that corrects for the fact that estimates of tau^2 are uncertain. The resulting distributions are t dist instead of Z
#viechtbauer recommends using them as default
  
  assign(paste0("m.rma.REML.",x),
         
         rma(yi=zi, vi = vi, method="REML", data = dfinc[dfinc[x] ==1,],  test = "knha"
         ),
  )

  assign(paste0("m.rma.EB.",x),
         
         rma(yi=zi, vi = vi, method="EB", data = dfinc[dfinc[x] ==1,],  test = "knha"
         ),
  )
    
  assign(paste0("m.rma.mv.",x),

         rma.mv(zi, vi, random = ~ as.factor(InterventionID) | `document ID`, data = dfinc[dfinc[x] ==1,], method = "REML"
                # mods = ~ SocialComparison:Feedback:MonetaryIncentives:Motivation:Information+0
         ),
  )
  
}

# RVE method implemented in metafor
# check using approximating the V matrix and using cluster robust inference methods does not change the estimate

V <- impute_covariance_matrix(dfinc$vi, dfinc$`document ID`, r=0.5, return_list=FALSE)
res <- rma.mv(zi, V, random = ~ as.factor(InterventionID) | `document ID`, data = dfinc, method = "REML")
robust(res, cluster=dfinc$`document ID`)
coef_test(res, vcov="CR2", cluster=dfinc$`document ID`)
conf_int(res, vcov="CR2", cluster=dfinc$`document ID`)

rma.DL.All <- tidy2.rma(m.rma.DL.All)
rma.REML.All <- tidy2.rma(m.rma.REML.All) 
rma.MV.All <- tidy3.rma(m.rma.mv.All)

screenreg(list(rma.DL.All, rma.REML.All, rma.MV.All), 
          custom.model.names = c("RE-DL","RE-REML", "Multilevel"))

# htmlreg(list(rma.DL.All,rma.REML.All, rma.mv.All, vre.All),
#         custom.model.names = c("RE-DL", "RE-REML", "Multilevel"),
#         caption = "Average Treatment Effect",
#         file = "Output/avgeffect.doc",
#         inline.css = FALSE, doctype = TRUE, html.tag = TRUE, head.tag = TRUE, body.tag = TRUE)

# Prediction/Credibility intervals 
predict(m.rma.DL.All)
predict(m.rma.REML.All)
predict(m.rma.mv.All)

# Comparing interventions in a meta regression model and interactions

dfinc$Feedback <- factor(dfinc$Feedback)
dfinc$Information <- factor(dfinc$Information)
dfinc$SocialComparison <- factor(dfinc$SocialComparison)
dfinc$Motivation <- factor(dfinc$Motivation)
dfinc$MonetaryIncentives <- factor(dfinc$MonetaryIncentives)

EffectSizeSummary_MV <- data_frame()

for (x in "All") {
  

  EffectSizeSummary_MV <- rbind(EffectSizeSummary_MV, 
                                tidy3(rma.mv(zi, vi, random = ~ as.factor(InterventionID) | `document ID`, data = dfinc[dfinc[x] ==1,], method = "REML",
                                             mods = ~ SocialComparison:Feedback:MonetaryIncentives:Motivation:Information+0 
                                )
                                      )
                                )
  
  }

EffectSizeSummary_MV<- separate(EffectSizeSummary_MV, model, sep = "[^[:alnum:]]+", into = c("SocialComparison","Feedback","MonetaryIncentives", "Motivation", "Information" ))


EffectSizeSummary_MV <- mutate(EffectSizeSummary_MV,
                               Feedback = as.numeric(recode(Feedback, "Feedback1" = 1, "Feedback0" = 0)),
                               MonetaryIncentives = as.numeric(recode(MonetaryIncentives, "MonetaryIncentives1" = 1, "MonetaryIncentives0" = 0)),
                               Motivation = as.numeric(recode(Motivation, "Motivation1" = 1, "Motivation0" = 0)),
                               Information = as.numeric(recode(Information, "Information1" = 1, "Information0" = 0)),
                               SocialComparison = as.numeric(recode(SocialComparison, "SocialComparison1" = 1, "SocialComparison0" = 0)),
                               All=1
)

EffectSizeSummary_MV <- mutate(EffectSizeSummary_MV,
                               order = case_when(rowSums(EffectSizeSummary_MV[,c( "SocialComparison","Feedback","MonetaryIncentives", "Motivation", "Information" )])==1 ~ 1,
                                                 rowSums(EffectSizeSummary_MV[,c( "SocialComparison","Feedback","MonetaryIncentives", "Motivation", "Information" )])==2 ~ 2,
                                                 rowSums(EffectSizeSummary_MV[,c( "SocialComparison","Feedback","MonetaryIncentives", "Motivation", "Information" )])==3 ~ 3,
                                                 rowSums(EffectSizeSummary_MV[,c( "SocialComparison","Feedback","MonetaryIncentives", "Motivation", "Information" )])==4 ~ 4),
                               label = recode(model1, 
                                              "SocialComparison1:Feedback0:MonetaryIncentives0:Motivation0:Information0" = "Social Comparison", 
                                              "SocialComparison0:Feedback1:MonetaryIncentives0:Motivation0:Information0" = "Feedback",
                                              "SocialComparison0:Feedback0:MonetaryIncentives1:Motivation0:Information0" = "Monetary Incentives",
                                              'SocialComparison1:Feedback1:MonetaryIncentives0:Motivation0:Information0'='Feedback+Social',
                                              'SocialComparison1:Feedback1:MonetaryIncentives0:Motivation0:Information0'='Feedback+Social',
                                              'SocialComparison1:Feedback0:MonetaryIncentives1:Motivation0:Information0'='Social+Monetary',
                                              'SocialComparison0:Feedback1:MonetaryIncentives1:Motivation0:Information0'='Feedback+Monetary',
                                              'SocialComparison1:Feedback1:MonetaryIncentives1:Motivation0:Information0'='Social+Monetary+Feedback',
                                              'SocialComparison0:Feedback0:MonetaryIncentives0:Motivation1:Information0'='Motivation',
                                              'SocialComparison0:Feedback1:MonetaryIncentives0:Motivation1:Information0'='Feedback+Motivation',
                                              'SocialComparison1:Feedback1:MonetaryIncentives0:Motivation1:Information0'='Social+Feedback+Motivation',
                                              'SocialComparison0:Feedback0:MonetaryIncentives1:Motivation1:Information0'='Monetary+Motivation',
                                              'SocialComparison1:Feedback0:MonetaryIncentives1:Motivation1:Information0'='Social+Monetary+Motivation',
                                              'SocialComparison0:Feedback1:MonetaryIncentives1:Motivation1:Information0'='Feedback+Monetary+Motivation',
                                              'SocialComparison0:Feedback0:MonetaryIncentives0:Motivation0:Information1'='Information',
                                              'SocialComparison1:Feedback0:MonetaryIncentives0:Motivation0:Information1'='Social+Information',
                                              'SocialComparison0:Feedback1:MonetaryIncentives0:Motivation0:Information1'='Feedback+Information',
                                              'SocialComparison1:Feedback1:MonetaryIncentives0:Motivation0:Information1'='Social+Feedback+Information',
                                              'SocialComparison0:Feedback0:MonetaryIncentives1:Motivation0:Information1'='Monetary+Information',
                                              'SocialComparison0:Feedback1:MonetaryIncentives1:Motivation0:Information1'='Feedback + Monetary + Information ',
                                              'SocialComparison1:Feedback1:MonetaryIncentives1:Motivation0:Information1'='Social + Feedback + Monetory + Information ',
                                              'SocialComparison0:Feedback0:MonetaryIncentives0:Motivation1:Information1'='Motivation + Information ',
                                              'SocialComparison0:Feedback1:MonetaryIncentives0:Motivation1:Information1'='Feedback + Motivation + Information ',
                                              'SocialComparison0:Feedback1:MonetaryIncentives1:Motivation1:Information1'='Feedback + Motivation + Monetary + Information ')

                               
)

# Ploting the effect sizes and CI

ggplot(EffectSizeSummary_MV[EffectSizeSummary_MV$order==1,], #change data frame to get the correct graph
  aes(x = reorder(label,beta), y = beta)) +
  scale_y_continuous(limits = c(-0.05,.5))+
  geom_point(size = 2) +
  geom_linerange(aes(ymin = lower_lim, ymax = upper_lim),
                 size = .5, alpha = 1, colour = "red")  + 
  labs(x = "", y = "") +
  geom_hline(yintercept = m.rma.mv.All$beta, color = "grey") +
  theme_bw()

color <- c('Social+Monetary+Feedback', 'Feedback+Monetary+Motivation','Social+Feedback+Motivation', 'Feedback+Social',
           'Feedback+Motivation')

  ggplot(EffectSizeSummary_MV,
  aes(x =order, y = beta, color = label %in% color))+
  guides(color = FALSE, legend.position="bottom") +
  scale_y_continuous(limits = c(-0.05,.5))+
  geom_point(size =2) +
  geom_text_repel(data = subset(EffectSizeSummary_MV, pvalues < 0.05), aes(label=label),  size = 5, hjust=0, vjust=0)+
  labs(x = "Number of interventions used", y = "") +
  theme_bw()+
  scale_color_brewer(palette="Dark2")+  
  theme(legend.position="bottom")
  


################################################################################################

# DATA CLEANING AND CONTROLS

dfinc$study_design <- as.factor(dfinc$study_design)

dfinc <- mutate(dfinc, 
             ElectricityUse = `controls Electricity use`,
             ElectricityPrice = `controls Energy prices`,
             EnvironmentAttitudes =  `controls Environmental attitudes`,
             HouseholdType =  `controls HH controls (demographics)`,
             ResidenceType = `controls Residence controls`,
             Season =   `controls Seasonal controls`,
             Weather = `controls Weather controls`,
             
             Region = case_when(
               country %in% c("China") ~ 'Asia',
               country %in% c("India") ~ 'Asia',
               country %in% c("Japan") ~ 'Asia',
               country %in% c("South Korea") ~ 'Asia',
               country %in% c("Singapore") ~ 'Asia',
               
               country %in% c("Austria") ~ 'Europe excl. UK',
               country %in% c("Denmark") ~ 'Europe excl. UK',
               country %in% c("Finland") ~ 'Europe excl. UK',
               country %in% c("France") ~ 'Europe excl. UK',
               country %in% c("Germany") ~ 'Europe excl. UK',
               country %in% c("Italy") ~ 'Europe excl. UK',
               country %in% c("Latvia") ~ 'Europe excl. UK',
               country %in% c("Netherlands") ~ 'Europe excl. UK',
               country %in% c("Northern") ~ 'Europe excl. UK',
               country %in% c("Poland") ~ 'Europe excl. UK',
               country %in% c("Sweden") ~ 'Europe excl. UK',
               country %in% c("Switzerland") ~ 'Europe excl. UK',
               country %in% c("Ireland") ~ 'Europe excl. UK',
               country %in% c("Portugal") ~ 'Europe excl. UK',
               country %in% c("Spain") ~ 'Europe excl. UK',
               
               country %in% c("Australia") ~ 'Others',
               country %in% c("Colombia") ~ 'Others',
               country %in% c("Ecuador") ~ 'Others',
               country %in% c("Israel") ~ 'Others',
               country %in% c("South Africa") ~ 'Others',
               country %in% c("Canada") ~ 'Others',
               country %in% c("Mexico") ~ 'Others',
               country %in% c("Iran") ~ 'Others',
               
               country %in% c("United Kingdom") ~ 'United Kingdom',
               country %in% c("United States") ~ 'United States',
               
               TRUE ~ country),
             
             StatsMethod = case_when( statistical_technique %in% c("Difference of means", "ANOVA") ~ "MeansDifferences",
                                      statistical_technique %in% c("Household Fixed effects regression", "Time & Household Fixed effects regression", "Fixed effects regression",
                                                                   "Random effects regression") ~ "PanelEffects", 
                                      statistical_technique %in% c("OLS regression") ~ "OLS",
                                      TRUE ~ "Others"),
             
             StudyDesign = case_when( study_design  %in% c("Control-treatment") ~ "Control-treatment",
                                      study_design  %in% c("Difference in Difference") ~ "DID",
                                      study_design  %in% c("Pre-test/post-test") ~ "PrePost"),
             
             Control = if_else(control_sample_size>0, 1, 
                               if_else(is.na(control_sample_size), 0, 0),0),
             
             BaselineConsumption = `Baseline Consumption`, 
             InterventionDuration = interventionDuration,
             
             Randomisation = case_when(`Randomisation method`  %in% c("no", "No", "none", "None", "not done","Not randoized", "not randomized", "Not randomized", 
                                                                            "Not randomized, lower socio economic status population") ~ "No",
                                             # `Randomisation method`  %in% c(NA, "No details") ~ "NA",
                                             TRUE ~ "Yes"),  
             
             OptedIn = case_when(Opt_in == 1 ~ "Yes",
                                TRUE ~ "No"),
             
             After2000 = case_when(`document PY` > 2000 ~ "After",
                                   TRUE ~ "Before"),
             
             StudyYear = `document PY`,
             
             Followup = case_when(interventionFollowup>0 ~ 1,
                                  interventionFollowup == NA ~ 0, 
                                  TRUE ~ 0)
)

dfinc$StudyDesign <- factor(dfinc$StudyDesign, levels = c("DID", "Control-treatment", "PrePost"))
dfinc$StudyDesign <- relevel(dfinc$StudyDesign, ref =  "PrePost")

dfinc$StatsMethod <- as.factor(dfinc$StatsMethod)
dfinc$StatsMethod <- factor(dfinc$StatsMethod, levels = c("PanelEffects", "PanelEffects", "OLS", "Others"))
dfinc$StatsMethod <- relevel(dfinc$StatsMethod, ref = "MeansDifferences")

dfinc$Region <- as.factor(dfinc$Region)
dfinc$Region <- factor(dfinc$Region, levels = c("Asia", "United States", "United Kingdom","Europe excl. UK", "Others"))
dfinc$Region <- relevel(dfinc$Region, ref = "United States")


################################################################################

# META REGRESSION MODELS

for (x in c("All", "Feedback","Information", "SocialComparison", "Motivation", "MonetaryIncentives")) {

      assign(paste0("rma.",x),

                 tidy2.reg(rma(yi=zi, vi = vi, method="REML", data = dfinc[dfinc[x] ==1,], test = "knha"
                           ,mods = ~ 
                            
                             StudyDesign + Randomisation + 
                             HouseholdType + ResidenceType  + Weather + #Controls used in papers with base as no controls 
                             OptedIn + interventionTreatmentPeriod +
                             Region + StudyYear 

                           
                           )),
      )
}

#reg type = Study level cluster
for (x in c("All","Feedback","Information", "SocialComparison", "Motivation", "MonetaryIncentives")) {

  assign(paste0("rma.mv.",x),

         tidy3.reg(rma.mv(zi, vi, random = ~ as.factor(InterventionID) | `document ID`, data = dfinc[dfinc[x] ==1,], method = "REML",
                        mods = ~
                          
                          StudyDesign + Randomisation + 
                          HouseholdType + ResidenceType  + Weather + 
                          OptedIn + interventionTreatmentPeriod +
                          Region + StudyYear 

         )),
  )
}



# print regression outputs to console 

screenreg(list(rma.All, rma.Feedback, rma.Information, rma.MonetaryIncentives, rma.Motivation, rma.SocialComparison),
          custom.model.names = c("All", "Feedback","Information", "MonetaryIncentives","Motivation",  "SocialComparison"),
          caption = "Method: RE-REML")

screenreg(list(rma.mv.All, rma.mv.Feedback, rma.mv.Information, rma.mv.MonetaryIncentives, rma.mv.Motivation, rma.mv.SocialComparison),
          custom.model.names = c("All", "Feedback","Information", "MonetaryIncentives","Motivation",  "SocialComparison"),
          caption = "Method: Multilevel")


# output to word
# htmlreg(list(rma.All, rma.Feedback, rma.Information, rma.MonetaryIncentives, rma.Motivation, rma.SocialComparison),
#         custom.model.names = c("All", "Feedback","Information", "MonetaryIncentives","Motivation",  "SocialComparison"),
#         caption = "Method: RE-REML",
#         file = "Output/MetaRegression - REML.doc",
#         inline.css = FALSE, doctype = TRUE, html.tag = TRUE, head.tag = TRUE, body.tag = TRUE
# )
# 
# htmlreg(list(rma.mv.All, rma.mv.Feedback, rma.mv.Information, rma.mv.MonetaryIncentives, rma.mv.Motivation, rma.mv.SocialComparison),
#         custom.model.names = c("All", "Feedback","Information", "MonetaryIncentives","Motivation",  "SocialComparison"),
#         caption = "Method: Multilevel",
#         file = "Output/MetaRegression - MV.doc",
#         inline.css = FALSE, doctype = TRUE, html.tag = TRUE, head.tag = TRUE, body.tag = TRUE
# )


################################################################################################

# Funnel Plots for small-study bias

funnel(m.rma.REML.All, xlab="estimate", yaxis= "sei", cex=.5, ylim=(c(.5, 0)),
       level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), refline=0, legend=FALSE)

# Egger's test
regtest(m.rma.REML.All)

#Trim and fill method to do sensitivity of results due to publication bias
sensitivity <- trimfill(m.rma.REML.All) 

################################################################################################

# DESCRIPTIVE STATISTICS

# https://cran.r-project.org/web/packages/summarytools/vignettes/Introduction.html
# https://dabblingwithdata.wordpress.com/2018/01/02/my-favourite-r-package-for-summarising-data/  

# Plotting the increase in literature over time (all papers for which full text screening was done)
dfinc %>% filter(`document PY` >0) %>% distinct(`document ID`, `document PY`) %>%
  ggplot( aes(x = `document PY`)) +
  #geom_line(stat = "density") +
  geom_histogram(binwidth = 1)+
  labs(title = "Literature related to Household
             Energy Interventions",
       x = "Publication Year",
       y = "Number of Documents") +
  theme_minimal()


sum1 <- ggplot(dfinc, #change data frame to get the correct graph
  ) +
  geom_point(aes(y = zi, x = log(total_sample_size), color = "total_sample_size") , size =1) +
  geom_point(aes(y = zi, x = log(control_sample_size), color = "control_sample_size") , size =1) +
  labs(y = "Average Effect", x = "Log Sample Size") +
  theme_bw() + 
  scale_fill_continuous(guide = guide_legend()) + theme(legend.position="bottom")

sum2<- ggplot(dfinc, #change data frame to get the correct graph
) +
  geom_point(aes(y = zi, x = dfinc$interventionDuration, color = "interventionDuration") , size =1) +
  geom_point(aes(y = zi, x = dfinc$interventionTreatmentPeriod, color = "interventionTreatmentPeriod") , size =1) +
  labs(y = "Average Effect", x = "Duration of Intervention in Weeks") +
  theme_bw()+ 
  scale_fill_continuous(guide = guide_legend()) + theme(legend.position="bottom")


sum3<- ggplot(dfinc) +
  geom_point(aes(y = zi, x = dfinc$`document PY`, color = StudyDesign) , size =1) +
  labs(y = "Average Effect", x = "Year") +
  theme_bw()+ 
  scale_fill_continuous(guide = guide_legend()) + theme(legend.position="bottom")

dfinclong <- pivot_longer(dfinc, MonetaryIncentives:Motivation, names_to = "intervention", values_to = "response")

sum4<-
  ggplot(dfinclong,) +
  geom_point(aes(y = zi, x = `document PY`, color = intervention) , size =1) +
  facet_wrap(~intervention)  +
  labs(y = "Average Effect", x = "Year") +
  theme_bw()+ 
  scale_fill_continuous(guide = guide_legend()) + theme(legend.position="bottom")

