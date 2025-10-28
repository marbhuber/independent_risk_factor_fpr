##########################################################
##########################################################
# Independent risk factors and the risk of false positives
#
# 28 October 2025
# Markus Huber (markus.huber@insel.ch)
##########################################################
##########################################################

#######
# setup
#######

rm(list=ls())
library(tidyverse)
library(HDSinRdata)
library(pCalibrate)
library(flextable)
library(compareGroups)

######
# data
######

data("pain")

df <- pain %>% 
  transmute(
    outcome         = case_when(PAIN_INTENSITY_AVERAGE.FOLLOW_UP>=7~1,TRUE~0),
    pre             = case_when(PAIN_INTENSITY_AVERAGE>=7~1,TRUE~0),
    AGE_AT_CONTACT,
    MEDICAID_BIN    = factor(MEDICAID_BIN),
    PAT_SEX         = factor(PAT_SEX),
    BMI,
    CCI_BIN         = factor(CCI_BIN,levels = c("No comorbidity","Any comorbidity"))) %>%
  na.omit()

##########################################################################
# Illustrate FPR with a common P-value between 0.040 and 0.045 for CCI_BIN
# Use only 1'000 patients for illustration purposes
##########################################################################

stop = 0
for (myseed in 1:1e7){
  if (stop==0){
    set.seed(myseed)
    rm(df.sample,mod.sample,p.sample)
    df.sample <- df %>% sample_n(1000)
    p.sample  <- summary(glm(outcome~.,data=df.sample,family=binomial()))$coefficients["CCI_BINAny comorbidity","Pr(>|z|)"]
    if (p.sample>0.040 & p.sample<0.045){stop=1}
  }
}

# fit the final logistic regression
mymod = glm(outcome~.,data=df.sample,family=binomial())

# export model summary
save_as_docx(
  gtsummary::tbl_regression(mymod, exponentiate = T) %>% 
    gtsummary::as_flex_table(),
  path = "table1_model.docx")

# export summary measures
export2word(
  createTable(
    compareGroups(~.,data=df.sample,method=NA)
    ),
  file = "table1_summary.docx")

# get minimum Bayes Factors
pCalibrate(summary(mymod)$coefficients[,"Pr(>|z|)"]) %>% formatBF()
pCalibrate(summary(mymod)$coefficients[,"Pr(>|z|)"]) %>% round(3)
BF01.CCI =  pCalibrate(summary(mymod)$coefficients[,"Pr(>|z|)"])["CCI_BINAny comorbidity"] %>% as.numeric()
BF01.pre =  pCalibrate(summary(mymod)$coefficients[,"Pr(>|z|)"])["pre"] %>% as.numeric()

#########################################
# Calculate the FPR as function of Pr(H1)
#########################################

df.fpr <- data.frame()
for (pr.H1 in seq(0,1,by=0.001)){
  # CCI
  rm(myvar)
  myvar = BF01.CCI*(1-pr.H1)/pr.H1
  df.fpr <- rbind(
    df.fpr,
    data.frame(pr.H1,fpr = myvar/(1+myvar),type = "Comorbidities"))
  # Prepain
  rm(myvar)
  myvar = BF01.pre*(1-pr.H1)/pr.H1
  df.fpr <- rbind(
    df.fpr,
    data.frame(pr.H1,fpr = myvar/(1+myvar),type = "Current severe pain (NRS \u2265 7)"))
}
df.fpr <- df.fpr %>% na.omit()

df.fpr.cci        = filter(df.fpr,type == "Comorbidities")
reverse.cci       = which(abs(df.fpr.cci$fpr - 0.05) == min(abs(df.fpr.cci$fpr - 0.05)))
pr.H1.reverse.cci = df.fpr.cci$pr.H1[reverse.cci]
  
df.reverse.cci = data.frame(
  type  = "Comorbidities",
  fpr   = 0.05,
  pr.H1 = pr.H1.reverse.cci
)
df.equal.cci = data.frame(
  type  = "Comorbidities",
  fpr   = df.fpr.cci %>% filter(pr.H1==0.5) %>% pull(fpr),
  pr.H1 = 0.5
)

mycol = ggsci::pal_jama()(2)

# Figure
df.fpr %>% 
  ggplot(aes(x=pr.H1,y=fpr,color = type,fill=type))+
  geom_hline(yintercept = 0.05,linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = 0.5,linetype="dashed",color="darkgrey")+
  geom_line(size=1)+
  geom_point(data=df.reverse.cci,size=4,show.legend = F,shape=23, fill=mycol[1])+
  geom_point(data=df.equal.cci,size=4,show.legend = F)+
  theme_classic()+
  theme(
    text = element_text(size=18),
    legend.title = element_blank(),
    legend.position = c(0.75,0.8)
  )+
  ggsci::scale_color_jama()+
  scale_x_continuous(n.breaks = 10,labels = scales::percent, expand = c(0.02, 0.02))+
  scale_y_continuous(n.breaks = 10,labels = scales::percent, expand = c(0.02, 0.02))+
  xlab("Probability of a real effect: P(H1)")+
  ylab("False positive risk (FPR)")

ggsave("figure1.jpeg",dpi=600,width = 7.4,height=4.7)
