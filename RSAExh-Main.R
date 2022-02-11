library(tidyverse)
library(future.apply) # to run lapply in parallel
library(scales) # for percentages in graph scales
library(xtable) # to export latex tables
library(ggnewscale) # to have two different fill scales in production plots

plan(multisession) ## Run in parallel on local computer

options(tibble.width=Inf)

# Custom function to map sample size to alpha
custom_log_trans <- function()trans_new("custom_log",function(x){1+log(x)},function(x){exp(x-1)})

# Standard error function for graphs
se <- function(x){sd(x)/sqrt(length(x))}

# Load the speaker and listener functions for each model
source("RSAexh-ModelFunctions.R")

# Load the MLE functions:
source("RSAExh-MLE-functions.R")


##########################
# Load and clean the data
##########################

production_data <- rbind(read.csv("Speaker-98-anon.csv"),read.csv("Speaker-99-aborted-anon.csv"))
comprehension_data <- rbind(read.csv("Listener-98-anon.csv"),read.csv("Listener-99-aborted-anon.csv"))

# Code answers to the comprehension question:
production_data$Comprehension <- NA
comprehension_data$Comprehension <- NA
# Ugly loop, but it works:
for(i in 1:nrow(production_data)){
  production_data$Comprehension[i] <- production_data[i,8+20-production_data$Input.ApproxPrice[i]]=="true"
}
for(i in 1:nrow(comprehension_data)){
  comprehension_data$Comprehension[i] <- comprehension_data[i,8+20-comprehension_data$Input.ApproxPrice[i]]=="true"
}
# A bit more cleanup:
production_data <- production_data %>%
  select(Participant,WorkTimeInSeconds,Input.ApproxPrice,Input.World,Answer.prior,Answer.utterance,Message,Comprehension) %>%
  rename(Price=Input.ApproxPrice,
         World=Input.World,
         Prior=Answer.prior,
         Utterance=Answer.utterance) %>%
  mutate(World = if_else(World=="only the sandwich","a","ab"),
         Utt.Length = str_length(Utterance))

comprehension_data <- comprehension_data %>%
  select(Participant,WorkTimeInSeconds,Input.ApproxPrice,Input.Utterance,Answer.prior,Answer.posterior,Comprehension) %>%
  rename(Price=Input.ApproxPrice,
         Utterance=Input.Utterance,
         Prior=Answer.prior,
         Posterior=Answer.posterior)

# Remove participants who failed control question and compress priors:
clean_production_data <- production_data %>%
  filter(Comprehension&!is.na(Message)) %>%
  mutate(Prior=Prior/100,
         Prior=.005+.99*Prior)

clean_comprehension_data <- comprehension_data %>%
  filter(Comprehension) %>%
  mutate(Prior=Prior/100,
         Prior=.005+.99*Prior,
         Posterior=Posterior/100,
         Utterance=if_else(Utterance=="The sandwich","a","a&b")
  )

######################
# Priors distribution
######################

# Merge the data from production and comprehension, as the surveys are identical until participants submit their prior response:
prior_data <- bind_rows(
  clean_comprehension_data  %>% mutate(Participant=as.character(Participant)) %>% select(Participant,Price,Prior),
  clean_production_data  %>% mutate(Participant=as.character(Participant)) %>% select(Participant,Price,Prior)
)

mean(prior_data$Prior);sd(prior_data$Prior);

# Test the correlation between prior and price: no surprise, our manipulation worked.
PricePriorCorrel <- lm(Prior~Price,data=prior_data)
summary(PricePriorCorrel)
anova(PricePriorCorrel)

pdf(file="Prior_Price.pdf",height=5,width=7)
# Boxplot of prior by price with regression line superimposed:
prior_data %>% ggplot(aes(x=Price,y=Prior,group=Price))+
  geom_boxplot(color="blue")+
  ylab(expression("Prior"~P(w[ab])))+
  geom_smooth(aes(group=NULL),method=lm,color='red',fill="red",fullrange=T)+
  #geom_abline(intercept=PricePriorCorrel$coefficients[1],slope=PricePriorCorrel$coefficients[2],color="red",size=1)+
  scale_x_continuous(breaks=unique(prior_data$Price),name = "Price of the lunch deal ($)",limits=c(10,22))+
  scale_y_continuous(limits = c(0,1),labels = scales::percent)+
  coord_cartesian(xlim=c(11.5,19.5))+
  theme_bw()+
  theme(
    panel.grid.minor=element_blank(),
    axis.title = element_text(size=16),
    axis.text = element_text(size=14)
  )
dev.off()

################################
# Effect of prior on production
################################

# What is the mean length of the expression used to convey each message?
production_data_summary <- production_data %>%
  filter(Comprehension,!is.na(Message)) %>%
  mutate(Message=factor(Message,levels=c("a","a&b","a&!b","b","!b"))) %>%
  group_by(Message) %>%
  summarize(Mean=mean(Utt.Length),SD=sd(Utt.Length),SE=se(Utt.Length),n=n())

print(production_data_summary)

# Boxplot of utterance length by message:
pdf(file="Length_Message.pdf",height=5,width=7)
production_data %>%
  filter(!is.na(Message)) %>%
  mutate(Message = factor(Message,levels=c("a","a&b","a&!b","b","!b"),labels=c("A","A&B","A&!B","B","!B"))) %>%
  ggplot(aes(x=Message,y=Utt.Length)) +
  geom_boxplot(color="blue") + 
  geom_hline(yintercept=c(1,50),linetype=2) +
  scale_y_continuous(limits=c(0,50),name="Length (#char)")+
  theme_bw()+
  theme(
    axis.title = element_text(size=16),
    axis.text = element_text(size=14)
  )
dev.off()  


# Are participants more likely to give an exhaustive response when the prior for w_ab is higher?
prod_regression <- clean_production_data %>%
  filter(World=="a"&Message%in%c("a&!b","!b","a")) %>%
  mutate(Exhaustive = Message%in%c("!b","a&!b"),
         NormPrior = as.vector(scale(Prior))) %>%
  glm(Exhaustive~NormPrior,family=binomial(),data=.)
summary(prod_regression)
print(paste("chi2(1) = ",round(anova(prod_regression)[2,2],2),", p = ",round(pchisq(anova(prod_regression)[2,2],1,lower.tail = F),3),sep=""))


production_regression_data <- clean_production_data %>%
  filter((World=="a"&Message%in%c("a&!b","!b","a"))|World=="ab") %>%
  mutate(TargetMessage = (World=="a"&!Message%in%c("!b","a&!b"))|(World=="ab"&Message=="a&b"),
         World = factor(World,levels=c("a","ab"))
         )

# Plot the logit regression on the data:
production_plot_data <- clean_production_data %>%
  filter((World=="ab"&Message%in%c("a&b","a","b"))|(World=="a"&Message%in%c("a&!b","!b","a"))) %>%
  mutate(Prior=100*(Prior-0.005)/.99,
         bin=(Prior%/%25),
         Message=if_else(Message=="!b","a&!b",as.character(Message))) %>%
  group_by(World,bin,Message) %>%
  summarise(n=n()) %>%
  mutate(freq=n/sum(n),
         x=if_else(bin==4,100,25*bin+12),
         width=if_else(bin==4,1,25),
         alpha=log(n)+1,
         alpha=alpha/(log(19)+1)) %>%
  filter((Message=="a"&World=="a")|(Message=="a&b"&World=="ab"))
production_plot_data$World <- factor(production_plot_data$World,levels=c("a","ab"))

speaker_labels <- c(expression(S*"(A|"*w[a]*")"),expression(S*"(A&B|"*w[ab]*")"))
levels(production_plot_data$World) <-speaker_labels 
levels(production_regression_data$World) <-speaker_labels 


pdf(file="Production_Prior.pdf",width=7,height=5)
ggplot(data=production_plot_data,aes(x=x/100,y=freq,group=World,fill=n)) +
  facet_grid(.~World,labeller = label_parsed)+
  geom_point(data=filter(production_plot_data,x==100),
             aes(x=x/100,y=freq,fill=n,group=Message),
             inherit.aes = F,size=8,color="black",pch=21) +
  geom_bar(data=production_plot_data,
           aes(x=x/100,y=freq,fill=n,group=Message),width=production_plot_data$width/100,
           inherit.aes = F,stat="identity",color="transparent") +
  scale_fill_gradient(low=rgb(0,0,0,.25),high=rgb(0,0,0,1),trans="custom_log",breaks=c(2,5,15,25))+
  new_scale_fill()+
  geom_smooth(data=production_regression_data,aes(x=Prior,y=as.numeric(TargetMessage),group=World,col=World,fill=World),method="glm",method.args = list(family = "binomial"),inherit.aes = F,show.legend = F)+
  scale_x_continuous(name=expression("Prior"~P(w[ab])),labels=scales::percent)+
  scale_y_continuous(name="Frequency",labels=scales::percent)+
  scale_color_manual(values=c(rgb(.7,0,.7),rgb(1,.7,0)),labels=speaker_labels)+
  scale_fill_manual(values=c(rgb(.7,0,.7),"transparent"),labels=speaker_labels)+
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=16),
    axis.text = element_text(size=12),
    strip.text = element_text(size=16)
  )
dev.off()


###################################
# Effect of prior on comprehension
###################################

# Linear regression Posterior by Prior:
comp_regression_a <- clean_comprehension_data %>%
  filter(Utterance=="a") %>%
  lm(Posterior~Prior,data=.)

comp_regression_ab <- clean_comprehension_data %>%
  filter(Utterance=="a&b") %>%
  lm(Posterior~Prior,data=.)

print(paste("beta = ",round(coefficients(comp_regression_a)[["Prior"]],2),", F(1,",anova(comp_regression_a)[2,1],") = ",round(anova(comp_regression_a)[1,4],1),", p = ",round(anova(comp_regression_a)[1,5],3),sep=""))
print(paste("beta = ",round(coefficients(comp_regression_ab)[["Prior"]],2),", F(1,",anova(comp_regression_ab)[2,1],") = ",round(anova(comp_regression_ab)[1,4],1),", p = ",round(anova(comp_regression_ab)[1,5],3),sep=""))

# Organize the regression results for plotting:
Correls <- tibble(intercept=c(comp_regression_a$coefficients[1],comp_regression_ab$coefficients[1]),
                  slope=c(comp_regression_a$coefficients[2],comp_regression_ab$coefficients[2]),
                  Utterance=c("A","A&B"))

# Graph with correlation lines for each utterance
pdf(file="Comprehension_Prior.pdf",height=5,width=7)
clean_comprehension_data %>%
  mutate(Utterance=factor(Utterance,levels=c("a","a&b"),labels=c("A","A&B"))) %>%
  ggplot(aes(x=Prior,y=Posterior,col=Utterance,fill=Utterance))+
  geom_point()+
  #geom_abline(data=Correls,aes(slope=slope,intercept=intercept,col=Utterance))+
  geom_smooth(method=lm,fullrange=T)+
  geom_abline(slope=1,intercept=0,linetype=2)+
  scale_x_continuous(name=expression("Prior"~P(w[ab])),labels=scales::percent,limits=c(0,1))+
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(name=expression("Posterior"~P(w[ab]*"|"*u)),labels=scales::percent)+
  scale_color_manual(values=c(rgb(.7,0,.7),rgb(1,.7,0)),aesthetics = c("fill","color"))+
  theme_bw()+
  theme(
    axis.title = element_text(size=16),
    axis.text = element_text(size=12),
    strip.text = element_text(size=16)
  )
dev.off()

# Logistic regression predicting the probability to observe anti-exhaustivity:
anti_exh_model <- clean_comprehension_data %>%
  as_tibble %>%
  filter(Utterance=="a") %>%
  mutate(AntiExh=((.99*Posterior+.005)>Prior)) %>%
  glm(AntiExh~scale(Prior),data=.,family=binomial())
summary(anti_exh_model)


##########################################
# Testing the Rational speaker assumption
##########################################

asin_comprehension_data <- clean_comprehension_data %>%
  filter(Utterance=="a") %>%
  mutate(Prior=asin(sqrt(Prior)))

EmpiricalListener_a_1 <-  lm(asin(sqrt(Posterior))~Prior,data=asin_comprehension_data)
EmpiricalListener_a_2 <- lm(asin(sqrt(Posterior))~Prior+I(Prior^2),data=asin_comprehension_data)
EmpiricalListener_a_3 <- lm(asin(sqrt(Posterior))~Prior+I(Prior^2)+I(Prior^3),data=asin_comprehension_data)
EmpiricalListener_a_4 <- lm(asin(sqrt(Posterior))~Prior+I(Prior^2)+I(Prior^3)+I(Prior^4),data=asin_comprehension_data)
EmpiricalListener_a_5 <- lm(asin(sqrt(Posterior))~Prior+I(Prior^2)+I(Prior^3)+I(Prior^4)+I(Prior^5),data=asin_comprehension_data)

BIC(EmpiricalListener_a_1,EmpiricalListener_a_2,EmpiricalListener_a_3,EmpiricalListener_a_4,EmpiricalListener_a_5)

# Degree 1 is the best by BIC:
EmpiricalListener <- function(p){
  sin(predict.lm(EmpiricalListener_a_1,newdata=data.frame(Prior=asin(sqrt(p))),type="response"))^2
}

# # Plot if curious:
# plot(jitter(Posterior,amount=.0025)~jitter(sin(Prior)^2,amount=.0025),data=asin_comprehension_data,xlim=c(0,1),ylim=c(0,1),pch=16,col=rgb(1,.7,0,.5),las=1,ylab="Posterior",xlab="Prior")
# curve(EmpiricalListener(x),add=T,n=10000,col=rgb(.4,0,.7),lwd=2)

# Define the RSA speaker for S(A|wa) and S(A|wab):
RSAspeaker_a <- function(p,lambda,LdeltaAnB){ # S("a"|w_a)
  EL=1-EmpiricalListener(p)
  1/(1+exp(-lambda*log(EL)-LdeltaAnB)) # Note: we're using the product lambda*Delta as our second variable
}
RSAspeaker_ab <- function(p,lambda,LdeltaAB){ # S("a"|w_ab)
  EL=EmpiricalListener(p)
  1/(1+exp(-lambda*log(EL)-LdeltaAB)) # Note: we're using the product lambda*Delta as our second variable
}

# Define the neg-log-likelihood functions to minimize:
neg_log_lik = function(par){ # par = c(log(lambda),LDeltaAB,LDeltaAnB,qlogis(eps))
  lambda = exp(par[1])
  eps=plogis(par[4])
  tmp <- clean_production_data %>%
    mutate(Pred=if_else(World=="a",RSAspeaker_a(Prior,lambda,par[3]),RSAspeaker_ab(Prior,lambda,par[2])),
           loglik=log(case_when(
             Message=="a" ~ (Pred+eps)/(1+3*eps),
             World=="a"&Message%in%c("a&!b","!b") ~ (1-Pred+eps)/(1+3*eps),
             World=="ab"&Message=="a&b" ~ (1-Pred+eps)/(1+3*eps),
             T ~ eps/(1+3*eps)
             ))
           )
  return(-sum(tmp$loglik))
}
neg_log_lik_PQ = function(par){ # par = c(log(lambda),log(LDeltaAB),log(LDeltaAnB),qlogis(eps))
  lambda = exp(par[1])
  par[2:3] <- exp(par[2:3])
  eps=plogis(par[4])
  tmp <- clean_production_data %>%
    mutate(Pred=if_else(World=="a",RSAspeaker_a(Prior,lambda,par[3]),RSAspeaker_ab(Prior,lambda,par[2])),
           loglik=log(case_when(
             Message=="a" ~ (Pred+eps)/(1+3*eps),
             World=="a"&Message%in%c("a&!b","!b") ~ (1-Pred+eps)/(1+3*eps),
             World=="ab"&Message=="a&b" ~ (1-Pred+eps)/(1+3*eps),
             T ~ eps/(1+3*eps)
           ))
    )
  return(-sum(tmp$loglik))
}
neg_log_lik_EQ = function(par){ # par = c(log(lambda),LDelta,qlogis(eps))
  lambda = exp(par[1])
  eps=plogis(par[3])
  tmp <- clean_production_data %>%
    mutate(Pred=if_else(World=="a",RSAspeaker_a(Prior,lambda,par[2]),RSAspeaker_ab(Prior,lambda,par[2])),
           loglik=log(case_when(
             Message=="a" ~ (Pred+eps)/(1+3*eps),
             World=="a"&Message%in%c("a&!b","!b") ~ (1-Pred+eps)/(1+3*eps),
             World=="ab"&Message=="a&b" ~ (1-Pred+eps)/(1+3*eps),
             T ~ eps/(1+3*eps)
           ))
    )
  return(-sum(tmp$loglik))
}

par_free <- optim(c(1,1,1,qlogis(1e-3)),neg_log_lik)
par_pos <- optim(c(1,1,1,qlogis(1e-3)),neg_log_lik_PQ)
par_equal <- optim(c(1,1,qlogis(1e-3)),neg_log_lik_EQ)

# Parameters:
round(c(exp(par_free$par[1]),par_free$par[2]/exp(par_free$par[1]),par_free$par[3]/exp(par_free$par[1]),plogis(par_free$par[4])),2)
round(c(exp(par_pos$par[1]),exp(par_pos$par[2]-par_pos$par[1]),exp(par_pos$par[3]-par_pos$par[1]),plogis(par_pos$par[4])),2)
round(c(exp(par_equal$par[1]),par_equal$par[2]/exp(par_equal$par[1]),plogis(par_equal$par[3])),2)

# Log-likelihoods:
-par_free$value
-par_pos$value
-par_equal$value


# Plot the resulting fit (choose which set of parameters you're interested in).

# Positive costs:
FittedSpeaker <- expand.grid(Prior=seq(0,1,by=.002), World=c("a","ab"))
FittedSpeaker$Prediction <-  ifelse(FittedSpeaker$World=="a",(RSAspeaker_a(FittedSpeaker$Prior,exp(par_pos$par[1]),exp(par_pos$par[3]))+plogis(par_pos$par[4]))/(1+3*plogis(par_pos$par[4])),(1-RSAspeaker_ab(FittedSpeaker$Prior,exp(par_pos$par[1]),exp(par_pos$par[2]))+plogis(par_pos$par[4]))/(1+3*plogis(par_pos$par[4])))
FittedSpeaker$World <- factor(FittedSpeaker$World,levels=c("a","ab"))
levels(FittedSpeaker$World) <- c(expression(S*"(A|"*w[a]*")"),expression(S*"(A&B|"*w[ab]*")"))

# Equal costs:
FittedSpeaker <- expand.grid(Prior=seq(0,1,by=.002), World=c("a","ab"))
FittedSpeaker$Prediction <-  ifelse(FittedSpeaker$World=="a",(RSAspeaker_a(FittedSpeaker$Prior,exp(par_equal$par[1]),par_equal$par[2])+plogis(par_equal$par[3]))/(1+3*plogis(par_equal$par[3])),(1-RSAspeaker_ab(FittedSpeaker$Prior,exp(par_equal$par[1]),par_equal$par[2])+plogis(par_equal$par[3]))/(1+3*plogis(par_equal$par[3])))
FittedSpeaker$World <- factor(FittedSpeaker$World,levels=c("a","ab"))
levels(FittedSpeaker$World) <- c(expression(S*"(A|"*w[a]*")"),expression(S*"(A&B|"*w[ab]*")"))

# Free costs:
FittedSpeaker <- expand.grid(Prior=seq(0,1,by=.002), World=c("a","ab"))
FittedSpeaker$Prediction <-  ifelse(FittedSpeaker$World=="a",(RSAspeaker_a(FittedSpeaker$Prior,exp(par_free$par[1]),par_free$par[3])+plogis(par_free$par[4]))/(1+3*plogis(par_free$par[4])),(1-RSAspeaker_ab(FittedSpeaker$Prior,exp(par_free$par[1]),par_free$par[2])+plogis(par_free$par[4]))/(1+3*plogis(par_free$par[4])))
FittedSpeaker$World <- factor(FittedSpeaker$World,levels=c("a","ab"))
levels(FittedSpeaker$World) <- c(expression(S*"(A|"*w[a]*")"),expression(S*"(A&B|"*w[ab]*")"))

pdf(file="RSA_speaker_Empirical_listener.pdf",width=7,height=5)
ggplot(data=production_plot_data,aes(x=x/100,y=freq,group=World,fill=n)) +
  facet_grid(.~World,labeller = label_parsed)+
  geom_point(data=filter(production_plot_data,x==100),
             aes(x=x/100,y=freq,fill=n,group=Message),
             inherit.aes = F,size=8,color="black",pch=21) +
  geom_bar(data=production_plot_data,
           aes(x=x/100,y=freq,fill=n,group=Message,width=width/100),
           inherit.aes = F,stat="identity",color="transparent") +
  geom_line(data=FittedSpeaker,aes(x=Prior,y=Prediction,color=World),
            inherit.aes = F,size=1,lineend="round",show.legend = FALSE)+
  scale_fill_gradient(low=rgb(0,0,0,.25),high=rgb(0,0,0,1),trans="custom_log",breaks=c(2,5,15,25))+
  scale_x_continuous(name=expression("Prior"~P(w[ab])),labels=scales::percent)+
  scale_y_continuous(name="Frequency",labels=scales::percent)+
  scale_color_manual(values=c(rgb(.7,0,.7),rgb(1,.7,0)))+
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=16),
    axis.text = element_text(size=12),
    strip.text = element_text(size=16)
  )
dev.off()

# How much A can be used to convey w_ab when we impose positive costs?
RSAspeaker_ab(1,exp(par_pos$par[1]),exp(par_pos$par[2]))

#######################################
# Define initial values for each model
#######################################


init_RSA <- list(
  list(lambda=2,deltaAB=1,deltaAnB=1,sigma_a=.3,sigma_ab=.1),
  list(lambda=10,deltaAB=1,deltaAnB=.1,sigma_a=.3,sigma_ab=.1),
  list(lambda=1,deltaAB=1,deltaAnB=1,sigma_a=.3,sigma_ab=.1),
  list(lambda=5,deltaAB=.1,deltaAnB=1,sigma_a=.3,sigma_ab=.1),
  list(lambda=5,deltaAB=1,deltaAnB=1,sigma_a=.3,sigma_ab=.2),
  list(lambda=10,deltaAB=.5,deltaAnB=.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=3,deltaAB=1.5,deltaAnB=1.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=.2,deltaAB=.5,deltaAnB=3,sigma_a=.3,sigma_ab=.2),
  list(lambda=.2,deltaAB=5,deltaAnB=5,sigma_a=.3,sigma_ab=.2),
  list(lambda=.2,deltaAB=5,deltaAnB=.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=2,deltaAB=5,deltaAnB=.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=2,deltaAB=5,deltaAnB=5,sigma_a=.3,sigma_ab=.2),
  list(lambda=2,deltaAB=.1,deltaAnB=.1,sigma_a=.3,sigma_ab=.2)
)
init_wRSA <- list(
  list(lambda=5,deltaAB=1,deltaAnB=1,omega=.95,sigma_a=.5,sigma_ab=.2),
  list(lambda=5,deltaAB=1,deltaAnB=1,omega=.5,sigma_a=.5,sigma_ab=.2),
  list(lambda=5,deltaAB=1,deltaAnB=1,omega=.25,sigma_a=.5,sigma_ab=.2),
  list(lambda=.2,deltaAB=1,deltaAnB=1,omega=.25,sigma_a=.5,sigma_ab=.2),
  list(lambda=15,deltaAB=1,deltaAnB=1,omega=.25,sigma_a=.5,sigma_ab=.2),
  list(lambda=6,deltaAB=1,deltaAnB=1,omega=.5,sigma_a=.5,sigma_ab=.2),
  list(lambda=4,deltaAB=1,deltaAnB=1,omega=.5,sigma_a=.5,sigma_ab=.2),
  list(lambda=10,deltaAB=.1,deltaAnB=1,omega=.97,sigma_a=.2,sigma_ab=.2),
  list(lambda=15,deltaAB=.1,deltaAnB=1,omega=.97,sigma_a=.2,sigma_ab=.2),
  list(lambda=10,deltaAB=1,deltaAnB=1,omega=.98,sigma_a=.2,sigma_ab=.2),
  list(lambda=10,deltaAB=.1,deltaAnB=1,omega=.99,sigma_a=.2,sigma_ab=.2),
  list(lambda=1,deltaAB=.1,deltaAnB=1,omega=.01,sigma_a=.2,sigma_ab=.2),
  list(lambda=.2,deltaAB=.1,deltaAnB=1,omega=.01,sigma_a=.2,sigma_ab=.2),
  list(lambda=4,deltaAB=.1,deltaAnB=1,omega=.01,sigma_a=.2,sigma_ab=.2),
  list(lambda=.2,deltaAB=1,deltaAnB=5,omega=.5,sigma_a=.2,sigma_ab=.2),
  list(lambda=.2,deltaAB=1,deltaAnB=5,omega=.85,sigma_a=.2,sigma_ab=.2)
)
init_svRSA <- list(
  list(lambda=5,deltaAB=1,deltaAnB=1,q=.95,sigma_a=.3,sigma_ab=.2),
  list(lambda=5,deltaAB=1,deltaAnB=1,q=.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=5,deltaAB=1,deltaAnB=1,q=.25,sigma_a=.3,sigma_ab=.2),
  list(lambda=1,deltaAB=1,deltaAnB=1,q=.25,sigma_a=.3,sigma_ab=.2),
  list(lambda=15,deltaAB=1,deltaAnB=1,q=.25,sigma_a=.3,sigma_ab=.2),
  list(lambda=6,deltaAB=1,deltaAnB=1,q=.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=5,deltaAB=1,deltaAnB=1,q=.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=10,deltaAB=.1,deltaAnB=1,q=.97,sigma_a=.3,sigma_ab=.2),
  list(lambda=25,deltaAB=10,deltaAnB=10,q=.95,sigma_a=.3,sigma_ab=.2),
  list(lambda=10,deltaAB=.1,deltaAnB=1,q=.99,sigma_a=.3,sigma_ab=.2),
  list(lambda=4,deltaAB=3,deltaAnB=.09,q=.89,sigma_a=.3,sigma_ab=.2),
  list(lambda=5,deltaAB=3,deltaAnB=.07,q=.89,sigma_a=.3,sigma_ab=.2),
  list(lambda=3,deltaAB=3,deltaAnB=.05,q=.1,sigma_a=.3,sigma_ab=.2),
  list(lambda=5,deltaAB=3,deltaAnB=.01,q=.25,sigma_a=.3,sigma_ab=.2),
  list(lambda=.5,deltaAB=3,deltaAnB=3,q=.75,sigma_a=.3,sigma_ab=.2),
  list(lambda=.2,deltaAB=3,deltaAnB=3,q=.75,sigma_a=.3,sigma_ab=.2),
  list(lambda=.2,deltaAB=1,deltaAnB=1,q=.75,sigma_a=.3,sigma_ab=.2)
)
init_FREE_LU <- list(
  list(lambda=2,deltaAB=1,deltaAnB=1,sigma_a=.3,sigma_ab=.1),
  list(lambda=10,deltaAB=1,deltaAnB=.1,sigma_a=.3,sigma_ab=.1),
  list(lambda=1,deltaAB=1,deltaAnB=1,sigma_a=.3,sigma_ab=.1),
  list(lambda=5,deltaAB=.1,deltaAnB=1,sigma_a=.3,sigma_ab=.1),
  list(lambda=5,deltaAB=1,deltaAnB=1,sigma_a=.3,sigma_ab=.2),
  list(lambda=10,deltaAB=.5,deltaAnB=.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=3,deltaAB=1.5,deltaAnB=1.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=.2,deltaAB=.5,deltaAnB=3,sigma_a=.3,sigma_ab=.2),
  list(lambda=.2,deltaAB=5,deltaAnB=5,sigma_a=.3,sigma_ab=.2),
  list(lambda=.2,deltaAB=5,deltaAnB=.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=2,deltaAB=5,deltaAnB=.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=2,deltaAB=5,deltaAnB=5,sigma_a=.3,sigma_ab=.2),
  list(lambda=2,deltaAB=.1,deltaAnB=.1,sigma_a=.3,sigma_ab=.2),
  list(lambda=100,deltaAB=.01,deltaAnB=.01,sigma_a=.4,sigma_ab=.2),
  list(lambda=1e5,deltaAB=0,deltaAnB=0,sigma_a=.4,sigma_ab=.2)
)
init_RSA_LI <- list(
  list(lambda=2,deltaAB=1,deltaAnB=1,sigma_a=.3,sigma_ab=.1),
  list(lambda=10,deltaAB=1,deltaAnB=.1,sigma_a=.3,sigma_ab=.1),
  list(lambda=1,deltaAB=1,deltaAnB=1,sigma_a=.3,sigma_ab=.1),
  list(lambda=5,deltaAB=.1,deltaAnB=1,sigma_a=.3,sigma_ab=.1),
  list(lambda=5,deltaAB=1,deltaAnB=1,sigma_a=.3,sigma_ab=.2),
  list(lambda=10,deltaAB=.5,deltaAnB=.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=3,deltaAB=1.5,deltaAnB=1.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=.2,deltaAB=.5,deltaAnB=3,sigma_a=.3,sigma_ab=.2),
  list(lambda=.2,deltaAB=5,deltaAnB=5,sigma_a=.3,sigma_ab=.2),
  list(lambda=.2,deltaAB=5,deltaAnB=.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=2,deltaAB=5,deltaAnB=.5,sigma_a=.3,sigma_ab=.2),
  list(lambda=2,deltaAB=5,deltaAnB=5,sigma_a=.3,sigma_ab=.2),
  list(lambda=2,deltaAB=.1,deltaAnB=.1,sigma_a=.3,sigma_ab=.2),
  list(lambda=100,deltaAB=.01,deltaAnB=.01,sigma_a=.4,sigma_ab=.2),
  list(lambda=1e5,deltaAB=0,deltaAnB=0,sigma_a=.4,sigma_ab=.2)
)


#################
# Fit the models
#################


# Takes up to a few minutes:
t <- Sys.time()
parms_RSA <- complete_model_fit(RSA_Production_LogLikelihood,
                                RSA_Comprehension_LogLikelihood,
                                init_RSA,
                                clean_production_data,clean_comprehension_data,
                                model_name = "RSA")
parms_wRSA <- complete_model_fit(wRSA_Production_LogLikelihood,
                                 wRSA_Comprehension_LogLikelihood,
                                init_wRSA,
                                clean_production_data,clean_comprehension_data,
                                model_name = "wRSA")
parms_BwRSA <- complete_model_fit(BwRSA_Production_LogLikelihood,
                                  BwRSA_Comprehension_LogLikelihood,
                                 init_wRSA,
                                 clean_production_data,clean_comprehension_data,
                                 model_name = "BwRSA")
# svRSA with P(Qa)S2(|Qa)+P(Qb)S2(|Qb)
parms_svRSA1 <- complete_model_fit(svRSA1_Production_LogLikelihood,
                                     svRSA1_Comprehension_LogLikelihood,
                                 init_svRSA,
                                 clean_production_data,clean_comprehension_data,
                                 model_name = "svRSA1")
# svRSA with S2(|Qb)
parms_svRSA2 <- complete_model_fit(svRSA2_Production_LogLikelihood,
                                     svRSA2_Comprehension_LogLikelihood,
                                     init_svRSA,
                                     clean_production_data,clean_comprehension_data,
                                     model_name = "svRSA2")
# RSA LI with S1
parms_RSA_LI1 <- complete_model_fit(RSA_LI1_Production_LogLikelihood,
                                    RSA_LI1_Comprehension_LogLikelihood,
                                    init_RSA_LI,
                                    clean_production_data,clean_comprehension_data,
                                    model_name = "RSA_LI1")
# RSA LI with S2
parms_RSA_LI2 <- complete_model_fit(RSA_LI2_Production_LogLikelihood,
                                RSA_LI2_Comprehension_LogLikelihood,
                                init_RSA_LI,
                                clean_production_data,clean_comprehension_data,
                                model_name = "RSA_LI2")
parms_FREE_LU <- complete_model_fit(FREE_LU_Production_LogLikelihood,
                                    FREE_LU_Comprehension_LogLikelihood,
                                   init_FREE_LU,
                                   clean_production_data,clean_comprehension_data,
                                   model_name = "FREE_LU")
parms_EXH_LU <- complete_model_fit(EXH_LU_Production_LogLikelihood,
                                     EXH_LU_Comprehension_LogLikelihood,
                                    init_FREE_LU,
                                    clean_production_data,clean_comprehension_data,
                                    model_name = "EXH_LU")
print(Sys.time()-t)


summary_table <- bind_rows(parms_RSA,parms_wRSA,parms_BwRSA,parms_svRSA1,parms_svRSA2,parms_RSA_LI1,parms_RSA_LI2,parms_FREE_LU,parms_EXH_LU) %>%
  mutate(xi = pmax(omega,q,na.rm=T),
         xi_prod = pmax(omega_prod,q_prod,na.rm=T)
  ) %>%
  select(-matches("pEXH|omega|q.*")) %>%
  arrange(AIC) %>%
  select(-logLik,-dropped,-AIC,-BIC,logLik,dropped,AIC,BIC)

# Keeping only joint submodel or picking best one doesn't change the order much (and crucially not for the best models)
summary_table %>%
  filter(submodel=="joint") %>%
  select(model,submodel,logLik,AIC,BIC)

# Only difference is the position of BwRSA:
cbind(summary_table %>% filter(submodel=="joint") %>% pull(model),
      summary_table %>% group_by(model) %>% filter(AIC==min(AIC)) %>% ungroup() %>% pull(model))

###########################################################
# Create a Latex table and format the parameters for plots
###########################################################

# Keep only the joint parameters version of each model and print a decent latex table:
optimal_submodels <- summary_table %>%
  filter(submodel=="joint") %>%
  select(-ends_with("prod"),-dropped,-submodel)

optimal_submodels %>%
  mutate(deltaAB=round(deltaAB,3),
         deltaAnB=round(deltaAnB,3),
  ) %>%
  select(model,lambda,deltaAB,deltaAnB,xi,sigma_a,sigma_ab,eps,logLik,AIC,BIC) %>%
  rename_at(vars(names(.)),~paste("\\multicolumn{1}{c|}{",
                                  c("Model","$\\lambda$","$\\Delta_{ab}$","$\\Delta_{a\\neg b}$","$\\xi$","$\\sigma_a$","$\\sigma_{ab}$","$\\epsilon$","$\\ell$","AIC","BIC"),
                                  "}",sep="")
  ) %>%
  xtable(.,align="|l|c|c|c|c|c|c|c|l|l|c|l|",
         display=c("s","s","g","fg","fg","fg","f","f","f","f","f","f"),
         digits=c(1,1,2,2,2,2,2,2,3,0,0,0),
         caption="Joint parameter fit (production+comprehension) for each theoretical model") %>%
  print(.,sanitize.text.function=function(x){x},latex.environments = "center",include.rownames=FALSE)

# Save parameters for plots:
optimal_parameters <- summary_table %>%
  filter(submodel=="joint") %>%
  mutate(lambda_prod=if_else(is.na(lambda_prod),lambda,lambda_prod),
         deltaAB_prod=if_else(is.na(deltaAB_prod),deltaAB,deltaAB_prod),
         deltaAnB_prod=if_else(is.na(deltaAnB_prod),deltaAnB,deltaAnB_prod),
         xi_prod=if_else(is.na(xi_prod),xi,xi_prod),
         sigma_ab=if_else(is.na(sigma_ab),sigma_a,sigma_ab)
  ) %>%
  column_to_rownames("model") %>%
  select(-c(logLik,dropped,AIC,submodel,BIC))


#####################################################
# Plot the models fit against the comprehension data
#####################################################

ModelListenerPredictions <- expand_grid(Model=row.names(optimal_parameters),
                                             Utterance=c("a","a&b"),
                                             Prior=.005+.99*seq(0,1,.001)) %>%
  group_by(Model) %>%
  mutate(Prediction=if_else(Utterance=="a",
                            do.call(paste0(first(Model),"_L1_a"),unname(c(list(Prior),as.list(na.omit(as_vector(optimal_parameters[first(Model),c(1:3,10)])))))),
                            do.call(paste0(first(Model),"_L1_ab"),unname(c(list(Prior),as.list(na.omit(as_vector(optimal_parameters[first(Model),c(1:3,10)]))))))
                            ),
         Sigma = case_when(
           Utterance=="a" ~ optimal_parameters[first(Model),"sigma_a"],
           Utterance=="a&b" ~ optimal_parameters[first(Model),"sigma_ab"],
         )) %>%
  ungroup() %>%
  mutate(FirstQuart=case_when(
           qnorm(0.25,Prediction,Sigma) <= 0 ~ 0,
           qnorm(0.25,Prediction,Sigma) >= 1 ~ 1,
           T ~ qnorm(0.25,Prediction,Sigma)
         ),
         ThirdQuart=case_when(
           qnorm(0.75,Prediction,Sigma) <= 0 ~ 0,
           qnorm(0.75,Prediction,Sigma) >= 1 ~ 1,
           T ~ qnorm(0.75,Prediction,Sigma)
         ),
         Model=factor(Model,levels=c("RSA","wRSA","BwRSA","FREE_LU","EXH_LU","RSA_LI1","RSA_LI2","svRSA1","svRSA2"))
         )

# Plot for both utterances:
pdf(file="Comprehension_FreeCosts.pdf",width=10,height=6)
ModelListenerPredictions %>%
  ggplot(aes(x=Prior,y=Prediction,col=Utterance,ymin=FirstQuart,ymax=ThirdQuart,fill=Utterance)) +
  facet_wrap(.~Model)+
  geom_abline(slope = 1,intercept = 0,linetype=2)+
  geom_line()+
  geom_ribbon(alpha=.1)+
  geom_point(data=clean_comprehension_data,aes(x=Prior,y=Posterior,color=Utterance),shape=1,size=.2,inherit.aes = F)+
  xlab(expression("Prior"~P(w[ab])))+
  ylab(expression("Posterior"~P(w[ab]*"|"*"u")))+
  scale_x_continuous(limits=c(0,1),labels = scales::percent)+
  scale_y_continuous(limits=c(0,1),labels = scales::percent)+
  scale_color_manual(values = c(rgb(.7,0,.7),rgb(1,.7,0)),labels=c("A","A&B"),aesthetics = c("colour", "fill"))+
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(10, "mm")
  )
dev.off()

##################################################
# Plot the models fit against the production data
##################################################

ModelSpeakerPredictions <- expand_grid(Model=row.names(optimal_parameters),
                                             World=c("a","ab"),
                                             Prior=.005+.99*seq(0,1,.001)) %>%
  group_by(Model) %>%
  mutate(
    Prediction = if_else(World=="a",
                         do.call(paste0(first(Model),"_S2_a_wa"),unname(c(list(Prior),as.list(na.omit(as_vector(optimal_parameters[first(Model),c(7:9,11)])))))),
                         do.call(paste0(first(Model),"_S2_ab_wab"),unname(c(list(Prior),as.list(na.omit(as_vector(optimal_parameters[first(Model),c(7:9,11)]))))))
                         ),
    Prediction = (exp(Prediction)+optimal_parameters[first(Model),"eps"])/(1+4*optimal_parameters[first(Model),"eps"])
  ) %>%
  ungroup() %>%
  mutate(World = factor(World,levels=c("a","ab")),
         Model=factor(Model,levels=c("RSA","wRSA","BwRSA","FREE_LU","EXH_LU","RSA_LI1","RSA_LI2","svRSA1","svRSA2")))

production_plot_data <- clean_production_data %>%
  mutate(Prior=(Prior-0.005)/.99,
         bin=(Prior%/%.25),
         Message=if_else(Message=="!b","a&!b",as.character(Message))) %>%
  group_by(World,bin,Message) %>%
  summarise(n=n()) %>%
  mutate(freq=n/sum(n),
         x=if_else(bin==4,1,.25*bin+.12),
         width=if_else(bin==4,.01,.25),
         alpha=log(n)+1,
         alpha=alpha/(log(19)+1)) %>%
  filter((Message=="a"&World=="a")|(Message=="a&b"&World=="ab"))
production_plot_data$World <- factor(production_plot_data$World,levels=c("a","ab"))

speaker_labels <- c(expression(S[2]*"(A|"*w[a]*")"),expression(S[2]*"(A&B|"*w[ab]*")"))

pdf(file="Production_FreeCosts.pdf",width=10,height=6)
ModelSpeakerPredictions %>%
  ggplot(aes(x=Prior,y=Prediction,col=World)) +
  facet_wrap(.~Model)+
  geom_point(data=filter(production_plot_data,x==1&World=="a"),
             aes(x=x,y=freq,fill=n,group=Message),
             inherit.aes = F,size=8,color="black",pch=21) +
  geom_bar(data=filter(production_plot_data,World=="a"),
           aes(x=x,y=freq,fill=n,group=Message),
           width=filter(production_plot_data,World=="a")$width,
           inherit.aes = F,stat="identity",color="transparent") +
  geom_line(aes(x=Prior,y=Prediction,col=World),size=1,lineend="round")+
  scale_fill_gradient(low=rgb(0,0,0,.25),high=rgb(0,0,0,1),trans="custom_log",breaks=c(10,15,20,25))+
  scale_y_continuous(limits=c(0,1),expand = c(0, 0),labels = scales::percent,name=expression("Speaker probability"))+
  scale_x_continuous(breaks=seq(0,1,.25),minor_breaks=NULL,labels = scales::percent,name=expression("Prior"~P(w[ab])))+
  scale_color_manual(values=c(rgb(.7,0,.7),rgb(1,.7,0)),labels = speaker_labels,name="")+
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(),
    legend.text.align = 0,
    panel.spacing.x = unit(10, "mm")
  )
dev.off()

################################
# Posthoc analysis: Equal costs
################################

# Function to make the likelihood functions use deltaAB instead of deltaAnB:
equal_costs <- function(FUN){
  new_FUN <- FUN
  string_body <- gsub("deltaAnB","deltaAB",as.character(body(FUN)))[-1]
  body(new_FUN) <- as.call(c(as.name("{"),Vectorize(str2lang)(string_body)))
  return(new_FUN)
  }

# Takes about 1min in all:
t <- Sys.time()
parms_RSA_ec <- complete_model_fit(equal_costs(RSA_Production_LogLikelihood),
                                equal_costs(RSA_Comprehension_LogLikelihood),
                                init_RSA,
                                clean_production_data,clean_comprehension_data,
                                model_name = "RSA")
parms_wRSA_ec <- complete_model_fit(equal_costs(wRSA_Production_LogLikelihood),
                                 equal_costs(wRSA_Comprehension_LogLikelihood),
                                 init_wRSA,
                                 clean_production_data,clean_comprehension_data,
                                 model_name = "wRSA")
parms_BwRSA_ec <- complete_model_fit(equal_costs(BwRSA_Production_LogLikelihood),
                                  equal_costs(BwRSA_Comprehension_LogLikelihood),
                                  init_wRSA,
                                  clean_production_data,clean_comprehension_data,
                                  model_name = "BwRSA")
parms_svRSA1_ec <- complete_model_fit(equal_costs(svRSA1_Production_LogLikelihood),
                                        equal_costs(svRSA1_Comprehension_LogLikelihood),
                                        init_svRSA,
                                        clean_production_data,clean_comprehension_data,
                                        model_name = "svRSA1")
parms_svRSA2_ec <- complete_model_fit(equal_costs(svRSA2_Production_LogLikelihood),
                                        equal_costs(svRSA2_Comprehension_LogLikelihood),
                                        init_svRSA,
                                        clean_production_data,clean_comprehension_data,
                                        model_name = "svRSA2")
parms_RSA_LI1_ec <- complete_model_fit(equal_costs(RSA_LI1_Production_LogLikelihood),
                                       equal_costs(RSA_LI1_Comprehension_LogLikelihood),
                                       init_RSA_LI,
                                       clean_production_data,clean_comprehension_data,
                                       model_name = "RSA_LI1")
parms_RSA_LI2_ec <- complete_model_fit(equal_costs(RSA_LI2_Production_LogLikelihood),
                                       equal_costs(RSA_LI2_Comprehension_LogLikelihood),
                                       init_RSA_LI,
                                       clean_production_data,clean_comprehension_data,
                                       model_name = "RSA_LI2")
parms_FREE_LU_ec <- complete_model_fit(equal_costs(FREE_LU_Production_LogLikelihood),
                                      equal_costs(FREE_LU_Comprehension_LogLikelihood),
                                      init_FREE_LU,
                                      clean_production_data,clean_comprehension_data,
                                      model_name = "FREE_LU")
parms_EXH_LU_ec <- complete_model_fit(equal_costs(EXH_LU_Production_LogLikelihood),
                                       equal_costs(EXH_LU_Comprehension_LogLikelihood),
                                       init_FREE_LU,
                                       clean_production_data,clean_comprehension_data,
                                       model_name = "EXH_LU")
print(Sys.time()-t)



# We need to recompute the AIC/BIC to avoid counting the dummy deltaAnB in the degrees of freedom:
N = nrow(clean_comprehension_data)+nrow(clean_production_data)
summary_table_ec <- bind_rows(parms_RSA_ec,parms_wRSA_ec,parms_BwRSA_ec,
                              parms_svRSA1_ec,parms_svRSA2_ec,parms_RSA_LI1_ec,
                              parms_RSA_LI2_ec,parms_FREE_LU_ec,parms_EXH_LU_ec) %>%
  mutate(xi = pmax(omega,q,na.rm=T),
         xi_prod = pmax(omega_prod,q_prod,na.rm=T)
  ) %>%
  select(-matches("pEXH|omega|q.*")) %>%
  select(-logLik,-dropped,-AIC,-BIC,logLik,dropped,AIC,BIC) %>%
  mutate(AIC = 2*dropped+2*rowSums(!is.na(.[c(3:4,6:10,12:13)]))-2*logLik,
         BIC = log(N-dropped)*(dropped+rowSums(!is.na(.[c(3:4,6:10,12:13)])))-2*logLik) %>%
  arrange(AIC)


summary_table_ec %>%
  filter(submodel=="joint") %>%
  select(model,logLik,AIC,BIC)

# No difference at all this time:
cbind(summary_table_ec %>% filter(submodel=="joint") %>% pull(model),
      summary_table_ec %>% group_by(model) %>% filter(AIC==min(AIC)) %>% ungroup() %>% pull(model))



###########################################################
# Create a Latex table and format the parameters for plots
###########################################################

# Keep only the joint parameters version of each model and print a decent latex table:
optimal_submodels_ec <- summary_table_ec %>%
  filter(submodel=="joint") %>%
  select(-deltaAnB,-deltaAnB_prod)

optimal_submodels_ec %>%
  select(model,lambda,deltaAB,xi,sigma_a,sigma_ab,eps,logLik,AIC,BIC) %>%
  rename_at(vars(names(.)),~paste("\\multicolumn{1}{c|}{",
                                  c("Model","$\\lambda$","$\\Delta_{ab}$","$\\xi$","$\\sigma_a$","$\\sigma_{ab}$","$\\epsilon$","$\\ell$","AIC","BIC"),
                                  "}",sep="")
  ) %>%
  xtable(.,align="|l|c|c|c|c|c|c|l|l|c|l|",
         display=c("s","s","fg","fg","fg","f","f","f","f","f","f"),
         digits=c(1,1,2,2,2,2,2,3,0,0,0),
         caption="Best AIC submodel for each theoretical model") %>%
  print(.,sanitize.text.function=function(x){x},latex.environments = "center",include.rownames=FALSE)

# Save parameters for plots:
optimal_parameters_ec <- summary_table_ec %>%
  filter(submodel=="joint") %>%
  mutate(lambda_prod=if_else(is.na(lambda_prod),lambda,lambda_prod),
         deltaAB_prod=if_else(is.na(deltaAB_prod),deltaAB,deltaAB_prod),
         deltaAnB=deltaAB,
         deltaAnB_prod=deltaAB_prod,
         xi_prod=if_else(is.na(xi_prod),xi,xi_prod),
         sigma_ab=if_else(is.na(sigma_ab),sigma_a,sigma_ab)
  ) %>%
  column_to_rownames("model") %>%
  select(-c(logLik,dropped,AIC,submodel,BIC))


#####################################################
# Plot the models fit against the comprehension data
#####################################################

ModelListenerPredictionsEC <- expand_grid(Model=row.names(optimal_parameters_ec),
                                        Utterance=c("a","a&b"),
                                        Prior=.005+.99*seq(0,1,.001)) %>%
  group_by(Model) %>%
  mutate(Prediction=if_else(Utterance=="a",
                            do.call(paste0(first(Model),"_L1_a"),unname(c(list(Prior),as.list(na.omit(as_vector(optimal_parameters_ec[first(Model),c(1:3,10)])))))),
                            do.call(paste0(first(Model),"_L1_ab"),unname(c(list(Prior),as.list(na.omit(as_vector(optimal_parameters_ec[first(Model),c(1:3,10)]))))))
  ),
  Sigma = case_when(
    Utterance=="a" ~ optimal_parameters_ec[first(Model),"sigma_a"],
    Utterance=="a&b" ~ optimal_parameters_ec[first(Model),"sigma_ab"],
  )) %>%
  ungroup() %>%
  mutate(FirstQuart=case_when(
    qnorm(0.25,Prediction,Sigma) <= 0 ~ 0,
    qnorm(0.25,Prediction,Sigma) >= 1 ~ 1,
    T ~ qnorm(0.25,Prediction,Sigma)
  ),
  ThirdQuart=case_when(
    qnorm(0.75,Prediction,Sigma) <= 0 ~ 0,
    qnorm(0.75,Prediction,Sigma) >= 1 ~ 1,
    T ~ qnorm(0.75,Prediction,Sigma)
  ),
  Model=factor(Model,levels=c("RSA","wRSA","BwRSA","FREE_LU","EXH_LU","RSA_LI1","RSA_LI2","svRSA1","svRSA2"))
  )

# Plot for both utterances:
pdf(file="Comprehension_EqualCosts.pdf",width=10,height=6)
ModelListenerPredictionsEC %>%
  ggplot(aes(x=Prior,y=Prediction,col=Utterance,ymin=FirstQuart,ymax=ThirdQuart,fill=Utterance)) +
  facet_wrap(.~Model)+
  geom_abline(slope = 1,intercept = 0,linetype=2)+
  geom_line()+
  geom_ribbon(alpha=.1)+
  geom_point(data=clean_comprehension_data,aes(x=Prior,y=Posterior,color=Utterance),shape=1,size=.2,inherit.aes = F)+
  xlab(expression("Prior"~P(w[ab])))+
  ylab(expression("Posterior"~P(w[ab]*"|"*"u")))+
  scale_x_continuous(limits=c(0,1),labels = scales::percent)+
  scale_y_continuous(limits=c(0,1),labels = scales::percent)+
  scale_color_manual(values = c(rgb(.7,0,.7),rgb(1,.7,0)),labels=c("A","A&B"),aesthetics = c("colour", "fill"))+
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(10, "mm")
  )
dev.off()

##################################################
# Plot the models fit against the production data
##################################################

ModelSpeakerPredictionsEC <- expand_grid(Model=row.names(optimal_parameters_ec),
                                       World=c("a","ab"),
                                       Prior=.005+.99*seq(0,1,.001)) %>%
  group_by(Model) %>%
  mutate(
    Prediction = if_else(World=="a",
                         do.call(paste0(first(Model),"_S2_a_wa"),unname(c(list(Prior),as.list(na.omit(as_vector(optimal_parameters_ec[first(Model),c(7:9,11)])))))),
                         do.call(paste0(first(Model),"_S2_ab_wab"),unname(c(list(Prior),as.list(na.omit(as_vector(optimal_parameters_ec[first(Model),c(7:9,11)]))))))
    ),
    Prediction = (exp(Prediction)+optimal_parameters_ec[first(Model),"eps"])/(1+4*optimal_parameters_ec[first(Model),"eps"])
  ) %>%
  ungroup() %>%
  mutate(Prior=Prior,
         World = factor(World,levels=c("a","ab")),
         Model=factor(Model,levels=c("RSA","wRSA","BwRSA","FREE_LU","EXH_LU","RSA_LI1","RSA_LI2","svRSA1","svRSA2"))
  )

speaker_labels <- c(expression(S[2]*"(A|"*w[a]*")"),expression(S[2]*"(A&B|"*w[ab]*")"))

pdf(file="Production_EqualCosts.pdf",width=10,height=6)
ModelSpeakerPredictionsEC %>%
  ggplot(aes(x=Prior,y=Prediction,col=World)) +
  facet_wrap(.~Model)+
  geom_point(data=filter(production_plot_data,x==1&World=="a"),
             aes(x=x,y=freq,fill=n,group=Message),
             inherit.aes = F,size=8,color="black",pch=21) +
  geom_bar(data=filter(production_plot_data,World=="a"),
           aes(x=x,y=freq,fill=n,group=Message),
           width=filter(production_plot_data,World=="a")$width,
           inherit.aes = F,stat="identity",color="transparent") +
  geom_line(aes(x=Prior,y=Prediction,col=World),size=1,lineend="round")+
  scale_fill_gradient(low=rgb(0,0,0,.25),high=rgb(0,0,0,1),trans="custom_log",breaks=c(10,15,20,25))+
  scale_y_continuous(limits=c(0,1),expand = c(0, 0),labels = scales::percent,name=expression("Speaker probability"))+
  scale_x_continuous(breaks=seq(0,1,.25),minor_breaks=NULL,labels = scales::percent,name=expression("Prior"~P(w[ab])))+
  scale_color_manual(values=c(rgb(.7,0,.7),rgb(1,.7,0)),labels = speaker_labels,name="")+
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(),
    legend.text.align = 0,
    panel.spacing.x = unit(10, "mm")
  )
dev.off()




# Note that the few a&b messages produces in w_a are not coming from lower prior values:
t.test(
  clean_production_data %>% filter(World=="a"&Message=="a&b") %>% pull(Prior),
  clean_production_data %>% filter(World=="a"&Message!="a&b") %>% pull(Prior)
  )


