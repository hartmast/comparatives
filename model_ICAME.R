# model ICAME


# load packages -----------------------------------------------------------

# install CRAN packages (if not yet installed)
sapply(c("lattice", "latticeExtra", "brms"), function(x) 
  if(!is.element(x, installed.packages())) install.packages(x, dependencies = T))

library(brms)
library(lattice)
library(latticeExtra)



# data wrangling ----------------------------------------------------------

# Read in data
data_icame <- read.csv("../../Data/Data for ICAME/data_icame.csv")
str(data_icame)

# Data clean-up
#---------------------------------------------------

# remove participle forms
data_icame_no_participles <- subset(data_icame, morphology_groups != "ed_participle" & morphology_groups != "ing_participle")
data_icame_no_participles <- droplevels(data_icame_no_participles)

# Number of syllables: Reduce to factor with three levels: 1, 2, 3+
data_icame_no_participles$n_syl_factor2 <- factor(ifelse(data_icame_no_participles$n_syll == 1, "1", ifelse(data_icame_no_participles$n_syll == 2, "2", "3+")))
str(data_icame_no_participles)

# Extrapolated counts: Convert to integers (for modeling)
data_icame_no_participles$est_er_grad_freq   <- as.integer(round(data_icame_no_participles$est_er_grad_freq, 0))
data_icame_no_participles$est_more_grad_freq <- as.integer(round(data_icame_no_participles$est_more_grad_freq, 0))

# Gradation frequency: Transform to log scale, then standardize to z-scores
data_icame_no_participles$log_rate_z <- as.numeric(scale(data_icame_no_participles$log_est_grad_pmw))

# Create a new variable for the phonological subgroups (number of syllables x final segment)
data_icame_no_participles$combinations <- factor(paste(as.character(data_icame_no_participles$n_syl_factor2), as.character(data_icame_no_participles$final_segments_class), sep = "_") )


# Data inspection
#----------------------------------------------------

# Overview of phonological subgroups
xtabs(~n_syl_factor2+final_segments_class, data_icame_no_participles)

# Gradation frequency
histogram(data_icame_no_participles$log_rate_z)

# -er
plot1 = histogram(log(data_icame$est_er_grad_freq), 
                  nint=50, col="grey", 
                  #par.settings=temp_my_settings, axis=axis_L, 
                  type="count",
          scales=list(x=list(at=log(c(1, 10, 100, 1000, 10000)), label=c(1, 10, 100, 1000, 10000)),
                      y=list(at=c(0, 50, 100, 150))), xlim=c(-.5, log(40000)), ylim=c(0, NA), 
          xlab="Frequency in the BNC", ylab="Number of adjectives")
# cairo_pdf(file = "./figures/frequency_er_BNC.pdf", height=3, width=3)
print(plot1)
# dev.off()

# more
plot1 = histogram(log(data_icame$est_more_grad_freq), nint=50, col="grey", 
                  #par.settings=temp_my_settings, 
                  # axis=axis_L, 
                  type="count",
          scales=list(x=list(at=log(c(1, 10, 100, 1000, 10000)), label=c(1, 10, 100, 1000, 10000)),
                      y=list(at=c(0, 150, 500, 1000))), xlim=c(-.5, log(40000)), ylim=c(0, 1150), 
          xlab="Frequency in the BNC", ylab="Number of adjectives")
#cairo_pdf(file = "./figures/frequency_more_BNC.pdf", height=3, width=3)
print(plot1)
#dev.off()



# Model
#----------------------------------------------------


# Logistic regression model with (too) weak priors
# m_new_2 = brm(est_more_grad_freq | trials(est_er_grad_freq + 
#                                             est_more_grad_freq) ~ 
#       # in brms, "success | trials(success + failure)" is
#       # equivalent to "cbind(success + failure)" in glm(er)
#                 log_rate_z*combinations + (1|type), 
#             data=data_icame_no_participles, 
#             family=binomial,
#             warmup=500, iter=2000, chains=4, cores=4,
#             prior=c(set_prior("normal(0,5)", class="b"),
#                     set_prior("normal(0,5)", class="sd")))
# 
# saveRDS(m_new_2, "./models/m_new_2.rds")
m_new_2 <- readRDS("../brms model/m_new_2.rds")


# null model
# m_null = brm(est_more_grad_freq | trials(est_er_grad_freq +
#                                             est_more_grad_freq) ~ 
#                (1|type),
#             data=data_icame_no_participles,
#             family=binomial,
#             warmup=500, iter=2000, chains=4, cores=4)
# 

# Compute posterior estimates
#----------------------------------------------------

# Extarct posterior samples from model
post <- posterior_samples(m_new_2)
str(post)

# Posterior estimates of interest: Gradation frequency trend 
# for each phonological subgroup
# Array map:
#     6 --- n = 5 quantiles from the posterior distribution (.025, .25, .5, .75, .975), i.e. the median + 50% and 95% uncertainty interval, everything on the proportion (of more) scale - the sixth column is filled with the steps along the gradation frequency range (see next point)
#   100 --- steps along the gradation frequency range to obtain a sufficiently dense sequence of estimates
#    17 --- number of phonological subgroups

estimates_array = array(NA, dim = c(6, 100, nlevels(data_icame_no_participles$combinations)))
str(estimates_array)

n_combinations = nlevels(data_icame_no_participles$combinations)
range(data_icame_no_participles$log_rate_z)

# fill up array with estimates
for(i in 1:n_combinations){
  x_range <- range(subset(data_icame_no_participles, combinations == levels(data_icame_no_participles$combinations)[i])$log_rate_z)
  x_seq <- seq(x_range[1], x_range[2], length.out=100)
  if(i == 1) {
    estimates_array[1:5,,i] <- plogis(apply(as.matrix(post[,c(1,2)]) %*% rbind(1, x_seq), 2, quantile, c(.025, .25, .5, .75, .975)))
    estimates_array[  6,,i] <- x_seq
  }
  else {
    estimates_array[1:5,,i] <- plogis(apply(as.matrix(post[,c(1,2, 1+i, n_combinations+i)]) %*% rbind(1, x_seq, 1, x_seq), 2, quantile, c(.025, .25, .5, .75, .975)))
    estimates_array[  6,,i] <- x_seq
  }
}



# Visualization
#----------------------------------------------------

# source("./r_utils/my_settings.R")
# temp_my_settings = my_settings
# temp_my_settings$layout.heights$bottom.padding = 1

# Overview: degree marking by gradation frequency

# add jitter
set.seed(1138)
data_icame_no_participles$est_prop_more_jitter <- jitter(data_icame_no_participles$est_prop_more, 100)
set.seed(1985)
data_icame_no_participles$log_rate_z_jitter <- jitter(data_icame_no_participles$log_rate_z, 400)


library(ggplot2)
library(plotly)
library(htmlwidgets)

p1 <- ggplot(data_icame_no_participles,
       aes(x = log_rate_z_jitter,
           y = est_prop_more_jitter,
           color = combinations,
           label = type)) + geom_text() +
  theme_bw() +
  scale_fill_manual(values = terrain.colors(17)) +
  scale_x_continuous(breaks = (log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw),
                     labels = c(0.01, 0.1, 1, 10, 100)) +
  ylab("Proportion of more") +
  xlab("Log gradation frequency")

p2 <- ggplotly(p1, tooltip = "type")

# saveWidget(as_widget(p2), "comparatives.html")

plot(data_icame_no_participles$log_rate_z_jitter,
     data_icame_no_participles$est_prop_more_jitter,
     type = "n")
text(data_icame_no_participles$log_rate_z_jitter,
     data_icame_no_participles$est_prop_more_jitter,
     data_icame_no_participles$type)

xyplot(
  jitter(est_prop_more, 100) ~ jitter(log_rate_z, 400), 
  data=data_icame_no_participles, xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.2, cex=.5, 
  # par.settings=temp_my_settings, axis=axis_L, 
  xlab="Log gradation frequency", ylab="Proportion of more",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100)),
              y=list(at=c(0, .5, 1))),
  panel=function(x,y,...){
    panel.abline(v=(log(6/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), col="grey")
    panel.xyplot(x,y,...)
  })


# Overview of trends grouped by number of syllables
plot_1_syllable <- 
  xyplot(
    jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="1_cluster"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05), type="n",
    # par.settings=temp_my_settings, axis=axis_L, 
    xlab=list(cex=.8, label="Log gradation frequency"), ylab=list(cex=.8, label="Proportion of more"),
    scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                       label=c(0.01, 0.1, 1, 10, 100)),
                y=list(at=c(0, .5, 1))),
    panel=function(x,y, subscripts,...){
      panel.points(x=estimates_array[6,,1], y=estimates_array[3,,1], type="l", col=1)
      panel.points(x=estimates_array[6,,2], y=estimates_array[3,,2], type="l", col="darkgrey")
      panel.points(x=estimates_array[6,,3], y=estimates_array[3,,3], type="l", col=1, lty=2)
      panel.points(x=estimates_array[6,,4], y=estimates_array[3,,4], type="l", col="darkgrey")
      panel.points(x=c(-.43, .4, .08, -.3), y=c(.89, .32, .35, .22), pch=16, col="white", cex=1)
      panel.text(x=-.45, y=.9, label="cluster", adj=0, cex=.8)
      panel.text(x=.4, y=.32, label="l", col="darkgrey", adj=0, cex=.8)
      panel.text(x=0.1, y=.35, label="other", adj=1, cex=.8)
      panel.text(x=-.3, y=.22, label="r", col="darkgrey", adj=0, cex=.8)
      })

plot_2_syllable <- 
  xyplot(
    jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="1_cluster"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05), type="n",
    # par.settings=temp_my_settings, axis=axis_L, 
    xlab=list(cex=.8, label=" "), ylab=list(cex=.8, label=" "),
    scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                       label=c(0.01, 0.1, 1, 10, 100)),
                y=list(at=c(0, .5, 1), lab=NULL)),
    panel=function(x,y, subscripts,...){
      panel.points(x=estimates_array[6,,5], y=estimates_array[3,,5], type="l", col=1)
      panel.points(x=estimates_array[6,,6], y=estimates_array[3,,6], type="l", col=1)
      panel.points(x=estimates_array[6,,7], y=estimates_array[3,,7], type="l", col=1)
      panel.points(x=estimates_array[6,,8], y=estimates_array[3,,8], type="l", col=1)
      panel.points(x=estimates_array[6,,9], y=estimates_array[3,,9], type="l", col=1)
      panel.points(x=estimates_array[6,,10], y=estimates_array[3,,10], type="l", col=1, lty=2)
      panel.text(x=c(0, .5), y=c(.25, .6), label=c("-y", "-ly"))
      })

plot_3_syllable <- 
  xyplot(
    jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="1_cluster"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05), type="n",
    # par.settings=temp_my_settings, axis=axis_L, 
    xlab=list(cex=.8, label=" "), ylab=list(cex=.8, label=" "),
    scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                       label=c(0.01, 0.1, 1, 10, 100)),
                y=list(at=c(0, .5, 1), lab=NULL)),
    panel=function(x,y, subscripts,...){
      panel.points(x=estimates_array[6,,11], y=estimates_array[3,,11], type="l", col=1)
      panel.points(x=estimates_array[6,,12], y=estimates_array[3,,12], type="l", col=1)
      panel.points(x=estimates_array[6,,13], y=estimates_array[3,,13], type="l", col=1)
      panel.points(x=estimates_array[6,,14], y=estimates_array[3,,14], type="l", col=1)
      panel.points(x=estimates_array[6,,15], y=estimates_array[3,,15], type="l", col=1)
      panel.points(x=estimates_array[6,,16], y=estimates_array[3,,16], type="l", col=1, lty=2)
      panel.points(x=estimates_array[6,,17], y=estimates_array[3,,17], type="l", col=1, lty=2)
      panel.text(x=1.2, y=.8, label="-y")
      })

# cairo_pdf("./figures/ICAME40_results_1_syllable.eps", width=7, height=2)
print(plot_1_syllable, position=c(0,0,.36,1), more=T)
print(plot_2_syllable, position=c(.34, 0, .67, 1), more=T)
print(plot_3_syllable, position=c(.67, 0, 1, 1))
# dev.off()

# individual conditions
histogram(data_icame_no_participles$log_rate_z)


# 1 syllable
#--------------------------------------------------------------------------------------------------------------

# cluster
plot_1_cluster <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="1_cluster"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,1], 
                      rev(estimates_array[6,,1])), 
                  y=c(estimates_array[1,,1], 
                      rev(estimates_array[5,,1])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,1], 
                      rev(estimates_array[6,,1])), 
                  y=c(estimates_array[2,,1], 
                      rev(estimates_array[4,,1])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points(x=estimates_array[6,,1], y=estimates_array[3,,1], 
                 type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="1_cluster"))), adj=0, cex=.8, col="darkgrey")
  })

# l
plot_1_l <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="1_l"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,2], rev(estimates_array[6,,2])), y=c(estimates_array[1,,2], rev(estimates_array[5,,2])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,2], rev(estimates_array[6,,2])), y=c(estimates_array[2,,2], rev(estimates_array[4,,2])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points(x=estimates_array[6,,2], y=estimates_array[3,,2], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="1_l"))), adj=0, cex=.8, col="darkgrey")
  })

# other
plot_1_other <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="1_other"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,3], rev(estimates_array[6,,3])), y=c(estimates_array[1,,3], rev(estimates_array[5,,3])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,3], rev(estimates_array[6,,3])), y=c(estimates_array[2,,3], rev(estimates_array[4,,3])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points( x=  estimates_array[6,,3], y=  estimates_array[3,,3], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="1_other"))), adj=0, cex=.8, col="darkgrey")
  })

# r
plot_1_r <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="1_r"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,4], rev(estimates_array[6,,4])), y=c(estimates_array[1,,4], rev(estimates_array[5,,4])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,4], rev(estimates_array[6,,4])), y=c(estimates_array[2,,4], rev(estimates_array[4,,4])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points( x=  estimates_array[6,,4], y=  estimates_array[3,,4], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="1_r"))), adj=0, cex=.8, col="darkgrey")
  })


# 2 syllables
#--------------------------------------------------------------------------------------------------------------


# cluster
plot_2_cluster <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="2_cluster"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,5], rev(estimates_array[6,,5])), y=c(estimates_array[1,,5], rev(estimates_array[5,,5])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,5], rev(estimates_array[6,,5])), y=c(estimates_array[2,,5], rev(estimates_array[4,,5])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points(x=   estimates_array[6,,5], y=  estimates_array[3,,5], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="2_cluster"))), adj=0, cex=.8, col="darkgrey")
  })

# l
plot_2_l <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="2_l"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,6], rev(estimates_array[6,,6])), y=c(estimates_array[1,,6], rev(estimates_array[5,,6])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,6], rev(estimates_array[6,,6])), y=c(estimates_array[2,,6], rev(estimates_array[4,,6])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points(x=   estimates_array[6,,6], y=  estimates_array[3,,6], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="2_l"))), adj=0, cex=.8, col="darkgrey")
  })

# ly
plot_2_ly <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="2_ly"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,7], rev(estimates_array[6,,7])), y=c(estimates_array[1,,7], rev(estimates_array[5,,7])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,7], rev(estimates_array[6,,7])), y=c(estimates_array[2,,7], rev(estimates_array[4,,7])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points(x=   estimates_array[6,,7], y=  estimates_array[3,,7], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="2_ly"))), adj=0, cex=.8, col="darkgrey")
  })

# other
plot_2_other <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="2_other"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,8], rev(estimates_array[6,,8])), y=c(estimates_array[1,,8], rev(estimates_array[5,,8])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,8], rev(estimates_array[6,,8])), y=c(estimates_array[2,,8], rev(estimates_array[4,,8])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points( x=  estimates_array[6,,8], y=  estimates_array[3,,8], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="2_other"))), adj=0, cex=.8, col="darkgrey")
  })

# r
plot_2_r <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="2_r"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,9], rev(estimates_array[6,,9])), y=c(estimates_array[1,,9], rev(estimates_array[5,,9])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,9], rev(estimates_array[6,,9])), y=c(estimates_array[2,,9], rev(estimates_array[4,,9])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points( x=  estimates_array[6,,9], y=  estimates_array[3,,9], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="2_r"))), adj=0, cex=.8, col="darkgrey")
  })

# y
plot_2_y <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="2_y"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,10], rev(estimates_array[6,,10])), y=c(estimates_array[1,,10], rev(estimates_array[5,,10])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,10], rev(estimates_array[6,,10])), y=c(estimates_array[2,,10], rev(estimates_array[4,,10])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points( x=  estimates_array[6,,10], y=  estimates_array[3,,10], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="2_y"))), adj=0, cex=.8, col="darkgrey")
  })


# 3+ syllables
#--------------------------------------------------------------------------------------------------------------


# cluster
plot_3_cluster <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="3+_cluster"), xlim=c(-1.2, 4.9), ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,11], rev(estimates_array[6,,11])), y=c(estimates_array[1,,11], rev(estimates_array[5,,11])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,11], rev(estimates_array[6,,11])), y=c(estimates_array[2,,11], rev(estimates_array[4,,11])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points( x=  estimates_array[6,,11], y=  estimates_array[3,,11], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="3+_cluster"))), adj=0, cex=.8, col="darkgrey")
  })

# l
plot_3_l <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="3+_l"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,12], rev(estimates_array[6,,12])), y=c(estimates_array[1,,12], rev(estimates_array[5,,12])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,12], rev(estimates_array[6,,12])), y=c(estimates_array[2,,12], rev(estimates_array[4,,12])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points( x=  estimates_array[6,,12], y=  estimates_array[3,,12], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="3+_l"))), adj=0, cex=.8, col="darkgrey")
  })

# ly
plot_3_ly <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="3+_ly"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,13], rev(estimates_array[6,,13])), y=c(estimates_array[1,,13], rev(estimates_array[5,,13])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,13], rev(estimates_array[6,,13])), y=c(estimates_array[2,,13], rev(estimates_array[4,,13])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points( x=  estimates_array[6,,13], y=  estimates_array[3,,13], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="3+_ly"))), adj=0, cex=.8, col="darkgrey")
  })

# other
plot_3_other <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="3+_other"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,14], rev(estimates_array[6,,14])), y=c(estimates_array[1,,14], rev(estimates_array[5,,14])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,14], rev(estimates_array[6,,14])), y=c(estimates_array[2,,14], rev(estimates_array[4,,14])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points( x=  estimates_array[6,,14], y=  estimates_array[3,,14], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="3+_other"))), adj=0, cex=.8, col="darkgrey")
  })

# r
plot_3_r <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="3+_r"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,15], rev(estimates_array[6,,15])), y=c(estimates_array[1,,15], rev(estimates_array[5,,15])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,15], rev(estimates_array[6,,15])), y=c(estimates_array[2,,15], rev(estimates_array[4,,15])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points( x=  estimates_array[6,,15], y=  estimates_array[3,,15], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="3+_r"))), adj=0, cex=.8, col="darkgrey")
  })

# ry
plot_3_ry <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="3+_ry"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  #par.settings=my_settings, axis=axis_L,
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,16], rev(estimates_array[6,,16])), y=c(estimates_array[1,,16], rev(estimates_array[5,,16])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,16], rev(estimates_array[6,,16])), y=c(estimates_array[2,,16], rev(estimates_array[4,,16])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points( x=  estimates_array[6,,16], y=  estimates_array[3,,16], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="3+_ry"))), adj=0, cex=.8, col="darkgrey")
  })

# y
plot_3_y <- xyplot(
  jitter(est_prop_more, amount=.025) ~ jitter(log_rate_z, amount=.09), data=subset(data_icame_no_participles, combinations=="3+_y"), xlim=c(-1.2, 4.9),ylim=c(-.05, 1.05),
  col=1, alpha=.33, cex=.5, 
  # par.settings=my_settings, axis=axis_L, 
  xlab=" ", ylab=" ",
  scales=list(x=list(at=(log(c(0.01, 0.1, 1, 10, 100)/1000000)-mean(data_icame_no_participles$log_est_grad_pmw)) / sd(data_icame_no_participles$log_est_grad_pmw), 
                     label=c(0.01, 0.1, 1, 10, 100))),
  panel=function(x,y, subscripts,...){
    panel.polygon(x=c(estimates_array[6,,17], rev(estimates_array[6,,17])), y=c(estimates_array[1,,17], rev(estimates_array[5,,17])), border=F, col="grey90")
    panel.polygon(x=c(estimates_array[6,,17], rev(estimates_array[6,,17])), y=c(estimates_array[2,,17], rev(estimates_array[4,,17])), border=F, col="grey80")
    panel.xyplot(x,y,type="smooth", lty=2, col="grey50")
    panel.points( x=  estimates_array[6,,17], y=  estimates_array[3,,17], type="l", col=1)
    panel.xyplot(x,y,...)
    panel.text(x=-1.2, y=1.15, label=paste0("n = ", nrow(subset(data_icame_no_participles, combinations=="3+_y"))), adj=0, cex=.8, col="darkgrey")
  })



comb = c(plot_3_ly, plot_3_y, plot_3_other, plot_3_l, plot_3_cluster, plot_3_r, plot_3_ry, 
         plot_2_ly, plot_2_y, plot_2_other, plot_2_l, plot_2_cluster, plot_2_r,   
                              plot_1_other, plot_1_l, plot_1_cluster, plot_1_r)


comb = update(comb, layout=c(7,3), between=list(x=.5, y=1.5),
              par.settings=list(layout.heights=list(top.padding=3)))
comb = update(comb, skip=c(F,F,F,F,F,F,F,  F,F,F,F,F,F,T,  T,T,F,F,F,F,T ))
comb

# cairo_pdf("./figures/results_phonological_patterns_details.eps", width=9, height=4.5)
print(comb)
# dev.off()




# Candidates for emergent schemas
#--------------------------------

plot1 = histogram(~log_est_grad_pmw, subset(data_icame, morphology_detail == "ly_suffix"), nint=30, col="grey", 
                  # par.settings=temp_my_settings, axis=axis_L, 
                  type="count",
                  scales=list(x=list(at=log(c(1, 10, 100, 1000, 10000)/100000000), label=c(0.01, 0.1, 1, 10, 100)),
                              y=list(at=c(0, 50, 100, 150))), xlim=log(c(0.006, 200)/1000000), ylim=c(0, 135), 
                  xlab="Token frequency", ylab="Number of adjectives")

# cairo_pdf(file = "./figures/schema_ly_token_frequency_of_types.pdf", height=2, width=2)
print(plot1)
# dev.off()

plot1 = histogram(~log_est_grad_pmw, subset(data_icame, morphology_detail == "y_suffix"), nint=30, col="grey", 
                  # par.settings=temp_my_settings, axis=axis_L, 
                  type="count",
          scales=list(x=list(at=log(c(1, 10, 100, 1000, 10000)/100000000), label=c(0.01, 0.1, 1, 10, 100)),
                      y=list(at=c(0, 50, 100, 150))), xlim=log(c(0.006, 200)/1000000), ylim=c(0, 135), 
          xlab="Token frequency", ylab="Number of adjectives")

# cairo_pdf(file = "./figures/schema_y_token_frequency_of_types.pdf", height=2, width=2)
print(plot1)
# dev.off()

