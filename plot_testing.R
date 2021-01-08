###############################################################################
## Early testing of plots and models. #########################################
###############################################################################
library(corrplot)
library('reshape2')

with(df_cleaned_,cor(df_cleaned_[11:22]))
with(df_cleaned_,corrplot(cor(df_cleaned_[11:22]), method="number"))

Nibor_correlations = sort(sapply(names(df_cleaned_[11:22]),function(bank)
  summary(lm(as.formula(paste("I(Nibor6m - NOKOIS6m) ~",bank)), data = df_cleaned_))$adj.r.squared))
barplot(sqrt(Nibor_correlations))

Euribor_correlations = sort(sapply(names(df_cleaned_[11:22]),function(bank)
  summary(lm(as.formula(paste("I(Euribor6m - OIS6m) ~",bank)), data = df_cleaned_))$adj.r.squared))
barplot(sqrt(Euribor_correlations))

summary(lm(I(Euribor6m - OIS6m) ~ I(apply(df_cleaned_[,c("DZ.Bank","Swedbank")],1,mean)), 
           data = df_cleaned_))$adj.r.squared

summary(lm(I(Nibor6m - NOKOIS6m) ~ I(apply(df_cleaned_[,c("DZ.Bank","DNB")],1,mean)), 
           data = df_cleaned_))$adj.r.squared

summary(lm(I(Nibor6m - NOKOIS6m) ~ I(apply(df_cleaned_[11:16],1,min)), 
           data = df_cleaned_))$adj.r.squared

summary(lm(I(Nibor6m - NOKOIS6m) ~ I(apply(df_cleaned_[11:16],1,mean)), 
           data = df_cleaned_))$adj.r.squared

names(df)[11:16][apply(df_cleaned_[17:22],1,which.min)]

melted = melt(df_cleaned_[,c(1,11:16)], id.vars="Date")
ggplot() + geom_line(data=melted, aes(x=Date, y=I(log(value)), group=variable, colour=variable))

# Forward stepwise regression
Nibor_spread = (df_cleaned_$Nibor6m - df_cleaned_$NOKOIS6m)
Euribor_spread = (df_cleaned_$Euribor6m - df_cleaned_$OIS6m)
min.model = lm(Euribor ~ 1, data = df_cleaned_[11:22])
max.model = formula(lm(Euribor_spread ~ ., data = df_cleaned_[11:22]))
fwd.model = step(min.model, direction='forward', scope=max.model)

cds_pca = prcomp(log(df_cleaned_[,11:22]),
                 center = TRUE,
                 scale. = TRUE) 
###############################################################################
# Lasso and ridge #############################################################
###############################################################################
library(glmnet)
library(corrplot)

CDS_inner = apply((df_cleaned_[,11:16]),1,mean)
summary(lm(nois_diff ~ DZ.Bank + CDS_inner, data = df_cleaned_))$adj.r.squared
summary(lm(nois_diff ~ DZ.Bank + Nordea, data = df_cleaned_))$adj.r.squared
summary(lm(nois_diff ~ DZ.Bank + DNB, data = df_cleaned_))
coefs = coef(lm(nois_diff ~ DZ.Bank + DNB, data = df_cleaned_))
summary(gam(nois_diff ~ s(DZ.Bank) + s(DNB), data = df_cleaned_))
plot(-0.3792+0.007*df_cleaned_$DZ.Bank+0.00378*df_cleaned_$DNB,nois_diff)

# We first run lasso with NOK. 
nois_diff = with(df_cleaned_,Nibor6m - NOKOIS6m)

f_nois = as.formula(nois_diff ~ .)
mod_inner_nois = model.matrix(f_nois, df_cleaned_[11:16])[, -1]
mod_outer_nois = model.matrix(f_nois, df_cleaned_[17:22])[, -1]
mod_complete_nois = model.matrix(f_nois, df_cleaned_[11:22])[, -1]

lasso_complete_nois = glmnet(mod_complete_nois, nois_diff, alpha=1)
cv_lasso_complete_nois = cv.glmnet(mod_complete_nois, nois_diff, alpha=1)
coef(lasso_complete_nois, s = cv_lasso_complete_nois$lambda.1se)
plot(lasso_complete_nois)

# And then with euro.
eois_diff = with(df_cleaned_,Euribor6m - OIS6m)

f_eois = as.formula(eois_diff ~ .)
mod_inner_eois = model.matrix(f_eois, df_cleaned_[11:16])[, -1]
mod_outer_eois = model.matrix(f_eois, df_cleaned_[17:22])[, -1]
mod_complete_eois = model.matrix(f_eois, df_cleaned_[11:22])[, -1]

lasso_complete_eois = glmnet(mod_complete_eois, eois_diff, alpha=0.5)
cv_lasso_complete_eois = cv.glmnet(mod_complete_eois, eois_diff, alpha=0.5)
coef(lasso_complete_eois, s = cv_lasso_complete_nois$lambda.1se)
plot(lasso_complete_eois)

###############################################################################
# Colored gg-plot, Nibor ######################################################
###############################################################################
library("lubridate")
source("multiplot.R")
ndf_cleaned_ = ndf_cleaned[ndf_cleaned$Date>as.Date("2014-01-01"),]
col_date = decimal_date(ndf_cleaned_$Date)

pretty_plot = function(bank, subsetDate, guide = FALSE, title=NULL) {
	sndf = subset(ndf,DecimalDate > subsetDate)
	plotobj = ggplot(sndf, aes(x = eval(parse(text=bank)), 
					   y = Spread, colour=sndf$DecimalDate)) +
  			geom_point(shape = 18,size=1.6) +    
  			geom_smooth(method="lm",col="purple",size=0.5) +
			scale_colour_gradient2(low="red", high="blue", mid="pink",
      		                   midpoint=median(sndf$DecimalDate), name = "Date",
        		    	             guide = ifelse(guide,"colourbar",FALSE))

	if (!is.null(title)) {
  		plotobj = plotobj + labs(x = paste(bank, " CDS"), y = "Nibor-OIS 6m spread",
				     title = title) 
	}

	else {
		plotobj + labs(x = paste(bank, " CDS"), y = "Nibor-OIS 6m spread")
	}

}

subsetDate = "2014-01-01"
plot7 = ggplot(subset(ndf,DecimalDate > subsetDate), aes(x = CDSMean, y = Spread, 
                       colour = DecimalDate, label = Date)) +
  geom_point(shape = 20,size = 1/8*sqrt(subset(ndf,DecimalDate > subsetDate)$CDSVar)) +    
  geom_smooth(method="loess",size=0.5) +
  geom_smooth(method="lm",col="purple",size=0.5) +
  scale_colour_gradient2(low="red", high="blue", mid="green",
                         midpoint=median(subset(ndf,DecimalDate > subsetDate)$DecimalDate), name = "Date") +
  labs(x = "Mean CDS", y = "Nibor-OIS spread, 6m", 
       title = "Nibor-OIS spread, Nibor banks")

CDS_inner_median = apply((ndf_cleaned_[,nor_banks]),1,median)
CDS_inner_variance = apply((ndf_cleaned_[,nor_banks]),1,var)

plot8 = ggplot(ndf_cleaned_, aes(x = CDS_inner_median, y = I(Nibor6m - NOKOIS6m), 
                               colour = col_date, label = Date)) +
  geom_point(shape = 20,size = 1/2*sqrt(CDS_inner_variance)) +    
  #geom_text(aes(label=ifelse(col_date>2016,as.character(Date),'')),hjust=0,vjust=0) +
  geom_smooth(method="loess",size=0.5) +
  geom_smooth(method="lm",col="purple",size=0.5) +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date",
                         guide=FALSE) +
  labs(x = "Median CDS", y = "Nibor-OIS spread, 6m")

multiplot(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,cols=3)


# Bootstrapping

boot_couples = data.frame(cbind(Spread=df_cleaned$Nibor6m - df_cleaned$NOKOIS6m,
                               Median=apply((ndf_cleaned[,9:14]),1,median)))


m = dim(boot_couples)[1]
inc = 200
ks = (1:240)*5
res = sapply(ks,function(k) replicate(100,coef(lm(Spread ~ Median, 
                             data=boot_couples[sample(k:(k+inc),inc,replace=TRUE),])))[2,])

res_means = colMeans(res)
res_ci = apply(res,2,function(x) quantile(x,c(0.025,0.975)))

df_cleaned_ = df_cleaned[df_cleaned$Date>as.Date("2010-01-01"),]
col_date = date_decimal(decimal_date(df_cleaned_$Date) + inc/365)[ks]


library("Hmisc")
errbar(col_date, res_means, res_ci[2,], res_ci[1,], 
       add=F, pch=1, cap=.015, bty="l")
abline(h=0,lty=2,col="grey")



ggplot(df_cleaned_, aes(x = DZ.Bank, y = I(Euribor6m - OIS6m), colour=col_date)) +
  geom_point(shape = 18) +    
  geom_smooth(method="loess") +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date")


ggplot(df_cleaned_, aes(x = DNB, y = I(Euribor6m - OIS6m), colour=col_date)) +
  geom_point(shape = 18) +    
  geom_smooth(method="loess") +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date")

###############################################################################
# Colored gg-plot, Euribor ####################################################
###############################################################################
library("lubridate")
source("multiplot.R")

df_cleaned_ = df_cleaned[df_cleaned$Date>as.Date("2013-01-01"),]
col_date = decimal_date(df_cleaned_$Date)

plot1 = ggplot(df_cleaned_, aes(x = Credit.Agricole, y = I(Euribor6m - OIS6m), colour=col_date)) +
  geom_point(shape = 18,size=1.6) +    
  #geom_smooth(method="loess",size=0.5) +
  geom_smooth(method="lm",col="purple",size=0.5) +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date",
                         guide = FALSE) +
  labs(x = "Credit Agricole CDS", y = "Euribor-OIS 6m spread")

plot2 = ggplot(df_cleaned_, aes(x = HSBC, y = I(Euribor6m - OIS6m), colour=col_date)) +
  geom_point(shape = 18,size=1.6) +    
  #geom_smooth(method="loess",size=0.5) +
  geom_smooth(method="lm",col="purple",size=0.5) +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date",
                         guide = FALSE) +
  labs(x = "HSBC CDS", y = "Euribor-OIS 6m spread")

plot3 = ggplot(df_cleaned_, aes(x = DZ.Bank, y = I(Euribor6m - OIS6m), colour=col_date)) +
  geom_point(shape = 18,size=1.6) +    
  #geom_smooth(method="loess",size=0.5) +
  geom_smooth(method="lm",col="purple",size=0.5) +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date",
                         guide = FALSE) +
  labs(x = "DZ.Bank CDS", y = "Euribor-OIS spread")

plot4 = ggplot(df_cleaned_, aes(x = SG, y = I(Euribor6m - OIS6m), colour=col_date)) +
  geom_point(shape = 18,size=1.6) +    
  #geom_smooth(method="loess",size=0.5) +
  geom_smooth(method="lm",col="purple",size=0.5) +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date",
                         guide = FALSE) +
  labs(x = "SG CDS", y = "Euribor-OIS spread")

plot5 = ggplot(df_cleaned_, aes(x = Barclays, y = I(Euribor6m - OIS6m), colour=col_date)) +
  geom_point(shape = 18,size=1.6) +    
  #geom_smooth(method="loess",size=0.5) +
  geom_smooth(method="lm",col="purple",size=0.5) +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date",
                         guide = FALSE) +
  labs(x = "Barclays CDS", y = "Euribor-OIS spread")

plot6 = ggplot(df_cleaned_, aes(x = Deutsche.Bank, y = I(Euribor6m - OIS6m), colour=col_date)) +
  geom_point(shape = 18,size=1.6) +    
  #geom_smooth(method="loess",size=0.5) +
  geom_smooth(method="lm",col="purple",size=0.5) +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date",
                         guide = FALSE) +
  labs(x = "Deutsche Bank CDS", y = "Euribor-OIS spread")

CDS_inner_mean = apply((df_cleaned_[,19:24]),1,mean)
CDS_inner_variance = apply((df_cleaned_[,19:24]),1,var)

plot7 = ggplot(df_cleaned_, aes(x = CDS_inner_mean, y = I(Euribor6m - OIS6m), 
                               colour = col_date, label = Date)) +
  geom_point(shape = 20,size = 1/8*sqrt(CDS_inner_variance)) +    
  #geom_text(aes(label=ifelse(col_date>2016,as.character(Date),'')),hjust=0,vjust=0) +
  geom_smooth(method="loess",size=0.5) +
  geom_smooth(method="lm",col="purple",size=0.5) +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date") +
  labs(x = "Mean CDS", y = "Euribor-OIS spread, 6m", 
       title = "Euribor-OIS spread, Euribor banks")

CDS_inner_median = apply((df_cleaned_[,19:24]),1,median)
CDS_inner_variance = apply((df_cleaned_[,19:24]),1,var)

plot8 = ggplot(df_cleaned_, aes(x = CDS_inner_median, y = I(Euribor6m - OIS6m), 
                               colour = col_date, label = Date)) +
  geom_point(shape = 20,size = 1/8*sqrt(CDS_inner_variance)) +    
  #geom_text(aes(label=ifelse(col_date>2016,as.character(Date),'')),hjust=0,vjust=0) +
  geom_smooth(method="loess",size=0.5) +
  geom_smooth(method="lm",col="purple",size=0.5) +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date",
                         guide=FALSE) +
  labs(x = "Median CDS", y = "Euribor-OIS spread, 6m")

multiplot(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,cols=3)

###############################################################################
# Plot variances, means and years, Nibor.######################################
###############################################################################
CDS_inner_mean = apply((df_cleaned_[,11:16]),1,mean)
CDS_inner_variance = apply((df_cleaned_[,11:16]),1,var)

ggplot(df_cleaned_, aes(x = CDS_inner_mean, y = I(Nibor6m - NOKOIS6m), 
                       colour = col_date, label = Date)) +
  geom_point(shape = 20,size = 1/8*sqrt(CDS_inner_variance)) +    
  #geom_text(aes(label=ifelse(col_date>2016,as.character(Date),'')),hjust=0,vjust=0) +
  geom_smooth(method="loess",size=0.5) +
  geom_smooth(method="lm",col="purple",size=0.5) +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date") +
  labs(x = "Mean CDS", y = "Nibor-OIS spread, 6m", 
       title = "Nibor-OIS spread, Nibor banks")

CDS_mean = apply((df_cleaned_[,11:22]),1,mean)
CDS_variance = apply((df_cleaned_[,11:22]),1,var)

ggplot(df_cleaned_, aes(x = CDS_mean, y = I(Nibor6m - NOKOIS6m), colour=col_date)) +
  geom_point(shape = 20,size = 1/10*sqrt(CDS_variance)) +    
  geom_smooth(method="loess") +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date") +
  labs(x = "Mean CDS", y = "Nibor-OIS spread, 6m", 
       title = "Nibor-OIS spread, all banks")

###############################################################################
# Plot variances, means and years, Euribor.####################################
###############################################################################
CDS_outer_mean = apply((df_cleaned_[,17:22]),1,min)
CDS_outer_variance = apply((df_cleaned_[,11:22]),1,var)

ggplot(df_cleaned_, aes(x = CDS_outer_mean, y = I(Euribor6m - OIS6m), colour=col_date)) +
  geom_point(shape = 20,size = 1/10*sqrt(CDS_outer_variance)) +    
  geom_smooth(method="loess") +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date") +
  labs(x = "Mean CDS", y = "Euribor-OIS spread, 6m", 
       title = "Euribor-OIS spread, Euribor banks")

CDS_inner_mean = apply((df_cleaned_[,11:16]),1,mean)
CDS_inner_variance = apply((df_cleaned_[,11:16]),1,var)

ggplot(df_cleaned_, aes(x = CDS_inner_mean, y = I(Euribor6m - OIS6m), colour=col_date)) +
  geom_point(shape = 20,size = 1/10*sqrt(CDS_inner_variance)) +    
  geom_smooth(method="loess") +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date") +
  labs(x = "Mean CDS", y = "Euribor-OIS spread, 6m", 
       title = "Euribor-OIS spread, Nibor banks")

CDS_mean = apply((df_cleaned_[,11:22]),1,mean)
CDS_variance = apply((df_cleaned_[,11:22]),1,var)

ggplot(df_cleaned_, aes(x = CDS_mean, y = I(Euribor6m - OIS6m), colour=col_date)) +
  geom_point(shape = 20,size = 1/10*sqrt(CDS_variance)) +    
  geom_smooth(method="loess") +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date") +
  labs(x = "Mean CDS", y = "Euribor-OIS spread, 6m", 
       title = "Euribor-OIS spread, all banks")



###############################################################################
# How well does the Euribor spread explain the Nibor spread?###################
###############################################################################

Nibor_mean = apply((df_cleaned_[,11:16]),1,mean)
Euribor_mean = apply((df_cleaned_[,17:22]),1,mean)
Nibor_Euribor_diff = abs(Nibor_mean-Euribor_mean)

# We use the Nibor_Euribor difference as the bubble size. The relationship does
# not appear linear, but square-rooty.
ggplot(df_cleaned_, aes(x = Euribor_spread, y = Nibor_spread, colour=col_date)) +
  geom_point(shape = 20,size=1/10*Nibor_Euribor_diff) +    
  geom_smooth(method="loess") +
  geom_smooth(method="lm",col="purple") +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date") +
  labs(x = "Euribor-OIS spread, 6m", y = "Nibor-OIS spread, 6m", 
       title = "Nibor spread vs. Euribor spread") +
  geom_abline(slope = 1, col="pink", size=1.2, linetype = 2) +
  coord_fixed(ratio = 1)

# Taking the square root of the Euribor spread works wonders on the AIC.
# > AIC(lm(Nibor_spread~Euribor_spread))
# [1] -1284.39
# > AIC(lm(Nibor_spread~sqrt(Euribor_spread)))
# [1] -1395.085
# > AIC(lm(Nibor_spread~sqrt(Euribor_spread)+Euribor_spread))
# [1] -1393.207
#
# Some of the recent blue points are clear outliers.

ggplot(df_cleaned_, aes(x = sqrt(Euribor_spread), y = Nibor_spread, colour=col_date)) +
  geom_point(shape = 20,size=1/10*Nibor_Euribor_diff) +    
  geom_smooth(method="lm") +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(col_date), name = "Date") +
  labs(x = "Square root of Euribor-OIS spread, 6m", y = "Nibor-OIS spread, 6m", 
       title = "Nibor spread vs. Euribor spread") 