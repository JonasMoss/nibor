library("ggplot2")
library("lubridate")

nois = read.csv("data/nok_ois.csv")
nois$Date = as.Date(nois$Date,format="%d/%m/%Y")

nibor = read.csv("data/nibor.csv",sep=";")
nibor$Date = as.Date(nibor$Date-min(nibor$Date), origin="1986-01-01")

nibor_CDS = read.csv("data/nibor_panel_CDS.csv",sep=";")
nibor_CDS$Date = as.Date(nibor_CDS$Date-min(nibor_CDS$Date), origin="2010-01-04")

euribor_CDS = read.csv("data/euribor_panel_CDS.csv",sep=";")
euribor_CDS$Date = as.Date(euribor_CDS$Date-min(euribor_CDS$Date), origin="2010-01-04")

euribor = read.csv("data/euribor.csv",sep=";")
euribor$Date = as.Date(euribor$Date-min(euribor$Date), origin="1999-03-16")

eois = read.csv("data/eur_ois.csv",sep=";")
eois$Date = as.Date(eois$Date-min(eois$Date), origin="1999-01-26")

covariates = read.csv("data/covariates.csv")
covariates$Date = as.Date(covariates$Date-min(covariates$Date), origin="1988-06-23")

df = merge(nois[,c("Date","NOKOIS6m")],eois[,c("Date","OIS6m")],key="Date")
df = merge(df,euribor)
df = merge(df,nibor)
df = merge(df,nibor_CDS)
df = merge(df,euribor_CDS)
df_cleaned = df[complete.cases(df),]

eur_banks = c("Credit.Agricole","HSBC","DZ.Bank","SG","Barclays","Deutsche.Bank")

edf = merge(eois[,c("Date","OIS6m")],euribor[,c("Date","Euribor6m")],key="Date")
edf = merge(edf,euribor_CDS)
edf = merge(edf,covariates[,c("Date","VIX")], key = "Date")
edf = edf[complete.cases(edf),]
edf$Mean = rowMeans(edf[,eur_banks])
edf$Median = apply(edf[,eur_banks],1,median)
edf$Variance = apply(edf[,eur_banks],1,var)
edf$Spread = edf$Euribor6m - edf$OIS6m
edf$DecimalDate = decimal_date(edf$Date)

nor_banks = c("DNB","Nordea","Danske.Bank","Svenska.Handelsbanken","Swedbank","Skandinaviske.Enskildabanken")

ndf = merge(nois[,c("Date","NOKOIS6m")],nibor,key="Date")
ndf = merge(ndf,nibor_CDS)
ndf = merge(ndf,covariates[,c("Date","VIX","NOKUSD6m","NOKUSDSpot","Brent")], key = "Date")
ndf = ndf[complete.cases(ndf),]
ndf$Mean = rowMeans(ndf[,nor_banks])
ndf$Median = apply(ndf[,nor_banks],1,median)
ndf$Variance = apply(ndf[,nor_banks],1,var)
ndf$Spread = ndf$Nibor6m - ndf$NOKOIS6m
ndf$TermSpread = ndf$NOKUSD6m - ndf$NOKUSDSpot
ndf$DecimalDate = decimal_date(ndf$Date)
dim(ndf)

summary(lm(Spread ~ Mean, data=ndf))
summary(lm(Spread ~ VIX, data=ndf))
summary(lm(Spread ~ Mean + iTraxx.Senior, data=ndf))
summary(lm(Spread ~ iTraxx.Senior, data=ndf))
summary(lm(Spread ~ Mean + NOKUSD6m, data=ndf))
summary(lm(Spread ~ Mean + NOKUSDSpot, data=ndf))
summary(lm(Spread ~ Mean + I(NOKUSD6m-NOKUSDSpot), data=ndf))
summary(lm(Spread ~ Mean + I(NOKEUR6m-NOKEURSpot), data=ndf))
summary(lm(Spread ~ Mean + Brent + I(NOKUSD6m-NOKUSDSpot), data=ndf))

mod = lm(Spread ~ Mean + I(NOKUSD6m-NOKUSDSpot), data=ndf)
ggplot(ndf, aes(x = coef(mod)[2]*Mean + coef(mod)[3]*(NOKUSD6m-NOKUSDSpot), y = Spread, colour = DecimalDate, label = Date)) +
         geom_point(shape = 20, size = 1/10*sqrt(ndf$Variance)) +    
         scale_colour_gradient2(low="red", high="blue", mid="pink",
                                midpoint=median(ndf$DecimalDate), name = "Date") +
         labs(x = "Mean CDS", y = "Nibor-OIS spread, 6m", 
              title = "Nibor-OIS spread, Nibor banks")


mod = lm(Spread ~ Mean + I(NOKUSD6m-NOKUSDSpot) + Brent, data=ndf)
ggplot(ndf, aes(x = coef(mod)[2]*Mean + coef(mod)[3]*(NOKUSD6m-NOKUSDSpot) + coef(mod)[4]*Brent, 
                y = Spread, colour = DecimalDate, label = Date)) +
  geom_point(shape = 20, size = 1/10*sqrt(ndf$Variance)) +    
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(ndf$DecimalDate), name = "Date") +
  labs(x = "Mean CDS", y = "Nibor-OIS spread, 6m", 
       title = "Nibor-OIS spread, Nibor banks")


ggplot(ndf, aes(x = VIX, y = Spread, colour = DecimalDate, label = Date)) +
  geom_point(shape = 20, size = 1/10*sqrt(ndf$Variance)) +    
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint=median(ndf$DecimalDate), name = "Date") +
  labs(x = "Mean CDS", y = "Nibor-OIS spread, 6m", 
       title = "Nibor-OIS spread, Nibor banks")




