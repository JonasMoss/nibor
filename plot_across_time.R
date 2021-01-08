date = decimal_date(ndf$Date)
plot_k = 10
inc = floor(length(date)/plot_k)
stop_dates = ndf$DecimalDate[c(1,(1:(plot_k))*inc)]

rsqs = sapply(1:(plot_k), function(i) {
  model = lm(Spread ~ DNB, 
             data = subset(ndf, stop_dates[i] <= DecimalDate & 
					  stop_dates[i+1] >= DecimalDate))
  c(Rsq=summary(model)$r.squared,Slope=coef(model)[2])
 })

x = sapply(1:(length(stop_dates)-1), function(i) (stop_dates[i+1]+stop_dates[i])*0.5)
x = date_decimal(x)
plot(x = x, y = rsqs[1,], xlab = "Rsq", ylab = "Date")
plot(x = x, y = rsqs[2,], xlab = "Slope", ylab = "Date")
