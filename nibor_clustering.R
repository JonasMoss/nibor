###############################################################################
# Testing with VIX. ###########################################################
###############################################################################

k = 10
cols = rainbow(k+1)
s = new_estimator(ndf$VIX,ndf$Spread,k=k,maxIter=50)
pchs = rep(15:20,7)

indices = s$splits

with(ndf,plot(VIX,Spread,type="n",bty="l"))
for (i in (1:k)) with(ndf[indices[i]:(indices[i+1]-1),],points(VIX,Spread,col=cols[i],pch=pchs[i]))
legend("bottomright",
       legend = ndf$Date[indices[1:(k)]],
       pch = pchs,
       col = cols,
       bty = "n") 

###############################################################################
# Testing with NOKUSD. ########################################################
###############################################################################

k = 10
cols = rainbow(k+1)
s = new_estimator(ndf$TermSpread,ndf$Spread,k=k,maxIter=50)
pchs = rep(15:20,7)

indices = s[50,12:(12+(k))]

with(ndf,plot(TermSpread,Spread,type="n",bty="l"))
for (i in (1:k)) with(ndf[indices[i]:(indices[i+1]-1),],points(TermSpread,Spread,col=cols[i],pch=pchs[i]))
legend("bottomright",
       legend = ndf$Date[indices[1:(k)]],
       pch = pchs,
       col = cols,
       bty = "n") 

###############################################################################
# Testing with CDS DNB  #######################################################
###############################################################################

k = 10
cols = rainbow(k+1)
s = new_estimator(ndf$DNB,ndf$Spread,k=k,maxIter=50)
pchs = rep(15:20,7)

indices = s[50,12:(12+(k))]

with(ndf,plot(Median,Spread,type="n",bty="l"))
for (i in (1:k)) with(ndf[indices[i]:(indices[i+1]-1),],points(DNB,Spread,col=cols[i],pch=pchs[i]))
legend("bottomright",
       legend = ndf$Date[indices[1:(k)]],
       pch = pchs,
       col = cols,
       bty = "n") 


###############################################################################
# Testing with CDS median #####################################################
###############################################################################

k = 9
cols = rainbow(k+1)
s = new_estimator(ndf$Median,ndf$Spread,k=k,maxIter=50)
pchs = rep(15:20,7)

indices = s$splits

with(ndf,plot(Median,Spread,type="n",bty="l"))
for (i in (1:k)) with(ndf[indices[i]:(indices[i+1]-1),],points(Median,Spread,col=cols[i],pch=pchs[i]))
legend("bottomright",
       legend = ndf$Date[indices[1:(k)]],
       pch = pchs,
       col = cols,
       bty = "n") 


###############################################################################
# Testing with nothing      ###################################################
###############################################################################

with(ndf,plot(DecimalDate,Spread,type="n",bty="l"))
for (i in (1:k)) with(ndf[indices[i]:(indices[i+1]-1),],points(DecimalDate,Spread,col=cols[i],pch=pchs[i]))
legend("bottomright",
       legend = ndf$Date[indices[1:(k)]],
       pch = pchs,
       col = cols,
       bty = "n") 


###############################################################################
# Testing with Euribor CDS median #############################################
###############################################################################

k = 8
cols = rainbow(k+1)
s = new_estimator(Median,edf$Spread,k=k,maxIter=50)
pchs = rep(15:20,7)

indices = s[50,12:(12+(k))]

with(edf,plot(Median,Spread,type="n",bty="l"))
for (i in (1:k)) with(edf[indices[i]:(indices[i+1]-1),],points(Median,Spread,col=cols[i],pch=pchs[i]))
legend("bottomright",
       legend = edf$Date[indices[1:(k)]],
       pch = pchs,
       col = cols,
       bty = "n") 


###############################################################################
# Testing with nothing (EUR)     ##############################################
###############################################################################

with(edf,plot(DecimalDate,Spread,type="n",bty="l"))
for (i in (1:k)) with(edf[indices[i]:(indices[i+1]-1),],points(DecimalDate,Spread,col=cols[i],pch=pchs[i]))
legend("topright",
       legend = ndf$Date[indices[1:(k)]],
       pch = pchs,
       col = cols,
       bty = "n") 

###############################################################################
# Testing with CDS median, GG #################################################
###############################################################################


prettyplot = function(bank = "Median",k = 10, s = NULL) {
  
  if(is.null(s)) {
    s = new_estimator(ndf[,bank],ndf$Spread,k=k,maxIter=50)
  }
  
  indices = s[50,12:(12+(k))]
  pchs = rep(15:20,7)
  
  group_lengths = sapply(1:k, function(i) indices[i+1] - indices[i])
  groups = unlist(sapply(1:k,function(i) rep(indices[i],group_lengths[i])))
  
  DateGroups = as.factor(ndf$Date[groups]) 
  
  print(ggplot(ndf, aes(x = eval(parse(text=bank)), y = Spread, 
                  colour = DateGroups, shape = DateGroups, size = Variance)) +
    geom_point() +    
    scale_colour_brewer(palette = "Paired", name = "Dates") +
    scale_shape_manual(values = pchs, name = "Dates") +
    guides(size = "none") +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = paste(bank, " CDS"), y = "Nibor-OIS spread, 6m", 
         title = "Nibor-OIS spread"))
}


for (i in (1:k)) {
  subndf = ndf[indices[i]:(indices[i+1]-1),]
  break_indices = c(1,floor((indices[i+1]-indices[i])/4)*(1:4))
  breaks = subndf$DecimalDate[break_indices]
  plot_obj = ggplot(subndf, aes(x = Median, y = Spread, colour = DecimalDate, size = Variance)) +
    geom_point() +    
    guides(size = "none") +
    scale_colour_gradient(low = "red", high = "blue", guide = "colourbar", 
                          name = NULL, breaks = breaks,
                          labels = subndf$Date[break_indices]) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Median CDS", y = "Nibor-OIS spread, 6m", 
         title = paste("Nibor-OIS: ",ndf$Date[indices[i]], " to ", ndf$Date[indices[i+1]-1]))
  print(plot_obj)
}