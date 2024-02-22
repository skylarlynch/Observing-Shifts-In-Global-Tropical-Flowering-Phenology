
library(circglmbayes)
library(coda) #for mcmc tools to use with circglmbayes
library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)

par(mar = c(1, 1, 1, 1))

setwd("C:/Users/skly5321/Downloads/ChapterOne/ChapterOne/Data")


# The link function for a circular GLM is the tan half function:
#   tan_half(x) = tan(x/2)
# So the inverse link function is the inverse tan half

inv_tan_half <- function(x) {
  return(2 * atan(x))
}

# Converting day to circular

day2circ <- function(day_of_year) {
  return(2 * pi * day_of_year / 365 - pi)
}

circ2day <- function(circ) {
  return(365 * (circ + pi) / (2 * pi))
}

#day2circ-leap <- function(day_of_year) {
 # return(2 * pi * day_of_year / 366 - pi)
#}

#For a circular response y on -pi:pi and a single
# continuous predictor x, the model is:

# $$
# y = \beta_0 + inv_tan_half(\beta_1 * x)
# $$ 




traceplots_list <- list()
finalplot_list <- list()
coef_list <- data.frame()
preds_list <- list()
plant_list <- list()
sdat_list <- list()


df <- read.csv("DatasetByLocation.csv")
df$species <- factor(df$species)
index <- df %>% group_by(species) %>%
  summarise(n=n()) %>% filter(n>15)
df <- df %>% filter(species %in% index$species) %>% mutate(Day.of.the.year = as.numeric(Day.of.the.year))

spec <- unique(as.character(df$species))

for(j in 1:length(spec)){
  plant2 <- na.omit(df[df$species == spec[j],])
  year <- plant2$Year
  year_s <- scale(year) #year scaled and centered
  year_center <- attr(year_s, "scaled:center")
  year_scale <- attr(year_s, "scaled:scale")
  plant_list[[j]] <- plant2
  
  
  
  
  circ <- day2circ(plant2$Day.of.the.year)
  sdat <- data.frame(year, year_s, circ)
  sdat_list[[j]] <- sdat 
  nchains <- 4
  chains <- list()
  for (k in 1:nchains ) {
    fit <- circGLM(circ ~ year_s, data=sdat)
    chains[[k]] <- fit$all_chains
  }
  chains <- mcmc.list(chains)
  traceplots_list[[j]] <- chains
  #autocorr.plot(chains, ask=FALSE)
  #gelman.diag(chains[,"b0_chain"])     
  #gelman.diag(chains[,"kp_chain"])     
  #gelman.diag(chains[,"bt_chain"])  
  fit <- circGLM(circ ~ year_s, data=sdat, burnin=200, thin=30, Q=2500)
  #fit_all[spec[j]] <- fit
  samples <- fit$all_chains
  samplesdf <- data.frame(samples)
  names(samplesdf) <- names(samplesdf) %>% 
    gsub("_chain", "", .)
  # samplesdf %>% 
  #   select("b0","kp","bt") %>% 
  #   pivot_longer(cols=everything(), names_to="parameter", values_to="sample_value") %>% 
  #   ggplot() +
  #   geom_histogram(mapping=aes(x=sample_value, y=after_stat(density)), bins=75) +
  #   facet_wrap(vars(parameter), scales="free")
  # 
  # mean(samplesdf$b0)
  # mean(samplesdf$kp)
  # mean(samplesdf$bt)
  
  hpdi <- function (samp, prob = 0.95) {
    vals <- sort(samp)
    nsamp <- length(vals)
    gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
    init <- 1:(nsamp - gap)
    inds <- which.min(vals[init + gap,drop=FALSE] - vals[init, drop=FALSE])
    ans <- cbind(lower=vals[inds], upper=vals[inds + gap])
    return(ans)
  }
  # 
  # hpdi(samplesdf$b0, prob=0.95)
  # hpdi(samplesdf$kp, prob=0.95)
  # hpdi(samplesdf$bt, prob=0.95)
  # 
  # quantile(samplesdf$b0, prob=c(0.025,0.975))
  # quantile(samplesdf$kp, prob=c(0.025,0.975))
  # quantile(samplesdf$bt, prob=c(0.025,0.975))
  # 
  
  year <- seq(from=min(sdat$year), to=max(sdat$year), by=1) 
  year_s <- scale(year, center=year_center, scale=year_scale)
  n <- length(year)
  results <- matrix(NA, nrow=n, ncol=5) 
  colnames(results) <- c("mnmu","mulo95","muhi95","ppdlo95","ppdhi95")
  
  for ( i in 1:n ) {
    
    mu <- samplesdf$b0 + inv_tan_half(samplesdf$bt * year_s[i])
    
    ppd <- rvon_mises(n=length(mu), mu=mu, kappa=samplesdf$kp)
    
    results[i,1] <- mean(mu)
    #results[i,2:3] <- hpdi(mu, prob=0.95)
    results[i,2:3] <- quantile(mu, prob=c(0.025,0.975)) #CPI
    #results[i,4:5] <- hpdi(ppd, prob=0.95)
    results[i,4:5] <- quantile(ppd, prob=c(0.025,0.975)) #CPI
  }
  results <- circ2day(results) #transform to day of year
  preds <- data.frame(year, year_s, results)
  preds_list[[spec[j]]] <- preds
  rm(year, year_s, n, results, mu, ppd) #clean up
  coef_spec <- as.data.frame(coef(fit))
  coef_spec$Species <- spec[j]
  coef_spec$Effects <- rownames(coef_spec)
  coef_list <- rbind(coef_list,coef_spec)
}


df2 <- preds_list
maxx <- as.data.frame(unlist(lapply(lapply(df2,"[", , "mnmu"),max)))
names(maxx) <- "maxs"
maxx$Species <- rownames(maxx)
minn <- as.data.frame(unlist(lapply(lapply(df2,"[", , "mnmu"),min)))
names(minn) <- "mins"
minn$Species <- rownames(minn)
newdf <- left_join(maxx,minn)

newdf$DaysChange <- (newdf$maxs - newdf$mins)


Days <- newdf$DaysChange
mean(Days)

### Test Direction Shift ###

a <- coef_list[which(coef_list$Effects=="year_s"),]
a$sign <- "last"
a[which(a$Estimate<0),"sign"] <- "first"
print(a$Estimate)
b <- a$Estimate
length(b[b>0]) #pos
length(b[b<0]) #neg
length(b[b==0])
bigdf <- left_join(newdf,a)


pos <- b[b>0]
neg <- b[b<0]
absneg <- abs(b[b<0])

var(pos)
var(neg)



# Welch / pos and neg diff
t.test(pos,absneg, var.equal = F)

# One sample / diff from 0
t.test(b,mu=0)

# One sample / diff from 0
t.test(Days,mu=0)

### Plots ###

# Slope Plot
ggplot() +
  geom_point(data=a, mapping=aes(x=Estimate, y=Species)) +
  geom_errorbarh(data=a, mapping=aes(xmin= Estimate-SD , xmax= Estimate+SD, y= Species)) +
  labs(title= "Location name",
       x ="Circular Slope", y = "Species") 


# Date Plot  

ggplot()+
  geom_segment(data = bigdf, mapping=aes(x = mins, xend = maxs, y = Species, yend=Species), arrow = arrow(length=unit(0.1, "inches"),ends = bigdf$sign)) + 
  xlim(0,365) +
  labs(title= "Location name",
       x ="Julian Date", y = "Species") 



for(i in 1:length(spec)){
  finalplot_list[[i]] <- preds_list[[i]] %>%
    ggplot() +
    geom_ribbon(mapping=aes(x=year, ymin=mulo95, ymax=muhi95), alpha=0.2) +
    geom_point(data=plant_list[[i]], 
               mapping=aes(x=Year, y=Day.of.the.year)) +
    geom_line(mapping=aes(x=year, y=mnmu)) +
    geom_line(mapping=aes(x=year, y=ppdlo95), lty=2) +
    geom_line(mapping=aes(x=year, y=ppdhi95), lty=2) + 
    labs(title= spec[i],
         x ="Year", y = "Julian Date")
  ggsave(filename = paste("C:/Users/skly5321/Downloads/ChapterOne/ChapterOne/Plots/location_All_Plots/",spec[i],"_plot.png",sep=""),device = "png")
}

print(finalplot_list)

print(coef_list)

