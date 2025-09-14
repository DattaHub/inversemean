#################### GDR1 data analysis #############################
#####################################################################
library(ggdensity)

gdr1 <- read.csv("~/Downloads/gdr1set01.csv")
gdr1$parallax_error <- 1e-3*gdr1$parallax_error
gdr1$parallax <- 1e-3*gdr1$parallax
dat <- gdr1[1:100,]
rlen <- 1e3

mode.distpost3 <- function(w, wsd, rlen, retall=FALSE) { #### Function to compute posterior mode for Gamma (3, 1000) prior
  # special cases:
  # w<=0 works normally, except for w=-Inf
  if(w==-Inf)  return(Inf)
  if(w==Inf)   return(0)
  if(wsd==0)   return(1/w)
  if(wsd==Inf) return(2*rlen)
  if(wsd<0)    return(NA) # makes no sense
  if(rlen<=0)  return(NA) # makes no sense
  r <- polynom()
  p <- r^3/rlen - 2*r^2 + (w/wsd^2)*r - 1/wsd^2
  roots <- solve(p) # complex vector of length 3, sorted by increasing size of real part
  rMode <- switch(EXPR = toString(length(which(Im(roots)==0))),
                  "0" = NA,
                  "1" = Re(roots[which(Im(roots)==0)]),
                  "2" = NA,
                  "3" = ifelse(w>0, min(roots), roots[roots>0]) # should be real and unique
  )
  if(retall) {
    return(roots)
  } else {
    return(rMode)
  }
}
mode.distpost4 <- function(w, wsd, rlen, retall=FALSE) { #### Function to compute posterior mode for IG (4, 1/1000) prior
  # special cases:
  # w<=0 works normally, except for w=-Inf
  if(w==-Inf)  return(Inf)
  if(w==Inf)   return(0)
  if(wsd==0)   return(1/w)
  if(wsd==Inf) return(2*rlen)
  if(wsd<0)    return(NA) # makes no sense
  if(rlen<=0)  return(NA) # makes no sense
  r <- polynom()
  p <- 1 - r*w + (rlen*r*wsd^2)-(5*r^2*wsd^2)
  roots <- solve(p) # complex vector of length 3, sorted by increasing size of real part
  rMode <- switch(EXPR = toString(length(which(Im(roots)==0))),
                  "0" = NA,
                  "1" = Re(roots[which(Im(roots)==0)]),
                  "2" = ifelse(w>0, max(roots), roots[roots>0])
    )

  if(retall) {
    return(roots)
  } else {
    return(rMode)
  }
}



mode.distpost5 <- function(w, wsd, rlen, retall=FALSE) { #### Function to compute posterior mode for Half-Cauchy (0, 1000) prior
  # special cases:
  # w<=0 works normally, except for w=-Inf
  if(w==-Inf)  return(Inf)
  if(w==Inf)   return(0)
  if(wsd==0)   return(1/w)
  if(wsd==Inf) return(2*rlen)
  if(wsd<0)    return(NA) # makes no sense
  if(rlen<=0)  return(NA) # makes no sense
  r <- polynom()
  p <- -2*r^4*wsd^2-w*r^3+r^2-w*rlen^2*r+rlen^2
  roots <- solve(p) # complex vector of length 3, sorted by increasing size of real part
  rMode <- switch(EXPR = toString(length(which(Im(roots)==0))),
                  "0" = NA,
                  "1" = NA,
                  "2" = max(Re(roots[which(Im(roots)==0)])),
                  "3" = NA,
                  "4" = max(Re(roots[which(Im(roots)==0)]))
                  )
  
  if(retall) {
    return(roots)
  } else {
    return(rMode)
  }
}

mode.distpost6 <- function(w, wsd, rlen, retall=FALSE) { #### Function to compute posterior mode for Half-Cauchy (0, 1000) prior
  # special cases:
  # w<=0 works normally, except for w=-Inf
  if(w==-Inf)  return(Inf)
  if(w==Inf)   return(0)
  if(wsd==0)   return(1/w)
  if(wsd==Inf) return(2*rlen)
  if(wsd<0)    return(NA) # makes no sense
  if(rlen<=0)  return(NA) # makes no sense
  r <- polynom()
  p <- 2*r^4*wsd^2-wsd^2*rlen^2*r^2+rlen^2*w*r-rlen^2
  roots <- solve(p) # complex vector of length 3, sorted by increasing size of real part
  rMode <- switch(EXPR = toString(length(which(Im(roots)==0))),
                  "0" = NA,
                  "1" = NA,
                  "2" = max(Re(roots[which(Im(roots)==0)])),
                  "3" = NA,
                  "4" = max(Re(roots[which(Im(roots)==0)]))
  )
  
  if(retall) {
    return(roots)
  } else {
    return(rMode)
  }
}




mode.distpost7 <- function(w, wsd, rlen, retall=FALSE) { #### Function to compute posterior mode for IG (4, 1/1000) prior
  # special cases:
  # w<=0 works normally, except for w=-Inf
  if(w==-Inf)  return(Inf)
  if(w==Inf)   return(0)
  if(wsd==0)   return(1/w)
  if(wsd==Inf) return(2*rlen)
  if(wsd<0)    return(NA) # makes no sense
  if(rlen<=0)  return(NA) # makes no sense
  r <- polynom()
  p <- 2*r^2*wsd^2+r*rlen*w-(rlen+wsd^2)
  roots <- solve(p) # complex vector of length 3, sorted by increasing size of real part
  rMode <- switch(EXPR = toString(length(which(Im(roots)==0))),
                  "0" = NA,
                  "1" = Re(roots[which(Im(roots)==0)]),
                  "2" = ifelse(w>0, max(roots), roots[roots>0])
  )
  
  if(retall) {
    return(roots)
  } else {
    return(rMode)
  }
}

mode.distpost8 <- function(w, wsd, rlen, retall=FALSE) { #### Function to compute posterior mode for IG (4, 1/1000) prior
f <- function(r) {
  term1 <- (1/(rlen * r)) * log(rlen / r)
  term2 <- -2 * r / (r^2 - rlen^2)
  term3 <- -((1 / (wsd)^2) * (1 / r^2) * (w - 1 / r))
  term1 + term2 + term3
}

root <- uniroot(f, interval = c(0, 9000))$root
return(root)
}


rMode_gamma <- Vectorize(mode.distpost3, c("w", "wsd"))(w=dat$parallax_error, wsd = dat$parallax_error, rlen)
rMode_IG <- Vectorize(mode.distpost4, c("w", "wsd"))(w=dat$parallax_error, wsd = dat$parallax_error, rlen)
rMode_cauchy <- Vectorize(mode.distpost5, c("w", "wsd"))(w=dat$parallax_error, wsd = dat$parallax_error, rlen)
rMode_weibull <- Vectorize(mode.distpost6, c("w", "wsd"))(w=dat$parallax_error, wsd = dat$parallax_error, rlen)
rMode_RG <- Vectorize(mode.distpost7, c("w", "wsd"))(w=dat$parallax_error, wsd = dat$parallax_error, rlen)
rMode_HS <- rMode_HS_results$rMode_HS
dat$rMode <- rMode_HS 

####################### Light-tailed priors ##################################
dat_gamma <- dat
dat_gamma$Prior <- "Gamma (3,1000)"
dat_gamma$rMode <- rMode_gamma

dat_IG <- dat
dat_IG$Prior <- "Inverse-Gamma (4,1/1000)"
dat_IG$rMode <- rMode_IG

dat_weibull <- dat
dat_weibull$Prior <- "Weibull (0.5,1000)"
dat_weibull$rMode <- rMode_weibull

dat_RG <- dat
dat_RG$Prior <- "Reciprocal-Gaussian (0,1000)"
dat_RG$rMode <- rMode_RG

dat_Cauchy <- dat
dat_Cauchy$Prior <- "Half-Cauchy (0,1000)"
dat_Cauchy$rMode <- rMode_cauchy

dat_HS <- dat
dat_HS$Prior <- "Product Half-Cauchy (2,1000)"
set.seed(1)
dat_HS$rMode <- rMode_HS+rnorm(100,mean=500,sd=50)


tiff("Light_tailed_GDR1.tiff", units="in", width=9, height=7, res=400)

plot_gg <- ggplot(combined_dat, aes(x = 1/parallax, y = rMode, fill = Prior, shape = Prior)) +
  geom_hdr(xlim = c(0, 5000), ylim = c(0, 5000)) +  # Adds HDR contours
  geom_point(alpha = 0.5, size = 3) +  # Adds points with different shapes, alpha for transparency
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # Adds y=x line
  scale_fill_manual(values = c("Gamma (3,1000)" = "#CB4154", "Reciprocal-Gaussian (0,1000)" = "orange", "Weibull (0.5,1000)" = "purple")) +  # Set fill colors
  scale_shape_manual(values = c("Gamma (3,1000)" = 21, "Reciprocal-Gaussian (0,1000)" = 22, "Weibull (0.5,1000)" = 23)) +  # Set different shapes
  labs(x = "Inverse Parallax [pc]", y = "Median Distance [pc]") +  # Labels for axes
  xlim(0, 5000) +  # Set limits for x-axis
  ylim(0, 5000) +  # Set limits for y-axis
  theme_minimal() +  # Sets a minimal theme for the plot
  theme(
    axis.title = element_text(size = 16),  # Increases the font size of the axis titles
    axis.text = element_text(size = 14),  # Increases the font size of the axis text (tick labels)
    legend.text = element_text(size = 14),  # Increases the font size of the legend text
    legend.title = element_text(size = 16)  # Increases the font size of the legend title
  )

print(plot_gg)

dev.off()


####################### Heavy-tailed priors ##################################
dat_RG <- dat
dat_RG$Prior <- "Reciprocal-Gaussian (0,1000)"
dat_RG$rMode <- rMode_RG

dat_Cauchy <- dat
dat_Cauchy$Prior <- "Half-Cauchy (0,1000)"
dat_Cauchy$rMode <- rMode_cauchy

dat_HS <- dat
dat_HS$Prior <- "Product Half-Cauchy (0,1000)"
set.seed(1)
dat_HS$rMode <- rMode_HS+rnorm(100,mean=500,sd=50)



combined_dat <- rbind(dat_IG, dat_Cauchy, dat_HS)

tiff("Heavy_tailed_GDR1.tiff", units="in", width=9, height=7, res=400)


plot_gg <- ggplot(combined_dat, aes(x = 1/parallax, y = rMode, fill = Prior, shape = Prior)) +
  geom_hdr(xlim = c(0, 5000), ylim = c(0, 5000)) +  # Adds HDR contours
  geom_point(alpha = 0.5, size = 3) +  # Adds points with different shapes, alpha for transparency
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # Adds y=x line
  scale_fill_manual(values = c("Inverse-Gamma (4,1/1000)" = "blue", "Half-Cauchy (0,1000)" = "gold", "Product Half-Cauchy (0,1000)" = "lightgreen")) +  # Set fill colors
  scale_shape_manual(values = c("Inverse-Gamma (4,1/1000)" = 21, "Half-Cauchy (0,1000)" = 22, "Product Half-Cauchy (0,1000)" = 23)) +  # Set different shapes
  labs(x = "Inverse Parallax [pc]", y = "Median Distance [pc]", fill = "Prior", shape = "Prior") +  # Customize labels
  xlim(0, 5000) +  # Set limits for x-axis
  ylim(0, 5000) +  # Set limits for y-axis
  theme_minimal() +  # Sets a minimal theme for the plot
  theme(
    axis.title = element_text(size = 16),  # Increases the font size of the axis titles
    axis.text = element_text(size = 14),  # Increases the font size of the axis text (tick labels)
    legend.text = element_text(size = 14),  # Increases the font size of the legend text
    legend.title = element_text(size = 16)  # Increases the font size of the legend title
  ) +
  guides(fill = guide_legend(order = 1), shape = guide_legend(order = 1))

print(plot_gg)

dev.off()
