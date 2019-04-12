#   Ryan Breslawski, rbreslawski@smu.edu
#   April 9, 2019
#   
# Four plots are created: (1) densities
# for the diversity values of each geographic
# unit, as calculated from 20 of the dataset
# permutations, (2) posterior diversity values
# aggregated across all dataset permutations,
# (3) a multinomial posterior sampled from
# one of the dataset permutations, and (4)
# scree plots for the proportions of artifact
# class, by region.
# 
# Output_and_Diagnostics.RData required
# 

# Load libraries
library(ggplot2)
library(patchwork)

# Load data
load("Output_and_Diagnostics.RData")

# Point type names
pt_names <- point_types

# Remove Superior Uplands if working with Phy_s data.
# There are no data for this region.
if("SUPERIOR_UPLAND" %in% Post_Div_Samples$Geog_Zone){
  Post_Div_Samples <- Post_Div_Samples[which(Post_Div_Samples$Geog_Zone!=
                                               "SUPERIOR_UPLAND"),]
}

# Set diversity index for plotting
if(div_arg=="EH"){
  plotdiv <- "italic(E[H])"
  plotdivex <- expression(italic(E[H]))
}else if(div_arg=="H"){
  plotdiv <- "italic(H)"
  plotdivex <- expression(italic(H))
}else if(div_arg=="D"){
  plotdiv <- "1-italic(D)"
  plotdivex <- expression("1-"*italic(D))
}else if(div_arg=="Inv_D"){
  plotdiv <- "1/italic(D)"
  plotdivex <- expression("1/"*italic(D))
}

# FUNCTION: Estimate beta distribution parameters
# from sample mean and variance of posterior samples.
Bpars <- function(samples_vector) {
  
  m <- mean(samples_vector)
  v <- var(samples_vector)
  
  alpha <- ((1-m)/v-1/m)*m^2
  beta <- alpha*(1/m-1)
  
  return(c(alpha, beta))
}

# FUNCTION: Estimate gamma distribution parameters
# from mean and variance of posterior samples.
Gpars <- function(samples_vector){
  
  m <- mean(samples_vector)
  v <- var(samples_vector)
  
  scale <- v/m
  shape <- m/scale
  
  return(c(shape, scale))
  
}

# FUNCTION: Break plot axis lines
ax_lab_break <- function(x,...){
  gsub('\\s','\n',x)
}

# Get geographic zones
Geog_Zones <- unique(Post_Div_Samples$Geog_Zone)

# Reformat geographic zone names
Geog_ZonesF <- sapply(Geog_Zones, function(x){
  
  n <- paste0(toupper(substr(x, 1, 1)),
              tolower(substr(x, 2, nchar(x))))
  
  if(grepl("_", n)){
    n <- strsplit(n, "_")[[1]]
  }else if(grepl("'.'", n)){
    n <- strsplit(n, "'.'")[[1]]
  }else if(grepl(" ", n)){
    n <- strsplit(n, " ")
  }
  
  if(length(n)> 1){
    n <- sapply(n, function(y){
      paste0(toupper(substr(y, 1, 1)),
             tolower(substr(y, 2, nchar(y))))
    })
    return(paste(n, collapse=" "))
  }else{
    return(n)
  }
  
})



####[[[[[[___ _ _ _ .   SUMMARIZE ZONES   . _ _ _ ___]]]]]]###
##############################################################

# Create matrix of diversity values
pm <- sapply(Geog_Zones, function(x){
  Post_Div_Samples$EH[which(Post_Div_Samples$Geog_Zone==x)]
})
pm <- unname(pm)

summary_stats <- lapply(1:ncol(pm), function(x){
  
  HPDI <- round(quantile(pm[, x], probs=c(0.025, 0.975)), 3)
  HPDI <- paste0(HPDI[1], "-", HPDI[2])
  
  mins <- apply(pm, 1, which.min)
  maxs <- apply(pm, 1, which.max)
  
  probs <- lapply(1:ncol(pm), function(y){
    
    if(x==y){
      r <- NA
    }else{
      maxranks <- apply(pm[, c(x, y)], 1, which.max)
      r <- length(which(maxranks==1))/nrow(pm)
    }
    
    return(r)
  })
  
  probs <- as.data.frame(do.call("cbind", probs))
  colnames(probs) <- sapply(Geog_Zones, function(y){
    paste0("ProbDivGreaterThan_", y)
  })
  
  returndf <- data.frame(zone=Geog_Zones[x],
                         post_med=round(median(pm[, x]), 3),
                         post_95HPDI=HPDI,
                         Prob_min_div=length(which(mins==x))/
                                               nrow(pm),
                         Prob_max_div=length(which(maxs==x))/
                                               nrow(pm),
                         stringsAsFactors=FALSE)
  returndf <- cbind(returndf, probs)
  
  return(returndf)
  
})

summary_stats <- do.call("rbind", summary_stats)
rm(pm)
write.csv(summary_stats, file="summary_stats.csv", row.names=FALSE)

# Create  data frame summarizing class probs by zone
class_probs_by_zone <- data.frame(class=point_types,
                                  order_from_plot_bottom=
                                    1:length(point_types),
                                  stringsAsFactors=FALSE)
class_probs <- sapply(1:length(Geog_Zones), function(x){
  apply(sample_model[[2]]$theta[ , x, ], 2, median)
})
class_probs_by_zone <- cbind(class_probs_by_zone,
                             class_probs)
colnames(class_probs_by_zone)[3:(length(Geog_Zones)+2)] <- Geog_Zones
write.csv(class_probs_by_zone, file="classes_by_zone.csv", row.names=FALSE)


####[[[[[[___ _ _ _ . . .   PLOT 1   . . . _ _ _ ___]]]]]]####
##############################################################

# Get maximum permutation ID
Max_iter <- max(Post_Div_Samples$Dataset_Permutation)

# Subset posteriors for up to 30 dataset permutations
if(Max_iter<31){
  PostSub <- Post_Div_Samples
}else{
  PostSub <- 
    Post_Div_Samples[which(Post_Div_Samples$Dataset_Permutation<31),]
}

# Get maximum permutation ID of subset posterior values
Max_iter_sub <- max(PostSub$Dataset_Permutation)

# Generate a list of beta density values for each permutation
# in each geographic zone
bds <- lapply(unique(PostSub$Geog_Zone), function(x){
  
  # Subset PostSub for zone x
  d <- PostSub[which(PostSub$Geog_Zone==x),]
  
  # For each permutation y in zone x
  iterBs <- lapply(1:Max_iter_sub, function(y){
    
    # Subset d for permutation y
    d1 <- d[which(d$Dataset_Permutation==y),]
    
    if(y<10){
      yp <- paste0("000", y)
    }else if(y<100){
      yp <- paste0("00", y)
    }else if(y<1000){
      yp <- paste0("0", y)
    }else{
      yp <- y
    }
    
    # Obtain distribution parameters from posterior
    # samples in permutation y, zone x. If Shannon
    # evenness or Simpson Index, get beta distribution
    # parameters. If Shannon diversity or inverse
    # Simpson, get gamma distribution parameters
    if(div_arg=="EH" | div_arg=="D"){
      
      pars <- Bpars(d1[,which(colnames(d1)==div_arg)])
      # Store sequence of beta density values for
      # permutation y, zone x, in a data frame.
      # Return this data frame.
      df <- data.frame(dens=dbeta(seq(0, 1, 0.00001), 
                                  pars[1], pars[2]),
                       x=seq(0, 1, 0.00001),
                       grp=rep(paste0(yp, x), 100001),
                       iter=rep(y, 100001),
                       stringsAsFactors=FALSE)
    }else{
      pars <- Gpars(d1[,which(colnames(d1)==div_arg)])
      
      # Store sequence of gamma density values for
      # permutation y, zone x, in a data frame.
      df <- data.frame(dens=dgamma(seq(0, 
              max(PostSub[,which(colnames(d1)==div_arg)]), 
              0.0001), shape=pars[1], scale=pars[2]),
              x=seq(0, max(PostSub[,which(colnames(d1)==div_arg)]), 
              0.0001), stringsAsFactors=FALSE)
      df$grp <- rep(paste0(y, x), nrow(df))
      df$iter=rep(y, nrow(df))
    }
    
    # Return data frame.
    return(df)
  })
  
  # Bind all density value data frames across
  # dataset permutation for geographic zone x.
  iterBsDF <- do.call("rbind", iterBs)
  # Add a variable to this data frame that indicates
  # it is associated with geographic zone x.
  iterBsDF$Geog_Zone <- rep(Geog_ZonesF[which(Geog_Zones==x)], 
                            nrow(iterBsDF))
  
  # Return the data frame.
  return(iterBsDF)
})

# Bind list of Beta densities into data frame
bds <- do.call("rbind", bds)
# Remove very large densities
bds <- bds[which(bds$dens!=Inf),]
# Standardize densities for plotting
if(max(bds$dens)<50){
  bds$PlotDens <- bds$dens/(max(bds$dens)*0.75)
}else if(max(bds$dens)<100){
  bds$PlotDens <- bds$dens/(max(bds$dens)*0.5)
}else{
  bds$PlotDens <- bds$dens/(max(bds$dens)*0.25)
}

# Remove very small densities for plotting
bds <- bds[which(bds$PlotDens>0.07),]
# Adjust standardized densities for y-axis placement
bds$PlotDens <- bds$PlotDens + bds$iter

# Arrange factor in bds for plotting
bds$grp <- factor(bds$grp, 
                  levels=unique(bds$grp[order(bds$grp, 
                                              decreasing=TRUE)]))

# Plot densities
plot_permutations <- ggplot(bds, aes(x=x, ymax=PlotDens))+
  annotate("segment", x=-Inf, xend=Inf, y=1:(max(bds$iter)+1),
           yend=1:(max(bds$iter)+1), color="white",
           size=0.4)+
  geom_ribbon(aes(fill=Geog_Zone, group=grp, ymin=iter), 
              color="white", alpha=0.6, size=0.3)+
  labs(x=plotdivex, y="Dataset permutation ID")

  if(div_arg=="EH" | div_arg=="D"){
    plot_permutations <- plot_permutations +
      scale_x_continuous(expand=c(0,0), 
                         breaks=seq(0, 1, 0.05),
                         sec.axis=dup_axis())
    
  }else if(div_arg=="H"){
    plot_permutations <- plot_permutations +
      scale_x_continuous(expand=c(0,0), 
                         breaks=seq(0, 6, 0.5),
                         sec.axis=dup_axis())
    
  }else if(div_arg=="Inv_D"){
    plot_permutations <- plot_permutations +
      scale_x_continuous(expand=c(0,0), 
                         breaks=seq(min(bds$x), 
                         max(bds$x), length.out=7),
                         sec.axis=dup_axis())
  }
  
  if(Max_iter_sub > 19){
    plot_permutations <- plot_permutations +
    scale_y_continuous(expand=c(0.01,0.01), 
                       labels=seq(5, Max_iter_sub, 5),
                       breaks=seq(5, Max_iter_sub, 5)+0.5,
                       limits=c(1, ceiling(max(bds$PlotDens))))+
      theme(panel.grid.major.y=element_blank(), 
            panel.grid.minor.y=element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.major.x=element_line(size=0.4),
            axis.ticks.y=element_blank())
  }else{
    plot_permutations <- plot_permutations +
      scale_y_continuous(expand=c(0.01,0.01), 
                         labels=1:Max_iter_sub,
                         breaks=(1:Max_iter_sub)+0.5,
                         limits=c(1, ceiling(max(bds$PlotDens))))+
      theme(panel.grid.major.y=element_blank(), 
            panel.grid.minor=element_blank(),
            axis.ticks.y=element_blank())
  }

# Display plot
plot(plot_permutations)



####[[[[[[___ _ _ _ . . .   PLOT 2   . . . _ _ _ ___]]]]]]####
##############################################################

# Thin diversity samples for plotting
PostThin <- Post_Div_Samples[seq(20, nrow(Post_Div_Samples), 20),]
# Set y-axis position for plotting
PostThin$y <- sapply(PostThin$Geog_Zone, function(x){
  which(Geog_Zones==x)-1 + rbeta(1, 5, 5)
})


# Get medians and 95% HDIs
sample_stats <- lapply(1:length(Geog_Zones), function(x){
  
  y <- x - 0.5
  inds <- which(Post_Div_Samples$Geog_Zone==Geog_Zones[x])
  m <- quantile(Post_Div_Samples$EH[inds], 
                probs=c(0.5, 0.025, 0.975, 0.25, 0.75))
  return(data.frame(y=y, med=m[1], l=m[2], u=m[3],
                    l2=m[4], u2=m[5], stringsAsFactors=FALSE))
})
sample_stats <- do.call("rbind", sample_stats)


# Plot densities
colnames(PostThin)[which(colnames(PostThin)==div_arg)] <- "x"
plot_aggregated <- ggplot(PostThin, aes(x=x, y=y))+
  annotate("segment", x=-Inf, xend=Inf,  color="white",
           y=seq(0, length(Geog_Zones), 1),
           yend=seq(0, length(Geog_Zones), 1))+
  geom_point(aes(color=Geog_Zone), alpha=0.05)+
  annotate("segment", x=sample_stats$l, xend=sample_stats$u,
           y=sample_stats$y, yend=sample_stats$y,
           color="white", size=0.5)+
  annotate("segment", x=sample_stats$l2, xend=sample_stats$u2,
           y=sample_stats$y, yend=sample_stats$y,
           color="white", size=1.25)+
  annotate("point", x=sample_stats$med, y=sample_stats$y, 
           color="white", size=2.25)+
  labs(x=plotdivex)+
  scale_y_continuous(expand=c(0.01,0.01), 
                     labels=ax_lab_break(Geog_ZonesF),
                     breaks=1:length(Geog_Zones)-0.5)
  
  if(div_arg=="EH" | div_arg=="D"){
    plot_aggregated <- plot_aggregated +
      scale_x_continuous(expand=c(0.02,0.02),
                         sec.axis=dup_axis())
    
  }else if(div_arg=="H"){
    plot_aggregated <- plot_aggregated +
      scale_x_continuous(expand=c(0.02,0.02), 
                         breaks=seq(0, 6, 0.5),
                         sec.axis=dup_axis())
    
  }else if(div_arg=="Inv_D"){
    plot_aggregated <- plot_aggregated +
      scale_x_continuous(expand=c(0.02,0.02), 
                         breaks=seq(min(PostThin$x), 
                         max(PostThin$x), length.out=7),
                         sec.axis=dup_axis())
  }

  plot_aggregated <- plot_aggregated +
    coord_flip()+
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          legend.position="none",
          axis.ticks.x=element_blank(), 
          axis.title.x=element_blank())

# Display plot
plot(plot_aggregated)


####[[[[[[___ _ _ _ . . .   PLOT 3   . . . _ _ _ ___]]]]]]####
##############################################################

# Number of point types
point_types <- ncol(sample_model[[2]]$theta[1, , ])

# Iterations to extract
i <- 500

# Get point diversity for each zone
ze <- lapply(1:length(Geog_Zones), function(z){
  
  # Extract posterior class probs for zone z
  dxz <- sample_model[[2]]$theta[1:i, z, ]

  ### Calculate diversity across posterior samples
  He <- sapply(1:nrow(dxz), function(i) calc_div(dxz[i,]))
  
  # Create list of data frames for plotting geometry
  df <- lapply(1:ncol(dxz), function(x){
    
    if(x==1){
      ymin <- rep(0, nrow(dxz))
    }else{
      ymin <- sapply(1:nrow(dxz), 
                     function(v) sum(dxz[v,1:(x-1)]))
    }

    ymax <- dxz[,x] + ymin
    pt <- rep(x, length(ymin))
    
    rdf <- data.frame(pt=pt, x=1:nrow(dxz), 
                      ymin=ymin, ymax=ymax)
    return(rdf)
  })
  
  df <- do.call("rbind", df)
  df$geog_zone <- rep(z, nrow(df))
  
  # Rescale diversity values for plotting
  He_scaled <- 1.1 + ((He - min(He))/(max(He) - min(He)))*0.3
  
  # Store as a data frame with geographic zone and dataset
  # permutation information.
  He_DF <- data.frame(Div_Index=He, x=1:length(He),
                      Div_Scaled=He_scaled,
                      Geog_Zone=rep(z, length(He)),
                      stringsAsFactors=FALSE)
  return(list(df, He_DF))
})

# Extract point type proportions and diversity values
ze_prop <- do.call("rbind", lapply(1:length(ze), function(x){
  ze[[x]][[1]]
}))
ze_H <- do.call("rbind", lapply(1:length(ze), function(x){
  ze[[x]][[2]]
}))

# Set plotting fill colors
c <- c("#f4cf70", "#f95457", "#b50a98", "#004ba8", "#007575")
c <- rep(c, point_types)[1:point_types]


# Loop through geographic zones to construct plot panels
for(j in 1:length(Geog_Zones)){
  
  # Subset data for geographic zone i
  ddata <- ze_H[which(ze_H$Geog_Zone==j),]
  mdata <- ze_prop[which(ze_prop$geog_zone==j),]
  
  # Set y min and max coordinates for diversity
  yminD <- round(min(ddata$Div_Index), 2)
  ymaxD <- round(max(ddata$Div_Index), 2)
  
  # Construct base plot
  plot_multi <- ggplot(data=mdata, aes(x=x, group=factor(pt)))+
    geom_ribbon(aes(fill=factor(pt), ymin=ymin, ymax=ymax))+
    annotate("text", label=Geog_ZonesF[j], x=i/2, y=1.02, vjust=0,
             color="grey", size=6)+
    scale_fill_manual(values=c)+
    annotate("point", x=ddata$x, y=ddata$Div_Scaled, size=0.2)+
    theme(legend.position="none", axis.title=element_blank(),
          axis.line=element_blank(), axis.ticks=element_blank(),
          panel.grid=element_blank(), axis.text=element_blank(),
          panel.background=element_blank())
  
  if(j%%2==0){
    # Finish plot for right panels
    plot_multi <- plot_multi +
      annotate("segment", y=1.1, yend=1.4, x=i*1.006, 
               xend=i*1.006, color="grey")+
      annotate("segment", y=c(seq(0, 1, 0.25), 1.1, 1.4), 
               yend=c(seq(0, 1, 0.25), 1.1, 1.4), 
               x=c(rep(i*1.002, 5), rep(i*1.006, 2)), 
               xend=rep(i*1.016, 7), color="grey")+
      annotate("text", x=i*1.11, y=1.25, angle=270, vjust=0,
               label=plotdiv, parse=TRUE)+
      annotate("text", x=i*1.11, y=0.5, angle=270, vjust=0,
               label="Point class proportion")+
      annotate("text", x=rep(i*1.09, 7), 
               y=c(seq(0, 1, 0.25), 1.1, 1.4),
               label=sprintf(round(c(0, 0.25, 0.50, 0.75, 
                     1.00, yminD, ymaxD), 2), fmt='%#.2f'),
               hjust=1, size=3, vjust=c(0.2, rep(0.4, 3), 
                                        0.8, 0.2, 0.8))+
      scale_x_continuous(expand=c(0,0), limits=c(-0.004*i, i*1.14))+
      scale_y_continuous(expand=c(0.02,0.02))
  }else{
    # Finish plot for left panels
    plot_multi <- plot_multi +
      annotate("segment", y=1.1, yend=1.4, x=-0.006*i, 
              xend=-0.006*i, color="grey")+
      annotate("segment", y=c(seq(0, 1, 0.25), 1.1, 1.4), 
               yend=c(seq(0, 1, 0.25), 1.1, 1.4), 
               x=c(rep(-0.002*i, 5), rep(-0.006*i, 2)), 
               xend=rep(-0.016*i, 7), 
               color="grey")+
      annotate("text", x=-0.11*i, y=1.25, angle=90, vjust=0,
               label=plotdiv, parse=TRUE)+
      annotate("text", x=-0.11*i, y=0.5, angle=90, vjust=0,
               label="Point class proportion")+
      annotate("text", x=rep(-0.09*i, 7), 
               y=c(seq(0, 1, 0.25), 1.1, 1.4),
               label=sprintf(round(c(0, 0.25, 0.50, 0.75, 
               1.00, yminD, ymaxD), 2), fmt='%#.2f'),
               hjust=0, size=3, vjust=c(0.2, rep(0.4, 3), 
                                        0.8, 0.2, 0.8))+
      scale_x_continuous(expand=c(0,0), limits=c(-0.14*i, 1.004*i))+
      scale_y_continuous(expand=c(0.02,0.02))
    
  }

  if(j==1){
    # If first panel, intialize plot...
    plot_multinomial <- plot_multi
  }else{
    # For subsequent panels, add to plot
    plot_multinomial <- plot_multinomial + plot_multi 
  }
  
}

# Format plot columns
plot_multinomial <- plot_multinomial + plot_layout(ncol=2)
# Display plot
plot(plot_multinomial)


####[[[[[[___ _ _ _ . . .   PLOT 4   . . . _ _ _ ___]]]]]]####
##############################################################

# Sort each geographiz zone by median point class proportion, and
# retrieve the median and 95-5% quantile proportion values.
zone_list <- lapply(Geog_Zones, function(x){
  
  x1 <- class_probs_by_zone[,which(colnames(class_probs_by_zone)==x)]
  sorted <- sapply(x1, function(y) which(sort(x1)==y))
  
  z <- which(Geog_Zones==x)
  
  l <- c(0.025, 0.1, 0.175, 0.25, 0.325, 0.4, 0.475)
  u <- c(0.975, 0.9, 0.825, 0.75, 0.675, 0.6, 0.525)
  
  dfl <- lapply(1:length(x1), function(y){
    
    l1 <- sapply(l, function(w){
      quantile(sample_model[[2]]$theta[,z,y], probs=w)
    })
    u1 <- sapply(u, function(w){
      quantile(sample_model[[2]]$theta[,z,y], probs=w)
    })
    
    dfu <- data.frame(class=rep(pt_names[y], length(l)),
                      x=rep(sorted[y], length(l)),
                      interval=c("x95", "x80", "x65", "x50",
                                 "x35", "x20", "x05"),
                      lower=l1, upper=u1, stringsAsFactors=FALSE)
    return(dfu)
  })
  
  dfl <- do.call("rbind", dfl)

  return(dfl) 
})

# Fix max y-axis value across plots
ymax <- max(sapply(zone_list, function(x) max(x$upper)))

# Create a list of ggplots, one for each region
m_plot_list <- lapply(1:length(zone_list), function(x){
  
  xdf <- zone_list[[x]]
  
  medy <- class_probs_by_zone[,which(colnames(class_probs_by_zone)==
                                                Geog_Zones[x])]
  medx <- sapply(medy, function(y) which(sort(medy)==y))
  
  
  p <- ggplot(xdf, aes(x=x, ymin=lower, ymax=upper))+
    geom_ribbon(aes(group=interval), fill="blue",
                alpha=1/7)+
    annotate("point", x=medx, y=medy, color="black",
             fill="white", size=1, shape=21)+
    labs(x="Ordered point class", y="Proportion", 
         title=Geog_ZonesF[x])+
    scale_y_continuous(limits=c(0, ymax))+
    theme(panel.background=element_blank(),
          panel.grid=element_blank(),
          axis.line=element_line(color="black"),
          axis.ticks=element_line(color="black"))
})

# Assemblage list of ggplots into a multi-panel plot
prop_plot <- m_plot_list[[1]]
for(i in 2:length(m_plot_list)){
  prop_plot <- prop_plot + m_plot_list[[i]]
}
prop_plot <- prop_plot + plot_layout(ncol=2)

# Display plot
plot(prop_plot)
  