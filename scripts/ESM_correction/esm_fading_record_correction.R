library(fields)
library(data.table)
library(PEcAn.allometry)

####################################### Read data in #######################################

file.in.normal = ".../ESM.rds"
file.in.fading = ".../ESM_F_survivalonly.rds"
file.in.model  = ".../HF_PREDS_tree_iter_25June2019.RDS"


# --- read in FIA DB query
# read in ESM 
q = as.data.table(readRDS(file.in.normal))
# read in ESM fading record scenario 
w = as.data.table(readRDS(file.in.fading))

sumNA  = function(x) sum(x,na.rm=T)
meanNA = function(x) mean(x,na.rm=T)
q[, PREVTPAsum     := sumNA(start1tpa), by=PLT_CN                 ]  # "startsumTPA"
q[, TPAsum         := sumNA(end1tpa), by=PLT_CN                   ]  # "endsumTPA"
# Compute plot mean diameters
q[, PREVDIAmean    := meanNA(DIA_BEGIN), by=PLT_CN                ]  # "DIAbeginmean"
q[, DIAmean        := meanNA(DIA_END),   by=PLT_CN                ]  # "DIAendmean"

w[, PREVTPAsum     := sumNA(start1tpa), by=PLT_CN                 ]  # "startsumTPA"
w[, TPAsum         := sumNA(end1tpa), by=PLT_CN                   ]  # "endsumTPA"
# Compute plot mean diameters
w[, PREVDIAmean    := meanNA(DIA_BEGIN), by=PLT_CN                ]  # "DIAbeginmean"
w[, DIAmean        := meanNA(DIA_END),   by=PLT_CN                ]  # "DIAendmean"

setnames(q, tolower(names(q)))
q[, plt_cn := as.character(plt_cn)]


setnames(w, tolower(names(w)))
w[, plt_cn := as.character(plt_cn)]


# --- Convert units
q[, diamean.m     := diamean     * 2.54      ] # inches to cm
q[, prevdiamean.m := prevdiamean * 2.54      ]
q[, tpasum.m      := tpasum      / 0.404686  ] # trees acre-1 to trees ha-1 0.404686
q[, prevtpasum.m  := prevtpasum  / 0.404686  ]

w[, diamean.m     := diamean     * 2.54      ]
w[, prevdiamean.m := prevdiamean * 2.54      ]
w[, tpasum.m      := tpasum      / 0.404686  ]
w[, prevtpasum.m  := prevtpasum  / 0.404686  ]

fia_to_paleon   <- read.table(".../fia_to_level3a_v0.5.csv", heade = TRUE, sep="\t")
chojnacky_table <- read.table(".../level3a_to_chojnacky_v0.2.csv", header = TRUE, sep = "\t")

# there isn't any, but filter out plots that don't have FIA spcd
q <- q[!(q$plt_cn %in% unique(q$plt_cn[is.na(q$spcd)])),]
# filter out plots tht don't have matching spp
plots_with_these_spcd <- fia_to_paleon$fia_spcd[!(fia_to_paleon$level3a %in% chojnacky_table$level3a)]
q <- q[!(q$plt_cn %in% unique(q$plt_cn[q$spcd %in% plots_with_these_spcd])), ]

# match it to PalEON PFTs
length(unique(fia_to_paleon$fia_spcd))
length(unique(fia_to_paleon$level3a))
q$paleon <- plyr::mapvalues(q$spcd, fia_to_paleon$fia_spcd, as.character(fia_to_paleon$level3a))

# read in PalEON tree-ring model for Harvard Forest
hf_rec <- readRDS(file.in.model)

years <- 2012:1961

site_id <- 3

hf_census <- hf_rec[hf_rec$model == "Model RW + Census" & hf_rec$plot == site_id, ]
hf_tronly <- hf_rec[hf_rec$model == "Model RW"          & hf_rec$plot == site_id, ]

big_or_small <- function(site_sub, 
                         years = rev(unique(site_sub$year)), 
                         niter = 250){
  
  site_sub$small <- 1
  
  collect_iter <- list()
  for(i in 1:niter){
    
    iter_sub <- site_sub[site_sub$iter == i, ]
    
    for(y in years){
      
      first_year <- iter_sub[iter_sub$year == y, ]
      first_year <- first_year[complete.cases(first_year),]
      
      big_trees <- first_year[first_year$dbh > 20,]
      
      if(nrow(big_trees)==0){
        print(y)
        next
      }
      # if it's a big tree it's always big in the past
      # so that they don't move into the inner circle as they get smaller in the past
      big_tree_ids <- big_trees$tree
      iter_sub$small[iter_sub$tree %in% big_tree_ids] <- 0
    }
    collect_iter[[i]] <- iter_sub
  }
  
  site_sub <- do.call("rbind", collect_iter)
  return(site_sub)
}

hf_tronly <- big_or_small(site_sub = hf_tronly)

inner_circle_area <- (pi*13^2) / 10000
outer_ring_area <- ((pi*20^2) / 10000) - inner_circle_area



########### explore 
# crec : census recrd
# frec : fading record

# in this loop we apply cutoff and clone small trees to the outer circle
frec <- crec <- outsmall <- list()
i <- y <- 1 
for( i in 1:250){
  frec[[i]] <- crec[[i]] <- matrix(NA,ncol= 3, nrow= length(years))
  outsmall[[i]] <- rep(NA, length(years))
  
  for( y in seq_along(years)){
    
    # extract current and previous year from the PalEON tree ring only model
    cur_f <- hf_tronly[hf_tronly$year == years[y] & hf_tronly$iter == i, ]
    
    # extract current and previous year from the PalEON tree ring + census model (validation)
    cur <- hf_census[hf_census$year == years[y] & hf_census$iter == i, ]
    
    cur_f <- cur_f[complete.cases(cur_f),]
    cur <- cur[complete.cases(cur),]
    
    #apply fia cutoff
    cur <- cur[cur$dbh >= 5*2.54, ]
    cur_f <- cur_f[cur_f$dbh >= 5*2.54, ]
    
    #print(sum((cur_f$dist_census > 13 & cur_f$dbh < 20)))
    # cur_f <- cur_f[!(cur_f$dist_census > 13 & cur_f$dbh < 20), ]
    
    iter_df <- cur_f
    small_trees <-  iter_df[(cur_f$dist_census < 13 & cur_f$small == 1),]
    #print(sum((cur_f$dist_census < 13 & cur_f$dbh < 20)))
    #small_trees <-  iter_df[iter_df$small == 1, ]
    x_trees <- ceiling(outer_ring_area* nrow(small_trees) /inner_circle_area)
    #print(x_trees)
    outsmall[[i]][y] <- sum(cur_f$dist_census > 13) - x_trees
    x_ind <- sample(1:nrow(small_trees), x_trees, replace = TRUE)
    
    iter_df <- rbind(iter_df, small_trees[x_ind,])
    
    fdbh <- mean(iter_df$dbh)
    cdbh <- mean(cur$dbh)
    
    fab <- sum(iter_df$ab) * 1e-03 / ((pi*20^2) / 10000)
    cab <- sum(cur$ab)* 1e-03 / ((pi*20^2) / 10000)
    
    fden <- nrow(iter_df) / ((pi*20^2) / 10000)
    cden <- nrow(cur) / ((pi*20^2) / 10000)
    
    frec[[i]][y,] <- cbind(fdbh,fden,fab)
    crec[[i]][y,] <- cbind(cdbh,cden,cab)
  }
  print(i)
}


## fading_arrow_length <- function(...){
##    ...
## }


# this function rescales the 5-yr ESM arrow and calculates AB
update_ab <- function(these_plots_cur, esm, sfc){
  rescaled_mean_dbhs <- rep(NA, nrow(these_plots_cur))
  rescaled_mean_dens <- rep(NA, nrow(these_plots_cur))
  rescaled_mean_agbs <- rep(NA, nrow(these_plots_cur))
  
  #recalculate ab allometrically for rescaled
  for(p in seq_len(nrow(these_plots_cur))){
    
    sub_plot <- esm[esm$plt_cn == these_plots_cur$plt_cn[p], ] 
    
    # first scale esm to per year rate
    sf <- 1/(sub_plot$invyr[1] - sub_plot$prevyr[1])
    new_dbh <- sub_plot$diamean.m[1] + (sub_plot$prevdiamean.m[1] - sub_plot$diamean.m[1]) * sf
    prcnt_dbh <- ((new_dbh- sub_plot$prevdiamean.m[1])/ sub_plot$prevdiamean.m[1])*100
    sub_plot$dia_begin <- sub_plot$dia_begin + (sub_plot$dia_begin * prcnt_dbh/100)
    
    new_den <- sub_plot$tpasum.m[1] + (sub_plot$prevtpasum.m[1] - sub_plot$tpasum.m[1]) * sf
    prcnt_den <- (new_den- sub_plot$prevtpasum.m[1])/ sub_plot$prevtpasum.m[1]*100
    sub_plot$start1tpa <- sub_plot$start1tpa + (sub_plot$start1tpa * prcnt_den/100)
    
    
    # if we want to update the rescaling the arrow with site specific magnitudes (using sfc) do it here
    
    #new_dbh <- sub_plot$diamean[1] + (mean(sub_plot$dia_begin,na.rm=TRUE)-sub_plot$diamean[1])*sfc
    #new_den <- sub_plot$tpasum[1] + (sum(sub_plot$start1tpa,na.rm=TRUE)-sub_plot$tpasum[1])*sfc
    
    #prcnt_dbh <- (new_dbh- mean(sub_plot$dia_begin,na.rm=TRUE))/ mean(sub_plot$dia_begin,na.rm=TRUE)*100
    #prcnt_den <- (new_den- sum(sub_plot$start1tpa,na.rm=TRUE))/ sum(sub_plot$start1tpa,na.rm=TRUE)*100
    
    #sub_plot$dia_begin <- sub_plot$dia_begin + (sub_plot$dia_begin * prcnt_dbh/100)
    #sub_plot$start1tpa <- sub_plot$start1tpa + (sub_plot$start1tpa * prcnt_den/100)
    
    beta0s <- chojnacky_table$beta0[match(sub_plot$paleon, as.character(chojnacky_table$level3a))]
    beta1s <- chojnacky_table$beta1[match(sub_plot$paleon, as.character(chojnacky_table$level3a))]
    
    ### exp(b0 + b1*log(dbh))
    agb           <- sum(exp(beta0s + beta1s*log(sub_plot$dia_begin* 2.54)), na.rm = TRUE)
    agb_Mg        <- agb * 1e-03 # kg to Mg
    agb_Mg_plt    <- agb_Mg / (4*7.3152^2*pi) # Mg/m2
    agb_Mg_ha_plt <- agb_Mg_plt / 1e-04  # Mg/m2 to Mg/ha
    
    
    rescaled_mean_dbhs[p] <- mean(sub_plot$dia_begin * 2.54, na.rm=TRUE)
    rescaled_mean_dens[p] <- sum(sub_plot$start1tpa / 0.404686, na.rm = TRUE)
    rescaled_mean_agbs[p] <- agb_Mg_ha_plt
  }
  
  #corrected_mat <- matrix(apply(cbind(rescaled_mean_dbhs, rescaled_mean_dens, rescaled_mean_agbs), 2, mean), ncol=3)
  corrected_mat <- cbind(rescaled_mean_dbhs, rescaled_mean_dens, rescaled_mean_agbs)
  
  return(corrected_mat)
}

pal_tronly <- matrix(NA,nrow=250,ncol=length(years))
pal_tr_cen <- matrix(NA,nrow=250,ncol=length(years))
radd <- 50
th <- 5
th2 <- 0.5

corrected_ab <- list()

#load ab lut table
load(".../ab_lut.Rdata")

#normalize
ab_lut_norm <- ab_lut
mins <- apply(ab_lut[,2:4], 2, min)
maxs <- apply(ab_lut[,2:4], 2, max)
ab_lut_norm[,2] <- (ab_lut_norm[,2] - mins[1])/ (maxs[1]-mins[1])
ab_lut_norm[,3] <- (ab_lut_norm[,3] - mins[2]) / (maxs[2]-mins[2])
ab_lut_norm[,4] <- (ab_lut_norm[,4] - mins[3]) / (maxs[3]-mins[3])


search_esm <- function(ab_lut, ab_lut_norm, current_pnt, th, mins, maxs, nneigh){
  
  current_pnt_norm <- current_pnt
  current_pnt_norm[,1] <- (current_pnt_norm[,1] - mins[1])/ (maxs[1]-mins[1])
  current_pnt_norm[,2] <- (current_pnt_norm[,2] - mins[2]) / (maxs[2]-mins[2])
  current_pnt_norm[,3] <- (current_pnt_norm[,3] - mins[3]) / (maxs[3]-mins[3])
  
  
  sub_norm <- ab_lut_norm[abs(ab_lut[,4] - current_pnt[1,3]) < th,]
  d <- rdist(sub_norm[,c(2:4)], current_pnt_norm) # Euclidean distance in 3D
  these_plots <- sub_norm[(order(d)[1:(nneigh*1)]),1] # use nneigh nearest-neighbours
  # little jitter
  #these_plots <- these_plots[sample(1:(nneigh*1), nneigh)]
  these_plots <- ab_lut[ab_lut$plt_cn %in% these_plots, ]
  
 
  return(these_plots)
}
years <- 2012:1962
th <- 1
nneigh <- 20
i <- y <- 1 
corrected_ab <- list()
for( i in 1:250){
  
  plot(12:40,seq(100,1200, length.out= length(12:40)),"n",xlab="dbh",ylab="dens",main=paste0("nneigh: ",nneigh, " i: ",i, " r: ", radd))
  arrows(meanB$diameanB,meanB$tpasumB,meanB$prevdiameanB,meanB$prevtpasumB, col="lightgray",
         length=0.065, angle=22,lwd=1) 
  legend("topright", legend =c("paleon tree_ring only", "paleon tr+census", "esm correction"),
         col = c(2,1,4), lty=1)
  corrected_ab[[i]] <- list()
  #corrected_ab[[i]] <- matrix(NA, ncol = 3, nrow = length(years)+1)
  
  for( y in seq_along(years)){
    
    # draw the arrows
    # non-fading arrow, we don't have this at every site
    arrows(crec[[i]][y,1], crec[[i]][y,2], crec[[i]][y+1,1], crec[[i]][y+1,2], col=1,
           length=0.065, angle=22)
    # fading arrow
    arrows(frec[[i]][y,1], frec[[i]][y,2], frec[[i]][y+1,1], frec[[i]][y+1,2], col=2,
           length=0.065, angle=22)
    
    
    if(y == 1){ 
      # in the most recent year there is no fading issue
      # use tr-only model to find the start
      #corrected_ab[[i]][y,]  <- frec[[i]][y,, drop=FALSE]
      
      # find the closest FIA plots to the beginning of the red arrow
      current_pnt <- frec[[i]][y,, drop=FALSE]
      these_plots <- search_esm(ab_lut, ab_lut_norm, current_pnt, th, mins, maxs, nneigh)
      corrected_ab[[i]][[y]] <- these_plots[,2:4]
      
      these_plots_past <- update_ab(these_plots, q, 1)
      corrected_ab[[i]][[y+1]]  <- these_plots_past
      
    }else{
      # use corrected arrow of the precious iteration
      # find the closest FIA plots to the beginning of the red arrow
      current_pnts <- corrected_ab[[i]][[y]]
      #current_pnt <- corrected_ab[[i]][y,,drop=FALSE]
       new_plots <- lapply(seq_len(nneigh), function(x){
         tmp_plots <- search_esm(ab_lut, ab_lut_norm, current_pnts[x,,drop=FALSE], th, mins, maxs, nneigh)
         return(tmp_plots)
       })
       for(nsl in seq_along(new_plots)){
         these_plots_past <- update_ab(new_plots[[nsl]], q, 1)
         new_plots[[nsl]] <- apply(these_plots_past, 2, mean, na.rm=TRUE)
       }
       new_plots <- do.call("rbind",new_plots)

       corrected_ab[[i]][[y+1]] <- new_plots

    }
    
    prev_arr <- apply(corrected_ab[[i]][[y]], 2, mean, na.rm=TRUE)
    next_arr <- apply(corrected_ab[[i]][[y+1]], 2, mean, na.rm=TRUE)
    # corrected arrow
    arrows(prev_arr[1], prev_arr[2], next_arr[1], next_arr[2], 
           col=4, length=0.065, angle=22)
    #arrows(corrected_ab[[i]][y,1], corrected_ab[[i]][y,2], 
    #       corrected_ab[[i]][y+1,1], corrected_ab[[i]][y+1,2], 
    #       col=4, length=0.065, angle=22)
    }
}

i <- 4
plot(corrected_ab[[i]][,3], crec[[i]][,3])
abline(0,1)

corrected_ab[[i]] <- NULL
y.list <- list()
for(l in 1:52){
  esm_correction <- list()
  for(tyi in seq_along(corrected_ab)){
    tmp <- sapply(corrected_ab[[tyi]], function(x) x[,3])
    #tmp <- corrected_ab[[tyi]][,3]
    esm_correction[[tyi]] <- tmp[,l]
  }
  y.list[[l]] <- quantile(unlist(esm_correction),  c(0.025, 0.5, 0.975))
  #y.list[[l]] <- mean(unlist(esm_correction))
}

agbCI_f <- do.call("cbind", y.list)
#agbCI_f <- apply(esm_correction, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
agb_poly <- 1:dim(agbCI_f)[2]

pal_tronly <-lapply(frec, function(x) x[,3])
pal_tronly <- do.call("rbind", pal_tronly)
and_CI_f <- apply(pal_tronly, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)

pal_tr_cen <-lapply(crec, function(x) x[,3])
pal_tr_cen <- do.call("rbind", pal_tr_cen)
and_CI <- apply(pal_tr_cen, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)



plot(agbCI_f[1,], ylim = c(50,  400), xaxt = "n",lwd = 3, xlab = "Years", ylab = "AB (Mg/ha)",
     main = "HF - site 1", type = "n", cex.lab=1.5)
polygon(c(agb_poly, rev(agb_poly)),c((agbCI_f[3,]), rev(agbCI_f[1,])),col=adjustcolor("lightgray",alpha.f=0.5),border="darkgray", lwd =2)
lines(agbCI_f[2,], col = "darkgray", lwd=2)
polygon(c(agb_poly, rev(agb_poly)),c((and_CI_f[3,]), rev(and_CI_f[1,])),col=adjustcolor("lightpink",alpha.f=0.5),border="lightpink", lwd =2)
lines(and_CI_f[2,], col = "lightpink", lwd=2)
polygon(c(agb_poly, rev(agb_poly)),c((and_CI[3,]), rev(and_CI[1,])),col=adjustcolor("lightblue",alpha.f=0.5),border="lightblue3", lwd =2)
lines(and_CI[2,], col = "lightblue3", lwd=2)
legend("topright", col = c("darkgray", "lightblue3","lightpink"),
       legend = c("corrected", "Raw + Census","Raw"),lty=1, lwd=3, cex=0.7)
yr <- 2012:(years[y])
axis(1, at=agb_poly, labels=yr)
