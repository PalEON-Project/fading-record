library(fields)
library(data.table)
library(PEcAn.allometry)

####################################### Read data in #######################################

file.in.normal = "../../data/ESM.rds"
file.in.fading = "../../data/ESM_F.rds"
file.in.model  = "../../data/HF_PREDS_tree_iter_25June2019_IF.RDS"
#file.in.model  = "../../data/HF_FR_input_corrected.RDS"


# --- read in FIA DB query
# read in ESM 
q = as.data.table(readRDS(file.in.normal))
# read in ESM fading record scenario 
w = as.data.table(readRDS(file.in.fading))

# q_in = q[which(q$COMPONENT == 'INGROWTH'),]
# tree_id = as.character(q_in[2, 'TRE_CN'])
# 
# q[which(q$TRE_CN == tree_id),]

sumNA  = function(x) sum(x,na.rm=T)
meanNA = function(x) mean(x,na.rm=T)

q = q[which(q$COMPONENT != 'INGROWTH'),]
w = w[which(w$COMPONENT != 'INGROWTH'),]

# convert diameter units before taking mean stand diameter (in case some trees get filtered)
q[, DIA_BEGIN     := DIA_BEGIN     * 2.54      ] # inches to cm
q[, DIA_END       := DIA_END       * 2.54      ]

w[, DIA_BEGIN     := DIA_BEGIN     * 2.54      ] # inches to cm
w[, DIA_END       := DIA_END       * 2.54      ]

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
# q[, diamean.m     := diamean     * 2.54      ] # inches to cm
# q[, prevdiamean.m := prevdiamean * 2.54      ]
q[, tpasum.m      := tpasum      / 0.404686  ] # trees acre-1 to trees ha-1 0.404686
q[, prevtpasum.m  := prevtpasum  / 0.404686  ]

# w[, diamean.m     := diamean     * 2.54      ]
# w[, prevdiamean.m := prevdiamean * 2.54      ]
w[, tpasum.m      := tpasum      / 0.404686  ]
w[, prevtpasum.m  := prevtpasum  / 0.404686  ]


fia_to_paleon   <- read.table("../../data/fia_to_level3a_v0.5.csv", heade = TRUE, sep="\t")
chojnacky_table <- read.table("../../data/level3a_to_chojnacky_v0.2.csv", header = TRUE, sep = "\t")

# there isn't any, but filter out plots that don't have FIA spcd
q <- q[!(q$plt_cn %in% unique(q$plt_cn[is.na(q$spcd)])),]
# filter out plots tht don't have matching spp
plots_with_these_spcd <- fia_to_paleon$fia_spcd[!(fia_to_paleon$level3a %in% chojnacky_table$level3a)]
q <- q[!(q$plt_cn %in% unique(q$plt_cn[q$spcd %in% plots_with_these_spcd])), ]

# match it to PalEON PFTs
length(unique(fia_to_paleon$fia_spcd))
length(unique(fia_to_paleon$level3a))
q$paleon <- plyr::mapvalues(q$spcd, fia_to_paleon$fia_spcd, as.character(fia_to_paleon$level3a))

mort_prop = sum(ifelse((q$prevtpasum - q$tpasum)>0, 1, 0))/nrow(q)


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
for( i in 1:50){#250){
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

# this function rescales the 5-yr ESM arrow and calculates AB
update_ab_current <- function(these_plots_cur, current_pnt, previous_pnt, esm, sfc){
  rescaled_mean_dbhs <- rep(NA, nrow(these_plots_cur))
  rescaled_mean_dens <- rep(NA, nrow(these_plots_cur))
  rescaled_mean_agbs <- rep(NA, nrow(these_plots_cur))
  
  #recalculate ab allometrically for rescaled
  for(p in seq_len(nrow(these_plots_cur))){
    
    sub_plot <- esm[esm$plt_cn == these_plots_cur$plt_cn[p], ] 
    
    # first scale esm to per year rate
    sf <- 1/(sub_plot$invyr[1] - sub_plot$prevyr[1])
    
    
    if ((sub_plot$tpasum.m[1] - sub_plot$prevtpasum.m[1])==0){
      rescaled_mean_dbhs[p] <- 0
      rescaled_mean_dens[p] <- 0
      rescaled_mean_agbs[p] <- 0
    } else {
      
      dbh_correction <- (sub_plot$diamean.m[1] - sub_plot$prevdiamean.m[1]) * sf
      den_correction <- (sub_plot$tpasum.m[1] - sub_plot$prevtpasum.m[1]) * sf
      
      beta0s <- chojnacky_table$beta0[match(sub_plot$paleon, as.character(chojnacky_table$level3a))]
      beta1s <- chojnacky_table$beta1[match(sub_plot$paleon, as.character(chojnacky_table$level3a))]
      
      # current plot
      ### exp(b0 + b1*log(dbh))
      agb_end           <- sum(exp(beta0s + beta1s*log(sub_plot$dia_end* 2.54)), na.rm = TRUE)
      agb_Mg_end        <- agb_end * 1e-03 # kg to Mg
      agb_Mg_plt_end    <- agb_Mg_end / (4*7.3152^2*pi) # Mg/m2
      agb_Mg_ha_plt_end <- agb_Mg_plt_end / 1e-04  # Mg/m2 to Mg/ha
      
      # previous plot
      ### exp(b0 + b1*log(dbh))
      agb           <- sum(exp(beta0s + beta1s*log(sub_plot$dia_begin* 2.54)), na.rm = TRUE)
      agb_Mg        <- agb * 1e-03 # kg to Mg
      agb_Mg_plt    <- agb_Mg / (4*7.3152^2*pi) # Mg/m2
      agb_Mg_ha_plt <- agb_Mg_plt / 1e-04  # Mg/m2 to Mg/ha
      
      
      agb_correction = (agb_Mg_ha_plt_end - agb_Mg_ha_plt)*sf
      
      # rescaled_mean_dbhs[p] <- mean(sub_plot$dia_begin * 2.54, na.rm=TRUE)
      # rescaled_mean_dens[p] <- sum(sub_plot$start1tpa / 0.404686, na.rm = TRUE)
      # rescaled_mean_agbs[p] <- agb_Mg_ha_plt
      rescaled_mean_dbhs[p] <- current_pnt[1] - dbh_correction
      rescaled_mean_dens[p] <- current_pnt[2] - den_correction
      rescaled_mean_agbs[p] <- current_pnt[3] - agb_correction
    }
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
  load("../../data/ab_lut.Rdata")
  
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
  
  search_esm2 <- function(ab_lut, ab_lut_norm, current_pnt, th, mins, maxs, nneigh){
    
    # normalize current position
    current_pnt_norm <- current_pnt
    current_pnt_norm[,1] <- (current_pnt_norm[,1] - mins[1])/ (maxs[1]-mins[1])
    current_pnt_norm[,2] <- (current_pnt_norm[,2] - mins[2]) / (maxs[2]-mins[2])
    current_pnt_norm[,3] <- (current_pnt_norm[,3] - mins[3]) / (maxs[3]-mins[3])
    
    # # get ring-width recon value from past year
    # previous_msd_norm <- (previous_pnt[,1]- mins[1]) / (maxs[1]-mins[1])
    # previous_ab_norm <- (previous_pnt[,3]- mins[3]) / (maxs[3]-mins[3])
    
    # find the fia plots with biomass similar to ring-width recon from past year
    sub_norm <- ab_lut_norm[abs(ab_lut[,4] - current_pnt[1,3]) < th,]
    # sub_norm <- ab_lut_norm[abs(ab_lut[,4] - previous_pnt[1,3]) < th,]
    
    # dold <- rdist(sub_norm[,c(2:4)], current_pnt_norm) # Euclidean distance in 3D
    
    # only use mean dbh and ab to find closest
    # d <- rdist(sub_norm[,c(2,4)], matrix(current_pnt_norm[c(1,3)],nrow=1)) # Euclidean distance in 3D
    
    d <- rdist(sub_norm[,c(2,4)], matrix(current_pnt_norm[c(1,3)],nrow=1)) # Euclidean distance in 3D
    
    these_plots <- sub_norm[(order(d)[1:(nneigh*1)]),] # use nneigh nearest-neighbours
    # these_plots <- these_plots[which(these_plots$ab >= previous_ab_norm),1]
    these_plots <- these_plots[,1]
    # little jitter
    #these_plots <- these_plots[sample(1:(nneigh*1), nneigh)]
    these_plots <- ab_lut[ab_lut$plt_cn %in% these_plots, ]
    
    return(these_plots)
  }
  
  
  search_esm2 <- function(ab_lut, ab_lut_norm, current_pnt, th, mins, maxs, nneigh){
    
    # normalize current position
    current_pnt_norm <- current_pnt
    current_pnt_norm[,1] <- (current_pnt_norm[,1] - mins[1])/ (maxs[1]-mins[1])
    current_pnt_norm[,2] <- (current_pnt_norm[,2] - mins[2]) / (maxs[2]-mins[2])
    current_pnt_norm[,3] <- (current_pnt_norm[,3] - mins[3]) / (maxs[3]-mins[3])
    
    # # get ring-width recon value from past year
    # previous_msd_norm <- (previous_pnt[,1]- mins[1]) / (maxs[1]-mins[1])
    # previous_ab_norm <- (previous_pnt[,3]- mins[3]) / (maxs[3]-mins[3])
    
    # find the fia plots with biomass similar to ring-width recon from past year
    sub_norm <- ab_lut_norm[abs(ab_lut[,4] - current_pnt[1,3]) < th,]
    # sub_norm <- ab_lut_norm[abs(ab_lut[,4] - previous_pnt[1,3]) < th,]
    
    # dold <- rdist(sub_norm[,c(2:4)], current_pnt_norm) # Euclidean distance in 3D
    
    # only use mean dbh and ab to find closest
    # d <- rdist(sub_norm[,c(2,4)], matrix(current_pnt_norm[c(1,3)],nrow=1)) # Euclidean distance in 3D
    
    d <- rdist(sub_norm[,c(2,4)], matrix(current_pnt_norm[c(1,3)],nrow=1)) # Euclidean distance in 3D
    
    these_plots <- sub_norm[(order(d)[1:(nneigh*1)]),] # use nneigh nearest-neighbours
    # these_plots <- these_plots[which(these_plots$ab >= previous_ab_norm),1]
    these_plots <- these_plots[,1]
    # little jitter
    #these_plots <- these_plots[sample(1:(nneigh*1), nneigh)]
    these_plots <- ab_lut[ab_lut$plt_cn %in% these_plots, ]
    
    return(these_plots)
  }
  
  ######################################################################################################################################
  ##
  ######################################################################################################################################
  get_correction <- function(these_plots, q, w){
    correction_dbhs <- rep(NA, nrow(these_plots))
    correction_dens <- rep(NA, nrow(these_plots))
    correction_agbs <- rep(NA, nrow(these_plots))
    
    for (plt in 1:nrow(these_plots)){
      q_sub <- q[which(q$plt_cn %in% these_plots$plt_cn[plt]),]
      w_sub <- w[which(w$plt_cn %in% these_plots$plt_cn[plt]),]
      
      sf <- 1/(q_sub$invyr[1] - q_sub$prevyr[1])
      
      
      if ((q_sub$tpasum[1] - q_sub$prevtpasum[1])==0){
        fr_dbh_correction = 0 
        fr_dbh_growth <-  (w_sub$diamean[1] - w_sub$prevdiamean[1]) * sf
        
        fr_den_correction = 0
        fr_den_growth  <- (w_sub$tpasum.m[1] - w_sub$prevtpasum.m[1]) * sf
        
        fr_agb_correction = 0
        
        w_beta0s <- chojnacky_table$beta0[match(q_sub$paleon, as.character(chojnacky_table$level3a))]
        w_beta1s <- chojnacky_table$beta1[match(q_sub$paleon, as.character(chojnacky_table$level3a))]
        
        # current plot
        ### exp(b0 + b1*log(dbh))
        w_agb_end           <- sum(exp(w_beta0s + w_beta1s*log(w_sub$dia_end)), na.rm = TRUE)
        w_agb_Mg_end        <- w_agb_end * 1e-03 # kg to Mg
        w_agb_Mg_plt_end    <- w_agb_Mg_end / (4*7.3152^2*pi) # Mg/m2
        w_agb_Mg_ha_plt_end <- w_agb_Mg_plt_end / 1e-04  # Mg/m2 to Mg/ha
        
        # previous plot
        ### exp(b0 + b1*log(dbh))
        w_agb           <- sum(exp(w_beta0s + w_beta1s*log(w_sub$dia_begin)), na.rm = TRUE)
        w_agb_Mg        <- w_agb * 1e-03 # kg to Mg
        w_agb_Mg_plt    <- w_agb_Mg / (4*7.3152^2*pi) # Mg/m2
        w_agb_Mg_ha_plt <- w_agb_Mg_plt / 1e-04  # Mg/m2 to Mg/ha
        
        
        fr_agb_growth = (w_agb_Mg_ha_plt_end - w_agb_Mg_ha_plt)*sf
      } else {  
        
        # q_dbh_correction <- (q_sub$diamean[1] - q_sub$prevdiamean[1]) * sf
        # w_dbh_correction <- (w_sub$diamean[1] - w_sub$prevdiamean[1]) * sf
        # 
        # # observe fading now, but in past want not faded
        # fr_dbh_correction <- (q_sub$diamean.m[1] - w_sub$prevdiamean.m[1]) * sf
        fr_dbh_correction <- (q_sub$prevdiamean[1] - w_sub$prevdiamean[1]) * sf
        fr_dbh_growth <-  (w_sub$diamean[1] - w_sub$prevdiamean[1]) * sf
        
        
        q_den_correction  <- (q_sub$tpasum.m[1] - q_sub$prevtpasum.m[1]) * sf
        fr_den_growth  <- (w_sub$tpasum.m[1] - w_sub$prevtpasum.m[1]) * sf
        fr_den_correction <- (q_sub$prevtpasum.m[1] - w_sub$prevtpasum.m[1]) * sf
        
        q_beta0s <- chojnacky_table$beta0[match(q_sub$paleon, as.character(chojnacky_table$level3a))]
        q_beta1s <- chojnacky_table$beta1[match(q_sub$paleon, as.character(chojnacky_table$level3a))]
        
        # current plot
        ### exp(b0 + b1*log(dbh))
        q_agb_end           <- sum(exp(q_beta0s + q_beta1s*log(q_sub$dia_end)), na.rm = TRUE)
        q_agb_Mg_end        <- q_agb_end * 1e-03 # kg to Mg
        q_agb_Mg_plt_end    <- q_agb_Mg_end / (4*7.3152^2*pi) # Mg/m2
        q_agb_Mg_ha_plt_end <- q_agb_Mg_plt_end / 1e-04  # Mg/m2 to Mg/ha
        
        # previous plot
        ### exp(b0 + b1*log(dbh))
        q_agb           <- sum(exp(q_beta0s + q_beta1s*log(q_sub$dia_begin)), na.rm = TRUE)
        q_agb_Mg        <- q_agb * 1e-03 # kg to Mg
        q_agb_Mg_plt    <- q_agb_Mg / (4*7.3152^2*pi) # Mg/m2
        q_agb_Mg_ha_plt <- q_agb_Mg_plt / 1e-04  # Mg/m2 to Mg/ha
        
        
        q_agb_correction = (q_agb_Mg_ha_plt_end - q_agb_Mg_ha_plt)*sf
        
        
        w_beta0s <- chojnacky_table$beta0[match(q_sub$paleon, as.character(chojnacky_table$level3a))]
        w_beta1s <- chojnacky_table$beta1[match(q_sub$paleon, as.character(chojnacky_table$level3a))]
        
        # current plot
        ### exp(b0 + b1*log(dbh))
        w_agb_end           <- sum(exp(w_beta0s + w_beta1s*log(w_sub$dia_end)), na.rm = TRUE)
        w_agb_Mg_end        <- w_agb_end * 1e-03 # kg to Mg
        w_agb_Mg_plt_end    <- w_agb_Mg_end / (4*7.3152^2*pi) # Mg/m2
        w_agb_Mg_ha_plt_end <- w_agb_Mg_plt_end / 1e-04  # Mg/m2 to Mg/ha
        
        # previous plot
        ### exp(b0 + b1*log(dbh))
        w_agb           <- sum(exp(w_beta0s + w_beta1s*log(w_sub$dia_begin)), na.rm = TRUE)
        w_agb_Mg        <- w_agb * 1e-03 # kg to Mg
        w_agb_Mg_plt    <- w_agb_Mg / (4*7.3152^2*pi) # Mg/m2
        w_agb_Mg_ha_plt <- w_agb_Mg_plt / 1e-04  # Mg/m2 to Mg/ha
        
        
        fr_agb_growth = (w_agb_Mg_ha_plt_end - w_agb_Mg_ha_plt)*sf
        
        fr_agb_correction = (q_agb_Mg_ha_plt - w_agb_Mg_ha_plt)*sf
      }
      
      correction_dbhs[plt] = fr_dbh_correction #- fr_dbh_growth
      correction_dens[plt] = fr_den_correction
      correction_agbs[plt] = fr_agb_correction #- fr_agb_growth
      
    }
    
    return(data.frame(dbh=correction_dbhs, dens=correction_dens, agb=correction_agbs))
  }
  
  apply_correction <- function(current_pnt, growth, correction){
    
    
    
    new_pnt = t(apply(correction, 1, function(x) x + current_pnt - growth))
    
    # new_pnt = current_pnt - growth + correction
    
    # new_pnt = current_pnt - growth + colMeans(correction) 
    
    return(new_pnt)
  }
  
  ################################################################################################
  
  years <- 2012:1962
  th <- 1
  nneigh <- 3 # was 20
  i <- y <- 1 
  corrected_ab <- list()
  
  pdf('../../figures/fr_test.pdf')
  # i can be 250
  for( i in 1:50){
    print(i)
    
    plot(12:40,seq(100,1200, length.out= length(12:40)),"n",xlab="dbh",ylab="dens",main=paste0("nneigh: ",nneigh, " i: ",i, " r: ", radd))
    # arrows(meanB$diameanB,meanB$tpasumB,meanB$prevdiameanB,meanB$prevtpasumB, col="lightgray",
    #        length=0.065, angle=22,lwd=1) 
    legend("topright", legend =c("paleon tree_ring only", "paleon tr+census", "esm correction"),
           col = c(2,1,4), lty=1)
    corrected_ab[[i]] <- list()
    #corrected_ab[[i]] <- matrix(NA, ncol = 3, nrow = length(years)+1)
    
    for( y in seq_along(years)){
      
      print(y)
      
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
        current_pnt  <- frec[[i]][y,, drop=FALSE]
        previous_pnt <- frec[[i]][y+1,, drop=FALSE]
        # these_plots <- search_esm(ab_lut, ab_lut_norm, current_pnt, th, mins, maxs, nneigh)
        these_plots <- search_esm2(ab_lut, ab_lut_norm, current_pnt, th, mins, maxs, nneigh)
        
        # # these_plots_past_old <- update_ab(these_plots, q, 1)
        # these_plots_past <- update_ab_current(these_plots, current_pnt, previous_pnt, q, 1)
        correction = get_correction(these_plots, q, w)
        these_plots = apply_correction(current_pnt, current_pnt - previous_pnt, correction)
        
        # these_plots =  these_plots[,2:4]
        colnames(these_plots) = colnames(correction)
        colnames(current_pnt) = colnames(correction)
        corrected_ab[[i]][[y]] <- current_pnt #these_plots # fix this???
        
        # corrected_ab[[i]][[y+1]]  <- these_plots
        corrected_ab[[i]][[y+1]]  <- these_plots[sample(seq(1,nrow(these_plots)), size=1, replace=FALSE),]
        
      } else {
        # use corrected arrow of the precious iteration
        # find the closest FIA plots to the beginning of the red arrow
        current_pnts <- corrected_ab[[i]][[y]]
        # previous_pnt <- frec[[i]][y+1,, drop=FALSE]
        #current_pnt <- corrected_ab[[i]][y,,drop=FALSE]
        
        current_data  <- frec[[i]][y,, drop=FALSE]
        previous_data <- frec[[i]][y+1,, drop=FALSE]
        
        # if there is only one current point
        if (nrow(current_pnts)==1){
          new_plots <- search_esm2(ab_lut, ab_lut_norm, current_pnts, th, mins, maxs, nneigh)
          # new_plots <- update_ab(new_plots, q, 1)
          # new_plots <- update_ab_current(new_plots, current_pnts, previous_pnt, q, 1)
          correction = get_correction(new_plots, q, w)
          new_plots = apply_correction(current_pnts, current_data-previous_data, correction)
          
          
          # if there are two or more current points
        } else if (nrow(current_pnts)>1){
          # new_plots <- lapply(nrow(current_pnts), function(x){
          #   tmp_plots <- search_esm2(ab_lut, ab_lut_norm, current_pnts[x,,drop=FALSE], th, mins, maxs, 
          #                            nneigh, previous_pnt)
          #   return(tmp_plots)
          # })
          new_plots = list()
          for (r in 1:nrow(current_pnts)){
            tmp_plots <- search_esm2(ab_lut, ab_lut_norm, current_pnts[r,,drop=FALSE], th, mins, maxs, 
                                     nneigh)
            correction = get_correction(tmp_plots, q, w)
            tmp_plots = apply_correction(current_pnts[r,,drop=FALSE],  current_data-previous_data, correction)
            new_plots[[r]] = tmp_plots
          }
          
          new_plots <- new_plots[which(sapply(new_plots, nrow)>0)]
          
          if (length(new_plots)==0){
            print("Cannot find an FIA plot near current value that satisfies assumptions.")
          }
          
          new_plots <- lapply(new_plots, function(x) x[sample(seq(1,nrow(x)), size=1, replace=FALSE),])
          # new_plots <- lapply(new_plots, function(x) apply(x, 2, median, na.rm=TRUE))
         # new_plots <- lapply(new_plots, colMeans, na.rm=TRUE)
          new_plots <- do.call("rbind",new_plots)
        } else {
          print("Missing current values for plot.")
        }
        
        corrected_ab[[i]][[y+1]] <- new_plots
        
      }
      
      # if (nrow(corrected_ab[[i]][[y]])==1){
        prev_arr = corrected_ab[[i]][[y]]
      # } else {
      #   prev_arr <- apply(corrected_ab[[i]][[y]], 2, mean, na.rm=TRUE)
      # }
      
      # if (nrow(corrected_ab[[i]][[y+1]])==1){
        next_arr = corrected_ab[[i]][[y+1]]
      # } else {
      #   next_arr <- apply(corrected_ab[[i]][[y+1]], 2, mean, na.rm=TRUE)
      # }
      
      # corrected arrow
      arrows(prev_arr[1], prev_arr[2], next_arr[1], next_arr[2], 
             col=4, length=0.065, angle=22)
      #arrows(corrected_ab[[i]][y,1], corrected_ab[[i]][y,2], 
      #       corrected_ab[[i]][y+1,1], corrected_ab[[i]][y+1,2], 
      #       col=4, length=0.065, angle=22)
    }
  }
  dev.off()
  
  
  
  i <- 5
  
  corrected_list = do.call(rbind, lapply(corrected_ab[[i]], function(x) colMeans(x)))
  plot(corrected_list[,3], crec[[i]][,3], xlim=c(140, 300), ylim=c(150, 300))
  abline(0,1)
  
  plot(years, corrected_list[1:length(years),3], type="l")
  lines(years, crec[[i]][1:length(years),3], col="blue")
  lines(years, frec[[i]][1:length(years),3], col="red")
  
  
  pdf('../../figures/fr_test_time.pdf')
  # corrected_ab[[i]] <- NULL
  y.list <- list()
  # for(l in 1:52){
    esm_correction <- list()
    for(tyi in seq_along(corrected_ab)){
      for (n in 1:nneigh){
        tmp <- sapply(corrected_ab[[tyi]], function(x) x[,3])
        tmp = t(do.call(rbind, tmp))
        #tmp <- corrected_ab[[tyi]][,3]
        # esm_correction = esm_correction
      }
      esm_correction[[tyi]] <- tmp#tmp[,l]
    }
    esm_correction = do.call(rbind, esm_correction)
    y.list <- apply(esm_correction, 2, function(x) quantile(x,  c(0.025, 0.5, 0.975)))
    #y.list[[l]] <- mean(unlist(esm_correction))
  # }
    
    
    
    
  ##################################################################################################
  agbCI_f <- y.list#do.call("cbind", y.list)
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
  dev.off()
  
  ######################################################################################################################################
  ##
  ######################################################################################################################################
  
  #load ab lut table
  load("../../data/ab_lut_yr.Rdata")
  
  #normalize
  ab_lut_norm <- ab_lut
  mins <- apply(ab_lut[,2:4], 2, min)
  maxs <- apply(ab_lut[,2:4], 2, max)
  ab_lut_norm[,2] <- (ab_lut_norm[,2] - mins[1])/ (maxs[1]-mins[1])
  ab_lut_norm[,3] <- (ab_lut_norm[,3] - mins[2]) / (maxs[2]-mins[2])
  ab_lut_norm[,4] <- (ab_lut_norm[,4] - mins[3]) / (maxs[3]-mins[3])
  
  
  # w fading
  # q not fading
  
  years <- 2012:1961
  th <- 1
  nneigh <- 10 # was 20
  i <- y <- 1 
  
  niter = length(frec)
  correct = list(niter)
  for (i in 1:niter){
    print(i)
    correct_these = list(length(years))
    for (y in 1:length(years)){
      current_pnt  <- frec[[i]][y,, drop=FALSE]
      # previous_pnt <- frec[[i]][y+1,, drop=FALSE]
      # these_plots <- search_esm(ab_lut, ab_lut_norm, current_pnt, th, mins, maxs, nneigh)
      these_plots <- search_esm2(ab_lut, ab_lut_norm, current_pnt, th, mins, maxs, nneigh)
      correct_these[[y]] <- get_correction(these_plots, q, w)
    }
    correct[[i]] = correct_these
  }
  ######################################################################################################################################
  ## METHOD 0: If add back, only add back if correction is positive, otherwise add back 0
  ######################################################################################################################################
  
  correct_cumsum <- matrix(NA, nrow=length(years), ncol=0)
  for (iter in 1:length(correct)){
    corrected_iter = matrix(NA, nrow=length(years), ncol=nneigh)
    for (y in 1:(length(years))){
      current_pnt = frec[[iter]][y,3]
      previous_pnt = frec[[iter]][y-1,3]
      
      if (y == 1){
        corrected_iter[y,1:nneigh] = rep(current_pnt, nneigh)
      } else {
        
        # the actual correction
        for (neigh in 1:nneigh){
          sum_correction = 0
          for (j in (y-1):1){
            add_success = 1#rbinom(1, 1, 0.3)
            if (add_success){
              add_back = ifelse(correct[[iter]][[j]][neigh,3]>=0, correct[[iter]][[j]][neigh,3], 0)
            } else {
              add_back = 0
            }
            
            sum_correction = sum_correction + add_back
            
          }
          
          corrected_iter[y,neigh] = current_pnt + sum_correction
          
        }
      }
      
    }
    correct_cumsum = cbind(correct_cumsum, corrected_iter)
  }
  
  corrected_ab = data.frame(year=years, correct_cumsum)
  corrected_ab_melt = melt(corrected_ab, id.vars=c('year'), factorsAsStrings = FALSE)
  
  agbCI_f <- apply(corrected_ab[,2:ncol(corrected_ab)], 1, function(x) quantile(x,  c(0.025, 0.5, 0.975), na.rm=TRUE))
  # agbCI_f <- do.call("cbind", y.list)
  #agbCI_f <- apply(esm_correction, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  agb_poly <- 1:dim(agbCI_f)[2]
  
  pal_tronly <-lapply(frec, function(x) x[,3])
  pal_tronly <- do.call("rbind", pal_tronly)
  and_CI_f <- apply(pal_tronly, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  
  pal_tr_cen <-lapply(crec, function(x) x[,3])
  pal_tr_cen <- do.call("rbind", pal_tr_cen)
  and_CI <- apply(pal_tr_cen, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  
  
  # pdf('../../figures/FADING-METHOD0-BIOMASS.pdf', width=12, height=10)
  png('../../figures/FADING-METHOD0-BIOMASS.png', width=620, height=480)
  plot(agbCI_f[1,], ylim = c(100,  320), xaxt = "n",lwd = 3, xlab = "Years", ylab = "AB (Mg/ha)",
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
  dev.off()
  
  
  ######################################################################################################################################
  ## METHOD 1: Draw from a bernoulli with 30% chance that we need to add something back
  ## If add back, only add back if correction is positive, otherwise add back 0
  ######################################################################################################################################
  
  correct_cumsum <- matrix(NA, nrow=length(years), ncol=0)
  for (iter in 1:length(correct)){
    corrected_iter = matrix(NA, nrow=length(years), ncol=nneigh)
    for (y in 1:(length(years))){
      current_pnt = frec[[iter]][y,3]
      previous_pnt = frec[[iter]][y-1,3]
      
      if (y == 1){
        corrected_iter[y,1:nneigh] = rep(current_pnt, nneigh)
      } else {
        
        # the actual correction
        for (neigh in 1:nneigh){
          sum_correction = 0
          for (j in (y-1):1){
            add_success = rbinom(1, 1, 0.3)
            if (add_success){
              add_back = ifelse(correct[[iter]][[j]][neigh,3]>=0, correct[[iter]][[j]][neigh,3], 0)
            } else {
              add_back = 0
            }
            
            sum_correction = sum_correction + add_back
            
          }
          
          corrected_iter[y,neigh] = current_pnt + sum_correction
          
        }
      }
      
    }
    correct_cumsum = cbind(correct_cumsum, corrected_iter)
  }
  
  corrected_ab = data.frame(year=years, correct_cumsum)
  corrected_ab_melt = melt(corrected_ab, id.vars=c('year'), factorsAsStrings = FALSE)
  
  agbCI_f <- apply(corrected_ab[,2:ncol(corrected_ab)], 1, function(x) quantile(x,  c(0.025, 0.5, 0.975), na.rm=TRUE))
  # agbCI_f <- do.call("cbind", y.list)
  #agbCI_f <- apply(esm_correction, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  agb_poly <- 1:dim(agbCI_f)[2]
  
  pal_tronly <-lapply(frec, function(x) x[,3])
  pal_tronly <- do.call("rbind", pal_tronly)
  and_CI_f <- apply(pal_tronly, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  
  pal_tr_cen <-lapply(crec, function(x) x[,3])
  pal_tr_cen <- do.call("rbind", pal_tr_cen)
  and_CI <- apply(pal_tr_cen, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  
  
  # pdf('../../figures/FADING-METHOD1-BIOMASS.pdf', width=12, height=10)
  png('../../figures/FADING-METHOD1-BIOMASS.png', width=620, height=480)
  plot(agbCI_f[1,], ylim = c(100,  320), xaxt = "n",lwd = 3, xlab = "Years", ylab = "AB (Mg/ha)",
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
  dev.off()
  
  ######################################################################################################################################
  ## METHOD 2: Draw from a bernoulli with 60% chance that we need to add something back
  ## If add back, and if correction is positive, add back something that is drawn from a distribution centered around correction
  ######################################################################################################################################
  
  # now try with many iters fixing sum_correction to be equal to or large than it was in previous year
  correct_cumsum <- matrix(NA, nrow=length(years), ncol=0)
  for (iter in 1:length(correct)){
    corrected_iter = matrix(NA, nrow=length(years), ncol=nneigh)
    sum_correction_iter = matrix(0, nrow=length(years), ncol=nneigh)
    for (y in 1:(length(years))){
      current_pnt = frec[[iter]][y,3]
      previous_pnt = frec[[iter]][y-1,3]
      
      if (y == 1){
        corrected_iter[y,1:nneigh] = rep(current_pnt, nneigh)
      } else {
        
        # the actual correction
        for (neigh in 1:nneigh){
          sum_correction = 0
          
          add_success = rbinom(1, 1, 0.6)
          if (add_success){
            # add_back = ifelse(correct[[iter]][[j]][neigh,3]>0, correct[[iter]][[j]][neigh,3], 0)
            add_back = ifelse(correct[[iter]][[j]][neigh,3]>0, runif(1, 1/2*correct[[iter]][[j]][neigh,3], 1.5*correct[[iter]][[j]][neigh,3]), 0)
            # add_back = ifelse(correct[[iter]][[j]][neigh,3]>0, rtruncnorm(1, 
            #                                                               a=0, 
            #                                                               b=2*correct[[iter]][[j]][neigh,3], 
            #                                                               mean=correct[[iter]][[j]][neigh,3],
            #                                                               sd=0.5), 0)
          } else {
            add_back = 0
          }
          
          sum_correction_iter[y,neigh] = sum_correction_iter[y-1,neigh] + add_back
          
          corrected_iter[y,neigh] = current_pnt + sum_correction_iter[y,neigh]
          
        }
      }
      
    }
    correct_cumsum = cbind(correct_cumsum, corrected_iter)
  }
  
  corrected_ab = data.frame(year=years, correct_cumsum)
  corrected_ab_melt = melt(corrected_ab, id.vars=c('year'), factorsAsStrings = FALSE)
  
  agbCI_f <- apply(corrected_ab[,2:ncol(corrected_ab)], 1, function(x) quantile(x,  c(0.025, 0.5, 0.975), na.rm=TRUE))
  # agbCI_f <- do.call("cbind", y.list)
  #agbCI_f <- apply(esm_correction, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  agb_poly <- 1:dim(agbCI_f)[2]
  
  pal_tronly <-lapply(frec, function(x) x[,3])
  pal_tronly <- do.call("rbind", pal_tronly)
  and_CI_f <- apply(pal_tronly, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  
  pal_tr_cen <-lapply(crec, function(x) x[,3])
  pal_tr_cen <- do.call("rbind", pal_tr_cen)
  and_CI <- apply(pal_tr_cen, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  
  
  # pdf('../../figures/FADING-METHOD2-BIOMASS.pdf', width=12, height=10)
  png('../../figures/FADING-METHOD2-BIOMASS.png', width=620, height=480)
  plot(agbCI_f[1,], ylim = c(100,  320), xaxt = "n",lwd = 3, xlab = "Years", ylab = "AB (Mg/ha)",
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
  dev.off()
  
  
  ######################################################################################################################################
  ## METHOD 3: If correction is positive, draw correction from uniform from 0 to correction
  ## 
  ######################################################################################################################################

  # now try with instead of binom, drawing from uniform from 0 to change
  correct_cumsum <- matrix(NA, nrow=length(years), ncol=0)
  for (iter in 1:length(correct)){
    corrected_iter = matrix(NA, nrow=length(years), ncol=nneigh)
    sum_correction_iter = matrix(0, nrow=length(years), ncol=nneigh)
    for (y in 1:(length(years))){
      current_pnt = frec[[iter]][y,3]
      previous_pnt = frec[[iter]][y-1,3]
      
      if (y == 1){
        corrected_iter[y,1:nneigh] = rep(current_pnt, nneigh)
      } else {
        
        # the actual correction
        for (neigh in 1:nneigh){
          sum_correction = 0
          
          # add_success = rbinom(1, 1, 0.5)
          # if (add_success){
          #   add_back = ifelse(correct[[iter]][[j]][neigh,3]>0, correct[[iter]][[j]][neigh,3], 0)
          # } else {
          #   add_back = 0
          # }
          if (correct[[iter]][[j]][neigh,3]>0){
           add_back = runif(1, min=0, max=correct[[iter]][[j]][neigh,3])
          } else {
            add_back=0
          }
          
          sum_correction_iter[y,neigh] = sum_correction_iter[y-1,neigh] + add_back
          
          corrected_iter[y,neigh] = current_pnt + sum_correction_iter[y,neigh]
          
        }
      }
      
    }
    correct_cumsum = cbind(correct_cumsum, corrected_iter)
  }
  
  corrected_ab = data.frame(year=years, correct_cumsum)
  corrected_ab_melt = melt(corrected_ab, id.vars=c('year'), factorsAsStrings = FALSE)
  
  agbCI_f <- apply(corrected_ab[,2:ncol(corrected_ab)], 1, function(x) quantile(x,  c(0.025, 0.5, 0.975), na.rm=TRUE))
  # agbCI_f <- do.call("cbind", y.list)
  #agbCI_f <- apply(esm_correction, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  agb_poly <- 1:dim(agbCI_f)[2]
  
  pal_tronly <-lapply(frec, function(x) x[,3])
  pal_tronly <- do.call("rbind", pal_tronly)
  and_CI_f <- apply(pal_tronly, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  
  pal_tr_cen <-lapply(crec, function(x) x[,3])
  pal_tr_cen <- do.call("rbind", pal_tr_cen)
  and_CI <- apply(pal_tr_cen, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  
  
  # pdf('../../figures/FADING-METHOD3-BIOMASS.pdf', width=12, height=10)
  png('../../figures/FADING-METHOD3-BIOMASS.png', width=620, height=480)
  plot(agbCI_f[1,], ylim = c(100,  320), xaxt = "n",lwd = 3, xlab = "Years", ylab = "AB (Mg/ha)",
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
  dev.off()
  
  ######################################################################################################################################
  ## apply cumulative sum correction
  ######################################################################################################################################
  
  
  # # works
  # correct_cumsum <- matrix(NA, nrow=length(years), ncol=0)
  # for (iter in 1:length(correct)){
  #   corrected_iter = matrix(NA, nrow=length(years), ncol=nneigh)
  #   for (i in 1:(length(years))){
  #     current_pnt = frec[[iter]][i,3]
  #     previous_pnt = frec[[iter]][i-1,3]
  # 
  #     if (i == 1){
  #       corrected_iter[i,1:nneigh] = rep(current_pnt, nneigh)
  #     } else {
  #       if (sign(correct[[iter]][[i-1]][1,3])==1){
  #         sum_correction = 0
  #         for (j in (i-1):1){
  #           if (sign(correct[[iter]][[j]][,3]) ==1){
  #             sum_correction = sum_correction + correct[[iter]][[j]][,3]
  #           }
  #         }
  #         corrected_iter[i,1:nneigh] = current_pnt + sum_correction
  #       } else {
  #         corrected_iter[i,1:nneigh] = NA
  #       }
  #     }
  # 
  #   }
  #   correct_cumsum = cbind(correct_cumsum, corrected_iter)
  # }

  library(reshape2)
  corrected_ab = data.frame(year=years, correct_cumsum)


  corrected_ab_melt = melt(corrected_ab, id.vars=c('year'), factorsAsStrings = FALSE)
  
  
  
  

  agbCI_f <- apply(corrected_ab[,2:ncol(corrected_ab)], 1, function(x) quantile(x,  c(0.025, 0.5, 0.975), na.rm=TRUE))
  # agbCI_f <- do.call("cbind", y.list)
  #agbCI_f <- apply(esm_correction, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  agb_poly <- 1:dim(agbCI_f)[2]

  pal_tronly <-lapply(frec, function(x) x[,3])
  pal_tronly <- do.call("rbind", pal_tronly)
  and_CI_f <- apply(pal_tronly, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)

  pal_tr_cen <-lapply(crec, function(x) x[,3])
  pal_tr_cen <- do.call("rbind", pal_tr_cen)
  and_CI <- apply(pal_tr_cen, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)


  # pdf('../../figures/correct_biomass_cumsum.pdf')
  plot(agbCI_f[1,], ylim = c(100,  320), xaxt = "n",lwd = 3, xlab = "Years", ylab = "AB (Mg/ha)",
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
  # dev.off()

  
  
  ## not sure if this works
  corrected_ab <- matrix(NA, nrow=length(years), ncol=0)
  for (iter in 1:length(correct)){
    corrected_iter = matrix(NA, nrow=length(years), ncol=nneigh)
    for (i in 1:(length(years))){
      current_pnt = frec[[iter]][i,3]
      previous_pnt = frec[[iter]][i-1,3]
      
      if (i == 1){
        corrected_iter[i,1:nneigh] = rep(current_pnt, nneigh)
      } else {
        if(sign(correct[[iter]][[i-1]][1,3])==1){
          corrected_iter[i,1:nneigh] = current_pnt + correct[[iter]][[i-1]][,3]
        } else {
          corrected_iter[i,1:nneigh] = NA
        }
      }
      
    }
    corrected_ab = cbind(corrected_ab, corrected_iter)
  }
  
  corrected_ab = data.frame(year=years, corrected_ab)
  
  
  library(reshape2)
  corrected_ab_melt = melt(corrected_ab, id.vars=c('year'), factorsAsStrings = FALSE)
  
  agbCI_f <- apply(corrected_ab[,2:ncol(corrected_ab)], 1, function(x) quantile(x,  c(0.025, 0.5, 0.975), na.rm=TRUE))
  # agbCI_f <- do.call("cbind", y.list)
  #agbCI_f <- apply(esm_correction, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  agb_poly <- 1:dim(agbCI_f)[2]
  
  pal_tronly <-lapply(frec, function(x) x[,3])
  pal_tronly <- do.call("rbind", pal_tronly)
  and_CI_f <- apply(pal_tronly, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  
  pal_tr_cen <-lapply(crec, function(x) x[,3])
  pal_tr_cen <- do.call("rbind", pal_tr_cen)
  and_CI <- apply(pal_tr_cen, 2, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
  
  
  pdf('../../figures/correct_biomass.pdf')
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
  dev.off()
  
  
  # 
  # 
  # # this function rescales the 5-yr ESM arrow and calculates AB
  # get_correction <- function(these_plots_cur, current_pnt, previous_pnt, esm, sfc){
  #   rescaled_mean_dbhs <- rep(NA, nrow(these_plots_cur))
  #   rescaled_mean_dens <- rep(NA, nrow(these_plots_cur))
  #   rescaled_mean_agbs <- rep(NA, nrow(these_plots_cur))
  #   
  #   #recalculate ab allometrically for rescaled
  #   for(p in seq_len(nrow(these_plots_cur))){
  #     
  #     sub_plot <- esm[esm$plt_cn == these_plots_cur$plt_cn[p], ] 
  #     
  #     # first scale esm to per year rate
  #     sf <- 1/(sub_plot$invyr[1] - sub_plot$prevyr[1])
  #     # new_dbh <- sub_plot$diamean.m[1] + (sub_plot$prevdiamean.m[1] - sub_plot$diamean.m[1]) * sf
  #     # prcnt_dbh <- ((new_dbh- sub_plot$prevdiamean.m[1])/ sub_plot$prevdiamean.m[1])*100
  #     # sub_plot$dia_begin <- sub_plot$dia_begin + (sub_plot$dia_begin * prcnt_dbh/100)
  #     dbh_correction <- (sub_plot$diamean.m[1] - sub_plot$prevdiamean.m[1]) * sf
  #     
  #     # new_den <- sub_plot$tpasum.m[1] + (sub_plot$prevtpasum.m[1] - sub_plot$tpasum.m[1]) * sf
  #     # prcnt_den <- (new_den- sub_plot$prevtpasum.m[1])/ sub_plot$prevtpasum.m[1]*100
  #     # sub_plot$start1tpa <- sub_plot$start1tpa + (sub_plot$start1tpa * prcnt_den/100)
  #     den_correction <- (sub_plot$tpasum.m[1] - sub_plot$prevtpasum.m[1]) * sf
  #     
  #     # if we want to update the rescaling the arrow with site specific magnitudes (using sfc) do it here
  #     
  #     #new_dbh <- sub_plot$diamean[1] + (mean(sub_plot$dia_begin,na.rm=TRUE)-sub_plot$diamean[1])*sfc
  #     #new_den <- sub_plot$tpasum[1] + (sum(sub_plot$start1tpa,na.rm=TRUE)-sub_plot$tpasum[1])*sfc
  #     
  #     #prcnt_dbh <- (new_dbh- mean(sub_plot$dia_begin,na.rm=TRUE))/ mean(sub_plot$dia_begin,na.rm=TRUE)*100
  #     #prcnt_den <- (new_den- sum(sub_plot$start1tpa,na.rm=TRUE))/ sum(sub_plot$start1tpa,na.rm=TRUE)*100
  #     
  #     #sub_plot$dia_begin <- sub_plot$dia_begin + (sub_plot$dia_begin * prcnt_dbh/100)
  #     #sub_plot$start1tpa <- sub_plot$start1tpa + (sub_plot$start1tpa * prcnt_den/100)
  #     
  #     beta0s <- chojnacky_table$beta0[match(sub_plot$paleon, as.character(chojnacky_table$level3a))]
  #     beta1s <- chojnacky_table$beta1[match(sub_plot$paleon, as.character(chojnacky_table$level3a))]
  #     
  #     # current plot
  #     ### exp(b0 + b1*log(dbh))
  #     agb_end           <- sum(exp(beta0s + beta1s*log(sub_plot$dia_end* 2.54)), na.rm = TRUE)
  #     agb_Mg_end        <- agb_end * 1e-03 # kg to Mg
  #     agb_Mg_plt_end    <- agb_Mg_end / (4*7.3152^2*pi) # Mg/m2
  #     agb_Mg_ha_plt_end <- agb_Mg_plt_end / 1e-04  # Mg/m2 to Mg/ha
  #     
  #     # previous plot
  #     ### exp(b0 + b1*log(dbh))
  #     agb           <- sum(exp(beta0s + beta1s*log(sub_plot$dia_begin* 2.54)), na.rm = TRUE)
  #     agb_Mg        <- agb * 1e-03 # kg to Mg
  #     agb_Mg_plt    <- agb_Mg / (4*7.3152^2*pi) # Mg/m2
  #     agb_Mg_ha_plt <- agb_Mg_plt / 1e-04  # Mg/m2 to Mg/ha
  #     
  #     
  #     agb_correction = (agb_Mg_ha_plt_end - agb_Mg_ha_plt)*sf
  #     
  #     # rescaled_mean_dbhs[p] <- mean(sub_plot$dia_begin * 2.54, na.rm=TRUE)
  #     # rescaled_mean_dens[p] <- sum(sub_plot$start1tpa / 0.404686, na.rm = TRUE)
  #     # rescaled_mean_agbs[p] <- agb_Mg_ha_plt
  #     rescaled_mean_dbhs[p] <- current_pnt[1] - dbh_correction
  #     rescaled_mean_dens[p] <- current_pnt[2] - den_correction
  #     rescaled_mean_agbs[p] <- current_pnt[3] - agb_correction
  #   }
  #   
  #   #corrected_mat <- matrix(apply(cbind(rescaled_mean_dbhs, rescaled_mean_dens, rescaled_mean_agbs), 2, mean), ncol=3)
  #   corrected_mat <- cbind(rescaled_mean_dbhs, rescaled_mean_dens, rescaled_mean_agbs)
  #   
  #   return(corrected_mat)
  # }