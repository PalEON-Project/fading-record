# explore data product
file.in.model  = "HF_PREDS_tree_iter_25June2019.RDS" 

# this file has these data:
# tree year iter      dbh       ab     model plot taxon census_id dist_census


# read in PalEON tree-ring model for Harvard Forest
hf_rec <- readRDS(file.in.model)

# this function tries to determine how to classify a tree, big or small
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

hf_rec <- big_or_small(hf_rec)


inner_circle_area <- (pi*13^2) / 10000
outer_ring_area <- ((pi*20^2) / 10000) - inner_circle_area



########### explore 
hf_census <- hf_rec[hf_rec$model == "Model RW + Census", ]
prec <- crec <- outsmall <- list()
years <- 2012:1960
i <- y <- 1 
for( i in 1:250){
  prec[[i]] <- crec[[i]] <- matrix(NA, ncol= 3, nrow= length(years))
  
  for( y in seq_along(years)){
    
    # will pretend paleon sampling
    cur_p <- hf_census[hf_census$year == years[y] & hf_census$iter == i, ]
    
    # census
    cur <- hf_census[hf_census$year == years[y] & hf_census$iter == i, ]
    
    cur_p <- cur_p[complete.cases(cur_p),]
    cur   <- cur[complete.cases(cur),]
    
    # apply FIA cutoff
    cur <- cur[cur$dbh >= 5*2.54, ]
    cur_p <- cur_p[cur_p$dbh >= 5*2.54, ]
    
    # remove small trees in outside ring (paleon sampling)
    cur_p <- cur_p[!(cur_p$small == 1 & cur_p$dist_census > 13), ]
    
    iter_df <- cur_p
    small_trees <-  iter_df[(cur_p$dist_census < 13 & cur_p$small == 1),]
    #print(sum((cur_f$dist_census < 13 & cur_f$dbh < 20)))
    #small_trees <-  iter_df[iter_df$small == 1, ]
    
    # approximate how many small trees there should be in the outer ring
    x_trees <- ceiling(outer_ring_area* nrow(small_trees) /inner_circle_area)
    
    # create small trees
    x_ind <- sample(1:nrow(small_trees), x_trees, replace = TRUE)
    iter_df <- rbind(iter_df, small_trees[x_ind,])
    
    #cur_p$dens <- ifelse(cur_p$dbh <20, 1/(pi*13^2),1/(pi*20^2))*1e+04
    
    # we will compare these two values
    # we first removed, then created small trees in the outer ring
    pdbh <- mean(iter_df$dbh) 
    cdbh <- mean(cur$dbh) # didn't do anything, this is the real plot census without paleon sampling
    
    pab <- sum(iter_df$ab) * 1e-03 / ((pi*20^2) / 10000)
    cab <- sum(cur$ab)* 1e-03 / ((pi*20^2) / 10000)
    
    pden <- nrow(iter_df) / ((pi*20^2) / 10000) #sum(cur_p$dens)
    cden <- nrow(cur) / ((pi*20^2) / 10000)
    
    prec[[i]][y,] <- cbind(pdbh,pden,pab)
    crec[[i]][y,] <- cbind(cdbh,cden,cab)
  }
  print(i) # this is slow!
}


# COMPARE - p: micmicking paleon sampling, c: census
prec_mat <- do.call("rbind",prec)
crec_mat <- do.call("rbind",crec)
plot(prec_mat[,1],crec_mat[,1], main="dbh", xlab="paleon-like", ylab="census")
abline(0,1,col="red",lwd=3)

plot(prec_mat[,2],crec_mat[,2],main="dens", xlab="paleon-like", ylab="census")
abline(0,1,col="red",lwd=3)

plot(prec_mat[,3],crec_mat[,3],main="ab", xlab="paleon-like", ylab="census")
abline(0,1,col="red",lwd=3)

bias.fit = lm(crec_mat[,3] ~ prec_mat[,3])
abline(bias.fit,col=3,lty=2,lwd=2)

plot(frec_mat[,3],crec_mat[,3],main="ab")
abline(0,1,col="red",lwd=3)
