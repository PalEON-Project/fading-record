rm(list=ls()) ##
library(data.table)
library(RPostgreSQL)

# This script contains functions written by Travis Andrews and Ryan Kelly that query FIA DB to create ESM phase space
#
#  First 200 lines are PSQL helper functions (mainly from PEcAn), if you have PEcAn.DB package in your library you can just load it library(PEcAn.DB)
#  Then, fia_extract function returns the resulting query
#  Here we call fia_xtract function twice:
#      1) to create ESM phase space
#      2) to create ESM phase space that corresponds to the PalEON fading record scenario (i.e. no mortality, no ingrowth)


#####################################   PSQL helper functions ################################
# These are PEcAn DB functions, copied for convenience to this project
#-------------------------------------------------------------------------------
# Copyright (c) 2012 University of Illinois, NCSA.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the 
# University of Illinois/NCSA Open Source License
# which accompanies this distribution, and is available at
# http://opensource.ncsa.illinois.edu/license.html
#-------------------------------------------------------------------------------

.db.utils <- new.env() 
.db.utils$created <- 0
.db.utils$queries <- 0
.db.utils$deprecated <- 0
.db.utils$showquery <- FALSE
.db.utils$connections <- list()

#---------------- Base database query function. ---------------------------------------------------#
##' Generic function to query database
##'
##' Given a connection and a query, will return a query as a data frame. Either con or params need
##' to be specified. If both are specified it will use con.
db.query <- function(query, con=NULL, params=NULL) {
  iopened <- 0
  if(is.null(con)){
    if (is.null(params)) {
      logger.error("No parameters or connection specified")
      stop()
    }
    con <- db.open(params)
    iopened <- 1
  }
  if (.db.utils$showquery) {
    logger.debug(query)
  }
  data <- dbGetQuery(con, query)
  res <- dbGetException(con)
  if (res$errorNum != 0 || (res$errorMsg != 'OK' && res$errorMsg != '')) {
    logger.severe(paste("Error executing db query '", query, "' errorcode=", res$errorNum, " message='", res$errorMsg, "'", sep=''))
  }
  .db.utils$queries <- .db.utils$queries+1
  if(iopened==1) {
    db.close(con)
  }
  invisible(data)
}

##' Generic function to open a database connection
##'
##' Create a connection to a database usign the specified parameters. If the paramters contain
##' driver element it will be used as the database driver, otherwise it will use PostgreSQL.
db.open <- function(params) {
  params$dbfiles <- NULL
  params$write <- NULL
  if (is.null(params$driver)) {
    args <- c(drv=dbDriver("PostgreSQL"), params, recursive=TRUE)
  } else {
    args <- c(drv=dbDriver(params$driver), params, recursive=TRUE)
    args[['driver']] <- NULL
  }
  c <- do.call(dbConnect, as.list(args))
  id <- sample(1000, size=1)
  while(length(which(.db.utils$connections$id==id)) != 0) {
    id <- sample(1000, size=1)
  }
  attr(c, "pecanid") <- id
  dump.log <- NULL
  dump.frames(dumpto="dump.log")
  .db.utils$created <- .db.utils$created+1
  .db.utils$connections$id <- append(.db.utils$connections$id, id)
  .db.utils$connections$con <- append(.db.utils$connections$con, c)
  .db.utils$connections$log <- append(.db.utils$connections$log, list(dump.log))
  invisible(c)
}

##' Generic function to close a database connection
##'
##' Close a previously opened connection to a database.
db.close <- function(con) {
  if (is.null(con)) {
    return
  }
  
  id <- attr(con, "pecanid")
  if (is.null(id)) {
    logger.warn("Connection created outside of PEcAn.db package")
  } else {
    deleteme <- which(.db.utils$connections$id==id)
    if (length(deleteme) == 0) {
      logger.warn("Connection might have been closed already.");
    } else {
      .db.utils$connections$id <- .db.utils$connections$id[-deleteme]
      .db.utils$connections$con <- .db.utils$connections$con[-deleteme]
      .db.utils$connections$log <- .db.utils$connections$log[-deleteme]
    }
  }
  dbDisconnect(con)
}

##' Debug method for db.open and db.close
##'
##' Prints the number of connections opened as well as any connections
##' that have never been closes.
db.print.connections <- function() {
  logger.info("Created", .db.utils$created, "connections and executed", .db.utils$queries, "queries")
  if (.db.utils$deprecated > 0) {
    logger.info("Used", .db.utils$deprecated, "calls to deprecated functions")
  }
  logger.info("Created", .db.utils$created, "connections and executed", .db.utils$queries, "queries")
  if (length(.db.utils$connections$id) == 0) {
    logger.debug("No open database connections.\n")
  } else {
    for(x in 1:length(.db.utils$connections$id)) {
      logger.info(paste("Connection", x, "with id", .db.utils$connections$id[[x]], "was created at:\n"))
      logger.info(paste("\t", names(.db.utils$connections$log[[x]]), "\n"))
      #      cat("\t database object : ")
      #      print(.db.utils$connections$con[[x]])
    }
  }
}

##' Test connection to database
##' 
##' Useful to only run tests that depend on database when a connection exists
db.exists <- function(params, write=TRUE) {
  # open connection
  con <- tryCatch({
    invisible(db.open(params))
  }, error = function(e) {
    logger.error("Could not connect to database.\n\t", e)
    invisible(NULL)
  })
  if (is.null(con)) {
    return(invisible(FALSE))
  }
  
  if(FALSE){
    # read a row from the database
    read.result <- tryCatch({
      invisible(db.query("SELECT * FROM users LIMIT 1", con))
    }, error = function(e) {
      logger.error("Could not query database.\n\t", e)
      db.close(con)
      invisible(NULL)
    })
    if (is.null(read.result)) {
      return(invisible(FALSE))
    }
    
    # if requested write a row to the database
    if (write) {
      result <- tryCatch({
        db.query(paste("UPDATE users SET created_at='", read.result$created_at, "' WHERE id=", read.result$id, sep=""), con)
        invisible(TRUE)
      }, error = function(e) {
        logger.error("Could not write to database.\n\t", e)
        invisible(FALSE)
      })
    } else {
      result <- TRUE
    }
    
  }
  result <- TRUE
  
  # close database, all done
  tryCatch({
    db.close(con)
  }, error = function(e) {
    logger.warn("Could not close database.\n\t", e)
  })
  
  invisible(result)
}

##' Sets if the queries should be shown that are being executed
db.showQueries <- function(show) {
  .db.utils$showquery <- show
}

##' Returns if the queries should be shown that are being executed
db.getShowQueries <- function() {
  invisible(.db.utils$showquery)
}





#################################### Query FIA DB for ESM data ####################################

##' Function to query ESM data from FIA database
##' 
##' 
##' @param file.pft lookup table to match FIA species codes to temperate PFTs (e.g. Early Hardwood - EH)
##' @param file.out full path to resulting file
##' @param fading_record set to TRUE to mimic PalEON tree-ring fading record scenario (no mortality, no ingrowth), FALSE by default
##' @author Ryan Kelly, Travis Andrews, Istem Fer



file.pft = "../../data/spcd2pft.csv"
file.normal = "../../data/ESM.rds"
file.fading = "../../data/ESM_F.rds"

fia_extract(file.normal, file.pft)
fia_extract(file.fading, file.pft, fading_record = TRUE)

fia_extract <- function(file.out, file.pft, fading_record = FALSE){
  
  
  sumNA  = function(x) sum(x,na.rm=T)
  meanNA = function(x) mean(x,na.rm=T)
  maxNA  = function(x) max(x,na.rm=T)
  
  dbsettings = list(
    user     = "bety",            # PSQL username
    password = "bety",                # PSQL password
    dbname   = "fia5",  # PSQL database name
    host     = "psql-pecan.bu.edu",       # PSQL server address (don't change unless server is remote)
    driver   = 'PostgreSQL',      # DB driver (shouldn't need to change)
    write    = FALSE              # Whether to open connection with write access. 
  )
  
  lon.bounds = c(-95,999)
  lat.bounds = c(-999,999)
  
  # -----------------------------
  # Open connection to database
  fia.con = db.open(dbsettings)
  
  # ---------- PLOT & COND DATA
  # --- Query PLOT
  cat("Query PLOT...\n")
  query = paste('SELECT 
                cn, statecd, prev_plt_cn, remper
                FROM plot WHERE remper>3 AND remper<9.5 AND designcd=1 AND statecd<=56 AND ',
                'lon>', min(lon.bounds),' AND lon<', max(lon.bounds), ' AND ',
                'lat>', min(lat.bounds),' AND lat<', max(lat.bounds))
  
  # ~10 sec
  PLOT = as.data.table(db.query(query, con=fia.con))
  setnames(PLOT, toupper(names(PLOT)))
  setnames(PLOT,"CN","PLT_CN")
  
  # Remove states that haven't been resurveyed
  PLOT[, REMPERTOT := sumNA(REMPER), by=STATECD]
  PLOT = PLOT[ REMPERTOT>10, ]
  
  # Remove this one miscellaneous plot, per TA
  PLOT = PLOT[ PLT_CN!= 134680578010854 ]
  
  # Store statecd and state names
  states = sort(unique(PLOT$STATECD))
  n.state = length(states)
  surv = db.query('SELECT statecd,statenm FROM survey', con=fia.con)
  state.names = surv$statenm[match(states,surv$statecd)]
  
  
  # --- Query COND
  cat("Query COND...\n")
  query = paste('SELECT 
                plt_cn, condid, stdorgcd
                FROM cond WHERE 
                stdorgcd=0 AND ',
                'statecd IN (', paste(states,collapse=','), ')')
  
  # ~ 15 sec
  COND = as.data.table(db.query(query, con=fia.con))
  setnames(COND, toupper(names(COND)))
  
  # Remove all plots with more than 1 condition
  COND[, CONmax := maxNA(CONDID), by=PLT_CN]
  # *** RK: This is slightly wrong. In a few cases plots have CONDID>1, but still only have a single condition. This would work better:
  #     COND[, CONmax2 := .N, by=PLT_CN]
  COND = COND[ CONmax==1,]
  
  # --- Merge PLOT and COND
  cat("Merge PLOT and COND ...\n")
  PC = merge(COND, PLOT, by="PLT_CN")
  
  
  # ---------- RESURVEY DATA
  # --- Query
  cat("Query TREE_GRM_ESTN...\n")
  query = paste('SELECT 
                plt_cn, subtyp_grm, invyr, tpagrow_unadj, dia_begin, dia_end, component, tre_cn, remper, statecd
                FROM tree_grm_estn WHERE ',
                #                 'dia_begin>5 AND ',
                'statecd IN (', paste(states,collapse=','),') AND ',
                'estn_type=\'AL\' AND land_basis=\'TIMBERLAND\'')
  
  GRM = as.data.table(db.query(query, con=fia.con))
  setnames(GRM, toupper(names(GRM)))
  
  # --- Filtering
  cat("Filtering TREE_GRM_ESTN...\n")
  
  # By plot/cond criteria
  GRM = GRM[ PLT_CN %in% PC$PLT_CN ]
  
  # IF: make sure we don't have the microplots
  GRM = GRM[ SUBTYP_GRM == 1 ]
  
  # Assign GRM$START + GRM$CUT and restrict to cut==0, start>0
  GRM[, START      := INVYR - REMPER                                  ]
  GRM[, REMPER := NULL]
  GRM[, CUT1TPA    := (COMPONENT=="CUT1") * TPAGROW_UNADJ             ]
  GRM[, CUT2TPA    := (COMPONENT=="CUT2") * TPAGROW_UNADJ             ]
  GRM[, CUT        := sumNA(CUT2TPA + CUT1TPA), by=PLT_CN             ]
  GRM = GRM[ START>0 & CUT==0, ]
  
  # Assign Reversion/Diversion, and exclude plots with either
  GRM[, DIVERSION1TPA  := (COMPONENT=="DIVERSION1") * TPAGROW_UNADJ   ]
  GRM[, DIVERSION2TPA  := (COMPONENT=="DIVERSION2") * TPAGROW_UNADJ   ]
  GRM[, REVERSION1TPA  := (COMPONENT=="REVERSION1") * TPAGROW_UNADJ   ]
  GRM[, REVERSION2TPA  := (COMPONENT=="REVERSION2") * TPAGROW_UNADJ   ]
  GRM[, REDIV          := sumNA(REVERSION2TPA+REVERSION1TPA+DIVERSION2TPA+DIVERSION1TPA), by=PLT_CN]
  GRM = GRM[ REDIV==0, ] 
  
  # Assign SURVIVORTPA, and remove records from any state with <1000 measured trees
  GRM[, SURVIVORTPA    := (COMPONENT=="SURVIVOR") * TPAGROW_UNADJ     ]
  GRM[, TPATOT         := sumNA(SURVIVORTPA), by=STATECD              ]
  GRM = GRM[ TPATOT>1000, ]
  
  
  # --- Assign additional variables
  cat("Calculating TPA and Diameter...\n")
  # Compute TPA
  GRM[, INGROWTHTPA    := (COMPONENT=="INGROWTH") * TPAGROW_UNADJ     ]
  GRM[, MORTALITY1TPA  := (COMPONENT=="MORTALITY1") * TPAGROW_UNADJ   ]
  GRM[, MORTALITY2TPA  := (COMPONENT=="MORTALITY2") * TPAGROW_UNADJ   ]
  GRM[, MORTALITYTPA   := MORTALITY1TPA + MORTALITY2TPA               ]
  
  if(!fading_record){
    
    # Initial number of trees is current survivors plus those that died during the resurvey period.
    GRM[, start1tpa      := SURVIVORTPA + MORTALITYTPA                  ]
    
    # Final number of trees is current survivors plus new trees that cross the 5" threshold
    GRM[, end1tpa        := SURVIVORTPA + INGROWTHTPA                   ]
    
  }else{
    
    # no mortality
    GRM[, start1tpa      := SURVIVORTPA                   ]
    GRM[COMPONENT!="SURVIVOR", DIA_BEGIN := NA] # Set non-survivor values to NA
    
    # no ingrowth
    GRM[, end1tpa        := SURVIVORTPA + INGROWTHTPA                   ]
  }
  
  
  GRM[, PREVTPAsum     := sumNA(start1tpa), by=PLT_CN                 ]  # "startsumTPA"
  GRM[, TPAsum         := sumNA(end1tpa), by=PLT_CN                   ]  # "endsumTPA"
  
  # Compute plot mean diameters
  GRM[, PREVDIAmean    := meanNA(DIA_BEGIN), by=PLT_CN                ]  # "DIAbeginmean"
  GRM[, DIAmean        := meanNA(DIA_END),   by=PLT_CN                ]  # "DIAendmean"
  
  
  # --- Subset for output
  GRM.out = GRM[, .(PLT_CN, TRE_CN, PREVTPAsum, TPAsum, PREVDIAmean, DIAmean)]
  
  
  # ---------- TREE
  cat("Query TREE...\n")
  # --- Query
  query = paste('SELECT 
                cn, prev_tre_cn, plt_cn, invyr, condid, dia, tpa_unadj, spcd, stocking, statuscd, 
                prevdia, prev_status_cd, p2a_grm_flg, reconcilecd
                FROM tree WHERE 
                (prevdia>5 OR dia>5) AND (statuscd=1 OR prev_status_cd=1) AND p2a_grm_flg!=\'N\'  AND
                statecd IN (', paste(states,collapse=','), ')')
  
  # ~ 5 min
  TREE = as.data.table(db.query(query, con=fia.con))
  setnames(TREE, toupper(names(TREE)))
  
  # --- Filter TREE
  cat("Filter TREE ...\n")
  # By plot/cond criteria
  TREE = TREE[ PLT_CN %in% PC$PLT_CN ]
  
  # CONDID ("Remove edge effects" --TA)
  TREE[, CONmax := maxNA(CONDID), by=PLT_CN]
  
  # STATUSCD
  # *** RK: Next line looks wrong. It's a sum, not max, despite the name. I did rewrite the line but this is equivalent to what Travis had so keeping for now.
  TREE[, STATUSCDmax := sumNA(3*as.integer(STATUSCD==3)), by=PLT_CN]
  
  # RECONCILECD
  TREE[is.na(RECONCILECD), RECONCILECD :=0] # Set NA values to 0 (unused)
  
  # Filter
  TREE = TREE[ CONmax==1 & INVYR<2014 & STATUSCDmax!=3 & STATUSCD!=0 & RECONCILECD<=4 ]
  
  
  # --- Merge in PFTS
  cat("Merge in PFTs...\n")
  # ~ 30 sec
  MCDPFT = as.data.table(read.csv(file.pft, header = TRUE))
  TREE = merge(TREE, MCDPFT, all.x=T, by = "SPCD")
  
  
  # --- Connect PREV_CN for each tree prior to subset
  cat("Connect consecutive observations...\n")
  TREE.prev = TREE[,.(CN, STOCKING, SPCD, TPA_UNADJ, PFT)]
  setnames(TREE.prev, paste0("PREV_TRE_",names(TREE.prev))) 
  
  # # Convert PREV_TRE_CN columns to integer64 (have experienced crashes otherwise. memory leak?)
  # TREE.prev[, PREV_TRE_CN := as.integer64(PREV_TRE_CN)]
  # TREE[, PREV_TRE_CN := as.integer64(PREV_TRE_CN)]
  
  TREE = merge(TREE, TREE.prev, all.x=T, by="PREV_TRE_CN")
  setnames(TREE,"CN","TRE_CN")
  
  
  # --- Define DIA and STOCKING columns for trees >5"
  cat("Calculate DIA and STOCKING...\n")
  # DIAmean of DIA>5
  TREE[DIA>=5 & STATUSCD==1,                                   DIA5alive := DIA      ]
  TREE[, DIA5meanalive     := meanNA(DIA5alive), by=PLT_CN                         ]
  TREE[PREVDIA>=5 & PREV_STATUS_CD==1,                     PREVDIA5alive := PREVDIA  ]
  TREE[, PREVDIA5meanalive := meanNA(PREVDIA5alive), by=PLT_CN                     ]
  
  #Stocking of plots for trees with DIA>5
  TREE[DIA5alive>0, STOCKING5 := STOCKING]
  TREE[, STOCKING5mid := sumNA(STOCKING5), by=PLT_CN]
  TREE[PREVDIA5alive>0, PREVSTOCKING5 := PREV_TRE_STOCKING]
  TREE[, PREVSTOCKING5mid := sumNA(PREVSTOCKING5), by=PLT_CN]
  
  
  
  # ---------- MERGE
  cat("Final merge...\n")
  ALL = merge(GRM, TREE, all.x=T, by='TRE_CN')
  ALL[, c("PLT_CN.x","INVYR.x") := list(NULL,NULL)]
  setnames(ALL, c("PLT_CN.y","INVYR.y"), c("PLT_CN","INVYR"))
  ALL = merge(ALL, PC, by='PLT_CN')
  ALL[, c("STATECD.x","CONmax.x") := list(NULL,NULL)]
  setnames(ALL, c("STATECD.y","CONmax.y"), c("STATECD","CONmax"))
  setnames(ALL, "START", "PREVYR")
  
  # --- Save outputs
  cat("Save...\n")
  saveRDS(ALL, file = file.out)
  
  db.close(fia.con)
  
} # end of function fia_extract




# ---------- DIAGNOSTICS by RK
# q = merge(GRM,TREE, by='TRE_CN', select=c("DIA_BEGIN","DIA_END","PREVDIAmean","DIAmean","DIA","PREVDIA","DIA5alive","DIA5meanalive","PREVDIA5alive","PREVDIA5meanalive"))
# 
# n=100000
# plot(q$DIA_BEGIN[1:n], q$PREVDIA[1:n])
# plot(q$DIA_END[1:n], q$DIA[1:n])
# plot(q$DIA_BEGIN[1:n], q$PREVDIA5alive[1:n])
# plot(q$DIA_END[1:n], q$DIA5alive[1:n])
# plot(q$PREVDIAmean[1:n], q$PREVDIA5meanalive[1:n])
# plot(q$DIAmean[1:n], q$DIA5meanalive[1:n])
# 
# n=100000
# 
# hist(q$DIA_BEGIN[1:n]-q$PREVDIA[1:n], 1000)
# hist(q$DIA_END[1:n]-q$DIA[1:n], 1000)
# hist(q$DIA_BEGIN[1:n]-q$PREVDIA5alive[1:n], 1000)
# hist(q$DIA_END[1:n]-q$DIA5alive[1:n], 1000)
# hist(q$PREVDIAmean[1:n]-q$PREVDIA5meanalive[1:n], 1000)
# hist(q$DIAmean[1:n]-q$DIA5meanalive[1:n], 1000)
# 
# mean(q$DIA_BEGIN[1:n]-q$PREVDIA[1:n], na.rm=T)
# mean(q$DIA_END[1:n]-q$DIA[1:n], na.rm=T)
# mean(q$DIA_BEGIN[1:n]-q$PREVDIA5alive[1:n], na.rm=T)
# mean(q$DIA_END[1:n]-q$DIA5alive[1:n], na.rm=T)
# mean(q$PREVDIAmean[1:n]-q$PREVDIA5meanalive[1:n], na.rm=T)
# mean(q$DIAmean[1:n]-q$DIA5meanalive[1:n], na.rm=T)
# 
# 
# 
# # Note: DIA_END from GRM is the same as DIA from TREE. However, DIA_BEGIN is not necessarily the same as PREVDIA
