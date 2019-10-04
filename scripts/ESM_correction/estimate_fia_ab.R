library(PEcAn.allometry)
library(data.table)
library(plyr)
library(methods)

file.in.normal = "../../data/ESM.rds"
file.in.fading = "../../data/ESM_F.rds"

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

############################ build ab_lut (loop over FIA plots) ############################ 

# Should I create 1 lookup table or 4 (normal ESM: present-past, fading record ESM: present-past)
# For now trying 1

fia_plots <- unique(q$plt_cn) # 64520 plots

ab_lut  <- as.data.frame(matrix(NA, nrow = length(fia_plots), ncol= 4))
colnames(ab_lut) <- c("plt_cn","dbh", "dens", "ab")

for(p in seq_along(fia_plots)){
  # develop/debug
  # p <- 5
  sub_plot <- q[q$plt_cn == fia_plots[p], ] 
  if(all(is.na(sub_plot$dia_end))) next
  
  beta0s <- chojnacky_table$beta0[match(sub_plot$paleon, as.character(chojnacky_table$level3a))]
  beta1s <- chojnacky_table$beta1[match(sub_plot$paleon, as.character(chojnacky_table$level3a))]
  
  ### exp(b0 + b1*log(dbh))
  agb           <- sum(exp(beta0s + beta1s*log(sub_plot$dia_end * 2.54)), na.rm = TRUE)
  agb_Mg        <- agb * 1e-03 # kg to Mg
  agb_Mg_plt    <- agb_Mg / (4*7.3152^2*pi) # Mg/m2
  agb_Mg_ha_plt <- agb_Mg_plt / 1e-04  # Mg/m2 to Mg/ha
  
  ab_lut[p, 1] <- sub_plot$plt_cn[1]    # plot no
  ab_lut[p, 2] <- sub_plot$diamean.m[1] # mean dbh (cm)
  ab_lut[p, 3] <- sub_plot$tpasum.m[1]
  ab_lut[p, 4] <- agb_Mg_ha_plt
  
  if(p %% 100 == 0) print(p)
}

#some plots had no dbh data
ab_lut <- ab_lut[complete.cases(ab_lut),]
save(ab_lut, file = "../../data/ab_lut.Rdata")


