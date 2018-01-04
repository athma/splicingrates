library(stats)

# read in attributes
args = commandArgs(trailingOnly = T)
path = args[1]
name = args[2]
sub = args[3]

# Root solver function
sumsqequationsolve <- function(atts, txnrate){
  # atts[1:3] are Dprimes
  # atts[4:6] are Rprimes
  D_prime = atts[1:3]
  R_prime = atts[4:6]
  hold.row <- c(NA, NA)
  f <- function(h){ ((h*(1 - 2^(-D_prime[1]/(h*txnrate)))) - R_prime[1])^2 + ((h*(1 - 2^(-D_prime[2]/(h*txnrate)))) - R_prime[2])^2 + ((h*(1 - 2^(-D_prime[3]/(h*txnrate)))) - R_prime[3])^2 }
  starth = 0
  if(sum(is.na(R_prime))==3){ return(hold.row) }
  try(fit.hold <- optim(starth, f))
  try(hold.row <- c(fit.hold$par, fit.hold$value))
  return(hold.row)
}

##### JUNCTION METHOD #####

rl=50-10
txnrate=1500

min5 <- read.table(paste0(path,'/',name,'_5min_merged_sub',sub,'_startsites_junc.combo'),sep="\t",header=T)
fullintronlist = unique(as.character(min5$intron))

## get d_dist information
full.intron.threelength <- read.table('/nobackup1/athma/genomes/dm3/threeprime_distance.txt',sep="\t",header=T)

print("reading in files...")

## read in simulation files
#5m
minute = 5
min5$ratio <- (min5$ie_count/rl)/(min5$ee_count/rl)
min5$D_dist = full.intron.threelength$threelength[match(min5$intron, full.intron.threelength$intron)]
min5$D_prime = (min5$D_dist) + (minute + txnrate)
min5$R_prime = (min5$D_prime*log(2))/(txnrate*((1/min5$ratio) +1))

min5_ord = min5[match(fullintronlist, min5$intron),]

#10m
min10 <- read.table(paste0(path,'/',name,'_10min_merged_sub',sub,'_startsites_junc.combo'),sep="\t",header=T)
min10$ratio <- (min10$ie_count/rl)/(min10$ee_count/rl)
min10$D_dist = full.intron.threelength$threelength[match(min10$intron, full.intron.threelength$intron)]
min10$D_prime = (min10$D_dist) + (minute + txnrate)
min10$R_prime = (min10$D_prime*log(2))/(txnrate+((1/min10$ratio) + 1))

min10_ord = min10[match(fullintronlist, min10$intron),]

#20m
min20 <- read.table(paste0(path,'/',name,'_20min_merged_sub',sub,'_startsites_junc.combo'),sep="\t",header=T)
min20$ratio <- (min20$ie_count/rl)/(min20$ee_count/rl)
min20$D_dist = full.intron.threelength$threelength[match(min20$intron, full.intron.threelength$intron)]
min20$D_prime = (min20$D_dist) + (minute + txnrate)
min20$R_prime = (min20$D_prime*log(2))/(txnrate*((1/min20$ratio) + 1))

min20_ord = min20[match(fullintronlist, min20$intron),]

#dataframe with intron name, D_dist, 5/10/20 ee, 5/10/20 ie, 5/10/20 ratio

print("making combo dataframe...")
combo.juncratio.data <- data.frame(intron = fullintronlist,
                                   D_dist = full.intron.threelength$threelength[match(fullintronlist, full.intron.threelength$intron)],
                                   ie_count_5 = min5_ord$ie_count, ie_count_10 = min10_ord$ie_count, ie_count_20 = min20_ord$ie_count,
                                   ee_count_5 = min5_ord$ee_count, ee_count_10 = min10_ord$ee_count, ee_count_20 = min20_ord$ee_count,
                                   ratio_5 = min5_ord$ratio, ratio_10 = min10_ord$ratio, ratio_20 = min20_ord$ratio,
                                   D_prime_5 = min5_ord$D_prime, D_prime_10 = min10_ord$D_prime, D_prime_20 = min20_ord$D_prime,
                                   R_prime_5 = min5_ord$R_prime, R_prime_10 = min10_ord$R_prime, R_prime_20 = min20_ord$R_prime)

print("fitting data...")
sumsqfit.data <- t(apply(combo.juncratio.data[,c('D_prime_5', 'D_prime_10', 'D_prime_20',
                                                 'R_prime_5', 'R_prime_10', 'R_prime_20')],
                         1, sumsqequationsolve, txnrate))
combo.juncratio.data$root = sumsqfit.data[,1]
combo.juncratio.data$zerodev = sumsqfit.data[,2]

print("writing final file...")
juncfile = paste0(path,'/juncrates_sub',sub,'.txt')
write.table(combo.juncratio.data, file=juncfile, sep="\t", quote=F, row.names=F, col.names=T)
