library(stats)

# read in attributes
args = commandArgs(trailingOnly = T)
path = args[1]
D_dist = args[2]

# list of half_lives
half_lives <- c(seq(0.2,0.9,by=0.1),seq(1,10,by=0.75),seq(11,100,by=2))
expression_levels <- seq(1,50,5)

# list of introns
introns = c(seq(0.04,0.09,by=0.02), seq(0.1,1,by=0.1), seq(1,50,by=2))
uplength = 500
downlength = 300

intronnames <- c()
for(intron in introns){
  #print(intron)
  chr = paste0("intron",(intron*1000))
  genelength = uplength + downlength + (intron*1000)
  downexon_start = uplength + (intron*1000)
  # gene row
  intronnames <- c(intronnames,
                   paste0("chr",chr,":1:",uplength,":+@chr",chr,":",downexon_start,":",genelength,":+"))
}

getLM_throughzero_tau_lowadj <- function(input){
  vals <- as.numeric(input[1:3])
  distance <- as.numeric(input[4])*1000
  # assume txn rate is 1.5kb/minute
  tau = distance/1500
  # get corrected time
  times <- c(5,10,20)
  times_with <- (times+tau)/2
  lm.row = c(NA,NA,NA,NA)
  try(lm.hold <- lm(log(vals)~0+times_with))
  try(lm.row <- c(as.numeric(lm.hold$coefficients[1]),summary(lm.hold)$coefficients[1,4],
                  summary(lm.hold)$adj.r.squared, "fit"))
  # correct for low half-life; if second effective timepoint should have PSI < 0.1
  k <- -log(vals[1])/times_with[1]
  meansecond <- exp(-k*times_with[2])
  try(if(meansecond <= 0.1){ lm.row <- c(-k, NA, NA, "corrected") })
  return(lm.row)
}

get_halflife <- function(x){
  #x = [constant]
  #half.val <- 1-x[3]
  #return((x[2] + log(half.val))/x[1])
  return(log(2)/(-x))
}

## go through all half-lives for given D_dist
psi.data.all <- c()

for(exp in expression_levels){
        for(halflife in half_lives){
                print(paste0(exp," - ",halflife))
                # name files
                min5.name <- paste0(path,'summary/m5/m5_D',D_dist,'_X',exp,'_H',halflife,'_reads_sorted_introns.psi.miso_summary')
                min10.name <- paste0(path,'summary/m10/m10_D',D_dist,'_X',exp,'_H',halflife,'_reads_sorted_introns.psi.miso_summary')
                min20.name <- paste0(path,'summary/m20/m20_D',D_dist,'_X',exp,'_H',halflife,'_reads_sorted_introns.psi.miso_summary')
                # read in files
                if(file.exists(min5.name) & file.exists(min10.name) & file.exists(min20.name)){
                        min5 <- read.table(min5.name,header=T,sep="\t")
                        min10 <- read.table(min10.name,header=T,sep="\t")
                        min20 <- read.table(min20.name,header=T,sep="\t")
                        # make dataframe of data
                        psi.data <- data.frame(intron = intronnames,
                                               intron_len = introns,
                                               D_dist = D_dist,
                                               expression = exp,
                                               sim_half_life = halflife,
                                               psi5 = min5$miso_posterior_mean[match(intronnames, min5$event_name)],
                                               psi10 = min10$miso_posterior_mean[match(intronnames, min10$event_name)],
                                               psi20 = min20$miso_posterior_mean[match(intronnames, min20$event_name)])
                        # get rates
                        psi.data$coef <- t(apply(psi.data[,c(6:8,3)], 1, getLM_throughzero_tau_lowadj))[,1]
                        psi.data$halflife <- sapply(as.numeric(psi.data$coef), get_halflife)
                        # add to overall dataframe
                        psi.data.all <- rbind(psi.data.all, psi.data)
                }
                else { print("... no files!") }
        }
}

psi.file <- paste0(path,'psiSIM_D',D_dist,'.txt')
write.table(psi.data.all, file=psi.file, quote=F, sep="\t",row.names=F,col.names=T)
