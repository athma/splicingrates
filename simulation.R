library(stats)

# read in attributes
args = commandArgs(trailingOnly = T)
DIR = args[1]
D_dists = as.numeric(args[2])
expression_levels = as.numeric(args[3])
#labelings = as.numeric(args[4])

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

# Fragment and Read Function
get_reads <- function(lengths, eta_val = 200, insertsize = c(200, 300))
{
  # Select a fragment from each transcript, size select, and return the starting position of the
  # resulting reads relative to the length of the transcript
  # Inputs: 
  #   lengths - the lengths of the transcript
  #   eta_val - the eta value input to the Weibull distribution
  #   insertsize - a list of length two with the lower and upper bound for fragments
  # Outputs:
  #   fragments - a data frame with the start positions of the filtered fragments and the index
  #               of the transcript that it came from (columns: transcript, start)

  # sample lengths from a weibull distribution for each transcript and transform them to the length 
  # of the transcripts
  deltas = log10(lengths)
  ns_minus_1 = pmax(round(lengths/eta_val/gamma(1/deltas + 1)) - 1, 0)
  xis = lapply(ns_minus_1, function(n) {diff(sort(c(runif(n), 0, 1)))})
  xis_transformed = mapply(function(x, d) {x^(1/d)}, xis, deltas, SIMPLIFY = F)
  delta_is = mapply(function(len, x_t) {round(len*x_t/sum(x_t))}, lengths, xis_transformed, SIMPLIFY = F)

  # get all the start and end points of the fragments 
  starts = lapply(delta_is, function(d) {
    if (length(d) > 1) {
      c(0, cumsum(d[1:(length(d)-1)]))+1
    } else{
      1
    }
  })
  ends = lapply(delta_is, function(d) { cumsum(d)
    #if (length(d) > 1) {
    #  c(cumsum(d[1:(length(d)-1)]), sum(d)-sample(min(insertsize[1], sum(d) - d[length(d)]), 1))
    #} else{
    #  d
    #}
  })
  
  # convert to a data frame of fragments and associated transcript index
  #fragments = data.frame(transcript = rep(1:length(deltas), lengths(delta_is)),
  fragments = data.frame(transcript = rep(1:length(deltas), unlist(lapply(delta_is, length))),
                         start = unlist(starts),
                         end = unlist(ends))
  fragments$length = fragments$end - fragments$start

  # Filter fragments by length and return
  fragments = fragments[fragments$length >= insertsize[1] & fragments$length <= insertsize[2],]
  return(fragments[c('transcript', 'start')])
}

# Function to format reads
formatreads <- function(reads, intron, start_pos, U_dist, rl){

  reads_formatted <- data.frame(name = reads$name, flag = '0', ref = paste0('intron',as.character(intron)),
                                pos = reads$start+U_dist, mapq = '50', cigar = '50M',
                                rnext = '*', pnext = '0', tlen = '0', seq = paste(rep('A',50),collapse=''), qual = '*',
                                tag1 = 'XA:i:1', tag2 = 'MD:Z:50', tag3 = 'NM:i:0', tag4 = 'NH:i:1', tag5 = 'XS:A:+')
  junction = (reads$start > -rl + 10) & (reads$start <= -9) & start_pos$spliced
  cigar.vec = as.character(reads_formatted$cigar)
  is.junction = which(junction==T)
  cigar.vec[is.junction] = paste(as.character(-reads$start[is.junction]), 'M', as.character(intron), 'N',
                                 as.character(rl+reads$start[is.junction]), 'M', sep="")
  reads_formatted$cigar <- cigar.vec
  return(reads_formatted)
}

# Simulation Function
simulate <- function(intron, exon, U_dist, D_dist, labeling, half_life,
                     expression_level, n_millions, transcription_rate)
{
  # Inputs: 
  #   intron - length of the intron in kb
  #   exon - length of downstream exon in kb
  #   U_dist - distance from the TSS to the beginning of the intron
  #   D_dist - distance from the end of the exon to the PAS in kb
  #   labeling - the length of the labeling period in minutes
  #   half_life - the half-life time in minutes
  #   expression_level - the expression level in TPM
  #   n_millions - the total number of millions of transcripts to consider
  #   transcription_rate - rate of transcription in kb/min
  # Outputs: (in list form)
  #   [[1]] data frame of all reads with the start position, the matching pattern, and a 
  #         string with all the associated parameters (columns: start, match, name)
  #   [[2]] the number of spliced junction reads
  #   [[3]] the number of unspliced junction reads
  #   [[4]] whether or not each read comes from a spliced or unspliced transcript

  intron = 1000*intron
  exon   = 1000*exon
  U_dist = 1000*U_dist
  D_dist = 1000*D_dist
  transcription_rate = 1000*transcription_rate

  # Generate expression_level*n_millions transcripts uniformly from the labeled region 
  #end_sites = round(seq(from = 1, to = intron + exon + D_dist + labeling*transcription_rate,
  #                     length.out = expression_level*n_millions))
  end_sites = U_dist + sample.int(intron + exon + D_dist + labeling*transcription_rate,
                                  expression_level*n_millions, replace = TRUE)

  # Determine if the transcripts are spliced or not and their resulting lengths 
  spliced = ((runif(length(end_sites))>2^(-(end_sites-intron)/half_life/transcription_rate))&
               ((end_sites) > intron))
  lengths = U_dist + pmin(end_sites, intron + exon + D_dist) - intron*spliced

  # Determine the number of reads to generate from each transcript using a negative binomial model
  # read_nums = lapply(lengths, function(length) { rnbinom(n=1, size = length/15, mu = length/5) })
  # lengths_all = c()
  # for (i in 1:length(lengths)){
  #   lengths_all = c(lengths_all, matrix(lengths[i], 1, read_nums[[i]]))
  # }

  # Get the reads from the transcripts and map them to the gene
  rl = 50
  start_pos = get_reads(lengths)
  toss_distance = 300
  start_pos = start_pos[start_pos$start > toss_distance,]
  start_pos$spliced = spliced[start_pos$transcript]
  #reads = t(data.frame(row.names = c('start', 'match')))
  reads = data.frame(start = start_pos$start - U_dist +
                       intron*start_pos$spliced*(start_pos$start > U_dist))
  reads$name = paste('i', as.character(intron), 'e', as.character(exon),
                     'u', as.character(U_dist), 'd', as.character(D_dist),
                     'L', as.character(labeling), 'hl', as.character(half_life),
                     'X', as.character(expression_level),
                     'M', as.character(n_millions), sep=':')
  reads_formatted = formatreads(reads, intron, start_pos, U_dist, rl)
  # reads$match = '5M'
  # #junction = (reads$start > -rl + 1) & (reads$start <= 0) & start_pos$spliced
  reads$junction = (reads$start > -rl + 10) & (reads$start <= -9) & start_pos$spliced
  # reads$match[junction] = 
  #   paste(as.character(-reads$start[junction]), 'M', as.character(intron), 'N', 
  #         as.character(rl+reads$start[junction]), 'M', sep = '')
  # reads$name = paste('intron', as.character(intron), 'exon', as.character(exon), 
  #                    'U_dist', as.character(U_dist), 'D_dist', as.character(D_dist),
  #                    'labeling', as.character(labeling), 'half_life', as.character(half_life), 
  #                    'expression_level', as.character(expression_level), 
  #                    'n_millions', as.character(n_millions), sep=':')
  spliced_num = sum(reads$junction)
  unspliced_num = sum(!start_pos$spliced &
                        (start_pos$start > intron + U_dist-rl+10) &
                        (start_pos$start <= intron + U_dist-9))
  intron_num = sum(!start_pos$spliced &
                     (start_pos$start >= 1) & (start_pos$start < 50))

  return(list(reads_formatted, spliced_num, unspliced_num, intron_num, start_pos$spliced))
}

##### Grid of paramters to simulate over #####
## keep constant
exon = 0.3
U_dist = 0.5
n_millions = 100
transcription_rate = 1.5

## variables to simulate over
labelings = c(5,10,20,60)
#D_dists = seq(0.5,5,by=0.5)
#expression_levels = seq(1,100,by=4)
half_lives = c(seq(0.2,0.9,by=0.1),seq(1,10,by=0.75),seq(11,100,by=2))
introns = c(seq(0.04,0.09,by=0.02), seq(0.1,1,by=0.1), seq(1,50,by=2))

## create header for SAM file - based on intron lengths
max_nonintron_length = (U_dist + exon + D_dists[length(D_dists)])*1000
sam_header = t(data.frame(row.names = c('at', 'sn', 'ln')))
for(intron in introns){
  sam_header = rbind(sam_header,
                     data.frame(at = '@SQ',
                                sn = paste0('SN:intron',as.character(intron*1000)),
                                ln = paste0('LN:',as.character(intron*1000+max_nonintron_length+1))))
}

##### Run the grid simulation #####
simulation_results = t(data.frame(row.names = c('spliced_num', 'unspliced_num', 'params')))
for(labeling in labelings){
  for (D_dist in D_dists){
    for (expression_level in expression_levels){
      for (half_life in half_lives){
        # create file to add reads across all introns and add header
        filename = paste0(DIR,'m',labeling,'/D',D_dist,'_X',expression_level,'_H',half_life,'_reads.sam')
        if(labeling<= 20){ write.table(sam_header, file=filename, sep="\t",quote=F,row.names=F,col.names=F) }
        print(paste0("minute ",labeling," - distance ",D_dist," - expression ",expression_level," - half-life ",half_life,"..."))
        for (intron in introns){
          cat(paste0(intron," - "))
          if(intron %% 1 == 0){ flush.console() }
          results = simulate(intron, exon, U_dist, D_dist, labeling, half_life,
                                     expression_level, n_millions, transcription_rate)
          reads = results[[1]]
          # reads: dataframe in read format - can use for MISO (for PSI decrease method) - chr is intron length
          if(labeling <= 20){ write.table(reads, file=filename, sep="\t",quote=F, row.names=F,col.names=F, append=T) }
          spliced_num = results[[2]]
          unspliced_num = results[[3]]
          intron_num = results[[4]]
          simulation_results = rbind(simulation_results,
                                     data.frame(spliced_num = spliced_num,
                                                unspliced_num = unspliced_num,
                                                intron_num = intron_num,
                                                intron = intron,
                                                D_dist = D_dist,
                                                labeling = labeling,
                                                half_life = half_life,
                                                exon = exon,
                                                U_dist = U_dist,
                                                expression_level = expression_level,
                                                n_millions = n_millions,
                                                transcription_rate = transcription_rate))
        }
        print("")
      }
    }
  }
}

##### JUNCTION DYNAMICS METHOD #####
# Use the results to get an estimate for the half life
simulation_results$intron <- simulation_results$intron*1000
simulation_results$D_dist <- simulation_results$D_dist*1000
simulation_results$exon <- simulation_results$exon*1000
simulation_results$U_dist <- simulation_results$U_dist*1000
simulation_results$transcription_rate <- simulation_results$transcription_rate*1000
simulation_results$ratio <- simulation_results$unspliced_num/simulation_results$spliced_num
simulation_results$D_prime <- (simulation_results$exon + simulation_results$D_dist) + (simulation_results$labeling*simulation_results$transcription_rate)
simulation_results$R_prime <- (simulation_results$D_prime*log(2)) / (simulation_results$transcription_rate*((1/simulation_results$ratio) + 1))

#simfile = paste0(DIR,"SIM_D",D_dists[1],"_X",expression_levels[1],"_L",labelings[1],".txt")
#write.table(simulation_results, file=simfile, sep="\t",quote=F, row.names=F, col.names=T)

simulation_results_5 <- subset(simulation_results, labeling==5)
simulation_results_10 <- subset(simulation_results, labeling==10)
simulation_results_20 <- subset(simulation_results, labeling==20)

# make one file across all timepoints
results_wide <- simulation_results_5[,c(4,5,7:12)]
results_wide <- data.frame(results_wide,
                           spliced_5 = simulation_results_5$spliced_num, spliced_10 = simulation_results_10$spliced_num, spliced_20 = simulation_results_20$spliced_num,
                           unspliced_5 = simulation_results_5$unspliced_num, unspliced_10 = simulation_results_10$unspliced_num, unspliced_20 = simulation_results_20$unspliced_num,
                           ratio_5 = simulation_results_5$ratio, ratio_10 = simulation_results_10$ratio, ratio_20 = simulation_results_20$ratio,
                           D_prime_5 = simulation_results_5$D_prime, D_prime_10 = simulation_results_10$D_prime, D_prime_20 = simulation_results_20$D_prime,
                           R_prime_5 = simulation_results_5$R_prime, R_prime_10 = simulation_results_10$R_prime, R_prime_20 = simulation_results_20$R_prime)


# calculate and add intron/exon density ratio to the data frame
sumsqfit.data <- t(apply(results_wide[,c('D_prime_5', 'D_prime_10', 'D_prime_20',
                                         'R_prime_5', 'R_prime_10', 'R_prime_20')],
                         1, sumsqequationsolve, results_wide$transcription_rate[1]))
results_wide$root = sumsqfit.data[,1]
results_wide$zerodev = sumsqfit.data[,2]

juncfile = paste0(DIR,"juncSIM_D",D_dists[1],"_X",expression_levels[1],".txt")
write.table(results_wide, file=juncfile, sep="\t", quote=F, row.names=F, col.names=T)

##### INTRON RATIO METHOD #####

simulation_results_60 <- subset(simulation_results, labeling==60)

results_ratio <- simulation_results_5[,c(4,5,7:12)]
results_ratio <- data.frame(results_ratio,
                            intron_5 = simulation_results_5$intron_num, intron_60 = simulation_results_60$intron_num)

results_ratio$splicingratio = results_ratio$intron_60 / results_ratio$intron_5

ratiofile = paste0(DIR,"ratioSIM_D",D_dists[1],"_X",expression_levels[1],".txt")
write.table(results_ratio, file=ratiofile, sep="\t", quote=F,row.names=F,col.names=T)


