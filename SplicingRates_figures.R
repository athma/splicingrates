setwd("~/Desktop/Dropbox (MIT)/Projects/Adelman/timecourse")
#load("finalcode_modelPSI_txpttxn.RData")

library(ggplot2)
library(wesanderson)
library(GGally)
library(corrplot)
library(cowplot)
library(scales)
library(ggrepel)
library(MatchIt)
library(seqLogo)
library(png)
library(grid)
library(RColorBrewer)

g.legend <- function(a.plot){
  tmp <- ggplot_gtable(ggplot_build(a.plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

##### FIGURE 1 - methods #####
#### A - progressive labeling & B - approaches ####
splicingschematic_img <- readPNG("Figures/revisedfigures/SplicingRate_schematics.png")
splicingschematic_grob <- rasterGrob(splicingschematic_img, interpolate=T)
f1_ab <- qplot(1:10, 1:10, geom="blank") + 
  annotation_custom(splicingschematic_grob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  labs(x="",y="") + theme_classic() + theme(axis.text=element_blank(), axis.ticks=element_blank())
#### C - simulation correlations ####
f1_c <- ggplot(subset(full.sim.cors.data, cor_type=="spearman"), aes(x=factor(type),y=meancor_hl20, fill=factor(type),color=factor(type))) + 
  geom_bar(stat="identity",position="dodge",color=NA) + geom_errorbar(aes(ymin=meancor_hl20-secor_hl20, ymax=meancor_hl20+secor_hl20),width=0.5) + 
  scale_fill_manual(values=c(brewer.pal(9,"RdPu")[3], brewer.pal(9,"Oranges")[3], brewer.pal(9,"BuPu")[6]),guide=F) +
  scale_color_manual(values=c(brewer.pal(9,"RdPu")[5], brewer.pal(9,"Oranges")[5], brewer.pal(9,"BuPu")[8]),guide=F) +
  scale_x_discrete(labels=c("intron ratio",expression(paste(Psi," decrease",sep="")),"junction dynamics")) + scale_y_continuous(limits=c(0,1.15),breaks=c(0,0.25,0.5,0.75,1)) +
  labs(x="approach",y=expression(paste("correlation btwn est and simulated t"["1/2"]))) + background_grid(major="y",minor="none") + 
  theme(axis.text.x=element_text(size=8,angle=45,hjust=1),axis.text.y=element_text(size=8),axis.title.x=element_blank(),axis.title.y=element_text(size=10)) 
#### D - simulation % errors ####
f1_d <- ggplot(fullsims.error, aes(x=factor(type),y=abslogpererror_hl20_mean, fill=factor(type),color=factor(type))) + 
  geom_bar(stat="identity",position="dodge",color=NA) + geom_text(x=1,y=0.05,label="NA",color=brewer.pal(9,"RdPu")[3],size=3) +
  geom_errorbar(aes(ymin=abslogpererror_hl20_mean-abslogpererror_hl20_se, ymax=abslogpererror_hl20_mean+abslogpererror_hl20_se),width=0.5) + 
  scale_fill_manual(values=c(brewer.pal(9,"RdPu")[3], brewer.pal(9,"Oranges")[3], brewer.pal(9,"BuPu")[6]),guide=F) +
  scale_color_manual(values=c(brewer.pal(9,"RdPu")[5], brewer.pal(9,"Oranges")[5], brewer.pal(9,"BuPu")[8]),guide=F) +
  scale_x_discrete(labels=c("intron ratio",expression(paste(Psi," decrease",sep="")),"junction dynamics")) + #scale_y_continuous(limits=c(0,10),breaks=c(-50,0,50,100)) +
  labs(x="approach",y=expression(paste("| log2(est t" ["1/2"], " / sim t" ["1/2"], " ) |"))) +  background_grid(major="y",minor="none") +
  theme(axis.text.x=element_text(size=8,angle=45,hjust=1),axis.text.y=element_text(size=8),axis.title.x=element_blank(),axis.title.y=element_text(size=10))
#### E - % relative truth ####
f1_e <- ggplot(fullsims.relative.error.hl20, aes(x=factor(type),y=abslogpererror_mean, fill=factor(type),color=factor(type))) + 
  geom_bar(stat="identity",position="dodge",color=NA) + 
  geom_errorbar(aes(ymin=abslogpererror_mean-abslogpererror_se, ymax=abslogpererror_mean+abslogpererror_se),width=0.5) + 
  scale_fill_manual(values=c(brewer.pal(9,"RdPu")[3], brewer.pal(9,"Oranges")[3], brewer.pal(9,"BuPu")[6]),guide=F) +
  scale_color_manual(values=c(brewer.pal(9,"RdPu")[5], brewer.pal(9,"Oranges")[5], brewer.pal(9,"BuPu")[8]),guide=F) +
  scale_x_discrete(labels=c("intron ratio",expression(paste(Psi," decrease",sep="")),"junction dynamics")) + #scale_y_continuous(limits=c(0,10),breaks=c(-50,0,50,100)) +
  labs(x="approach",y=expression(paste("| log2( est t" ["1/2"]^"j", "/t" ["1/2"]^"k", " / sim t" ["1/2"]^"j", "/t" ["1/2"]^"k", ") |"))) +  background_grid(major="y",minor="none") +
  theme(axis.text.x=element_text(size=8,angle=45,hjust=1),axis.text.y=element_text(size=8),axis.title.x=element_blank(),axis.title.y=element_text(size=10))
#### F - simulation - across half-lives ####
f1_f <- ggplot(fullsims, aes(x=half_life, y=sim_hl,fill=factor(type),color=factor(type))) + 
  geom_abline(color="lightgoldenrod1",size=3) + 
  geom_point(data=subset(fullsims, type=="psi"),alpha=0.05,shape=21,color=NA) + 
  geom_point(data=subset(fullsims, type=="junc"),alpha=0.05,shape=21,color=NA) + 
  stat_smooth(data=subset(fullsims, type=="psi"),method="lm",size=0.5) + stat_smooth(data=subset(fullsims, type=="junc"),method="lm",size=0.5) +
  scale_fill_manual(values=c(brewer.pal(9,"BuPu")[6], brewer.pal(9,"Oranges")[3]),labels=c("junction dynamics",expression(paste(Psi," decrease",sep="")))) +
  scale_color_manual(values=c(brewer.pal(9,"BuPu")[8], brewer.pal(9,"Oranges")[5]),labels=c("junction dynamics",expression(paste(Psi," decrease",sep="")))) +
  annotate("segment",x=12.9,xend=13.65,y=1.45,yend=1.45,color="lightgoldenrod1",size=3) +
  annotate("text",x=16.5,y=1.45,label="y = x line",size=3) +
  xlim(0,20) + ylim(0,20) + guides(fill=guide_legend(reverse=T),color=guide_legend(reverse=T)) +
  labs(x="half-life used in simulations (min)",y="estimated half-life (min)",fill="approach",color="approach") + 
  theme(legend.position=c(0.8,0.2),legend.text.align=0,legend.key.size=unit(3,"mm"),legend.title=element_text(size=10),legend.text=element_text(size=9),axis.text=element_text(size=9),axis.title=element_text(size=11)) + 
  background_grid(major="xy",minor="none")
#### G - example (short intron) ####
short.intron <- "chr3L:8968425:8968774:+@chr3L:8968837:8969499:+"
short.intron.filename <- "chr3L:8968425:8968774:+@chr3L:8968837:8969499:+.cov"
short.cov.file <- read.table(paste0("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/coverage_examples/short/coveragefiles/",short.intron.filename),header=T)
short.cov.file$time <- factor(short.cov.file$time, levels=c("5min","10min","20min","total"))
short.cov.file$type <- factor(short.cov.file$type, levels=c("upexon","intron","downexon"))
# get starting position of upstream exon
short.start.pos <- as.numeric(strsplit(as.character(short.intron),split=":")[[1]][2])
# get region for zoommed in plot
short.intron.min <- min(subset(short.cov.file, type=="intron")$position) - 100
short.intron.max <- max(subset(short.cov.file, type=="intron")$position) + 100

f1_g1 <- ggplot(short.cov.file, aes(x=position+short.start.pos,y=coverage_norm,fill=factor(time))) + geom_bar(data=subset(short.cov.file, type=="intron"),stat="identity") +  
  geom_bar(data=subset(short.cov.file, type!="intron"),stat="identity",size=0.75,aes(color=factor(time))) + scale_x_continuous(limits=c(short.intron.min+short.start.pos, short.intron.max+short.start.pos),labels=comma) + 
  scale_fill_manual(values=c(rev(brewer.pal(9,"Purples"))[2:4],"goldenrod"),guide=F) + scale_color_manual(values=rev(brewer.pal(9,"RdGy"))[1:4],guide=F) + scale_y_continuous(label=comma) + 
  labs(x="chromosome 3L, intron 2 of Srp68",y="coverage",fill="timepoint",color="timepoint") + 
  theme(legend.position=c(0.85,0.9),axis.text=element_text(size=7),axis.title=element_text(size=9),legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size=unit(3,"mm"))
legend.covplots1 <- ggplot(short.cov.file, aes(x=position+short.start.pos,y=coverage_norm,fill=factor(time))) + geom_bar(data=subset(short.cov.file, type=="intron"),stat="identity") +  
  geom_bar(data=subset(short.cov.file, type!="intron"),stat="identity",size=0.75,aes(color=factor(time))) + scale_x_continuous(limits=c(short.intron.min+short.start.pos, short.intron.max+short.start.pos),labels=comma) + 
  scale_fill_manual(values=c(rev(brewer.pal(9,"Purples"))[2:4],"goldenrod"),guide=guide_legend(nrow=1)) + scale_color_manual(values=rev(brewer.pal(9,"RdGy"))[1:4]) + scale_y_continuous(label=comma) + 
  labs(x="chromosome 3L, intron 2 of Srp68",y="coverage",fill="timepoint",color="timepoint") + 
  theme(legend.position="bottom",axis.text=element_text(size=7),axis.title=element_text(size=9),legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size=unit(3,"mm"))
f1_g1_legend <- g.legend(legend.covplots1)

short.example <- subset(combo.juncratio.data.parsed, intron==short.intron)
short.example.data <- data.frame(intron = rep(short.example$intron, 4),
                                 ee_count = c(short.example$ee_count_5, short.example$ee_count_10, short.example$ee_count_20, short.example$ee_count_total),
                                 ie_count = c(short.example$ie_count_5, short.example$ie_count_10, short.example$ie_count_20, short.example$ie_count_total),
                                 time = rep(c(5, 10, 20, 25), each=nrow(short.example)))
short.example.data$sum <- short.example.data$ee_count + short.example.data$ie_count

f1_g2 <- ggplot(short.example.data, aes(x=time)) + geom_point(aes(y=(ee_count/sum)*100,color="blue"),size=2) + geom_point(aes(y=(ie_count/sum)*100,color="maroon1"), size=2) + 
  scale_x_continuous(breaks=c(0,5,10,20,25),labels=c("","5","10","20","total")) + ylim(0,100) +
  scale_color_manual(values=c("blue"="blue","maroon1"="maroon1"),labels=c("exon-exon junctions","intron-exon junctions"),guide=F) +
  geom_vline(xintercept=22.5, color="grey17",linetype="dotted") + labs(x="labeling period (L)",y="junction reads (%)",color="") +
  theme(legend.position="bottom",axis.text=element_text(size=7),axis.title=element_text(size=9),legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size=unit(3,"mm"))
legend.covplots2 <- ggplot(short.example.data, aes(x=time)) + geom_point(aes(y=(ee_count/sum)*100,color="blue"),size=2) + geom_point(aes(y=(ie_count/sum)*100,color="maroon1"), size=2) + 
  scale_x_continuous(breaks=c(0,5,10,20,25),labels=c("","5","10","20","total")) + ylim(0,100) +
  scale_color_manual(values=c("blue"="blue","maroon1"="maroon1"),labels=c("exon-exon junctions","intron-exon junctions"),guide=guide_legend(nrow=1)) +
  geom_vline(xintercept=22.5, color="grey17",linetype="dotted") + labs(x="labeling period (L)",y="junction reads (%)",color="") +
  theme(legend.position="bottom",axis.text=element_text(size=7),axis.title=element_text(size=9),legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size=unit(3,"mm"))
f1_g2_legend <- g.legend(legend.covplots2)
#### H - example (long intron) ####

#long1 <- "chr3R:17093342:17093462:-@chr3R:17094493:17094782:-" #Rab1 - 2.44m
#long2 <- "chr3R:17646369:17646745:-@chr3R:17648205:17648396:-" #CASK - 3.15m

long.intron <- "chr3R:17093342:17093462:-@chr3R:17094493:17094782:-"
long.intron.filename <- paste(long.intron,"cov",sep=".")
long.cov.file <- read.table(paste0("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/coverage_examples/long_minus/coveragefiles/",long.intron.filename),header=T)
long.cov.file$time <- factor(long.cov.file$time, levels=c("5min","10min","20min","total"))
long.cov.file$type <- factor(long.cov.file$type, levels=c("upexon","intron","downexon"))
# get starting position of upstream exon
long.start.pos <- as.numeric(strsplit(as.character(long.intron),split=":")[[1]][2])
# get region for zoommed in plot
long.intron.min <- min(subset(long.cov.file, type=="intron")$position) - 100
long.intron.max <- max(subset(long.cov.file, type=="intron")$position) + 100

f1_h1 <- ggplot(long.cov.file, aes(x=position+long.start.pos,y=coverage_norm,fill=factor(time))) + geom_bar(data=subset(long.cov.file, type=="intron"),stat="identity") +  
  geom_bar(data=subset(long.cov.file, type!="intron"),stat="identity",size=0.75,aes(color=factor(time))) + scale_x_continuous(limits=c(long.intron.min+long.start.pos, long.intron.max+long.start.pos),labels=comma) + 
  scale_fill_manual(values=c(rev(brewer.pal(9,"Purples"))[2:4],"goldenrod"),guide=F) + scale_color_manual(values=rev(brewer.pal(9,"RdGy"))[1:4],guide=F) + scale_y_continuous(label=comma) + 
  labs(x="chromosome 3R, intron 1 of Rab1 (transcribed from Crick strand)",y="coverage",fill="timepoint",color="timepoint") + 
  theme(legend.position=c(0.85,0.9),axis.text=element_text(size=7),axis.title=element_text(size=9),legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size=unit(3,"mm"))

long.example <- subset(combo.juncratio.data.parsed, intron==long.intron)
long.example.data <- data.frame(intron = rep(long.example$intron, 4),
                                 ee_count = c(long.example$ee_count_5, long.example$ee_count_10, long.example$ee_count_20, long.example$ee_count_total),
                                 ie_count = c(long.example$ie_count_5, long.example$ie_count_10, long.example$ie_count_20, long.example$ie_count_total),
                                 time = rep(c(5, 10, 20, 25), each=nrow(long.example)))
long.example.data$sum <- long.example.data$ee_count + long.example.data$ie_count

f1_h2 <- ggplot(long.example.data, aes(x=time)) + geom_point(aes(y=(ee_count/sum)*100,color="blue"),size=2) + geom_point(aes(y=(ie_count/sum)*100,color="maroon1"), size=2) + 
  scale_x_continuous(breaks=c(0,5,10,20,25),labels=c("","5","10","20","total")) + ylim(0,100) +
  scale_color_manual(values=c("blue"="blue","maroon1"="maroon1"),labels=c("exon-exon junctions","intron-exon junctions"),guide=F) +
  geom_vline(xintercept=22.5, color="grey17",linetype="dotted") + labs(x="labeling period (L)",y="junction reads (%)",color="") +
  theme(legend.position="bottom",axis.text=element_text(size=7),axis.title=element_text(size=9),legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size=unit(3,"mm"))
#### COMBINE ####
pdf("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/fig1.pdf",width=8,height=10.5, useDingbats = F)
ggdraw() + # 
  draw_plot(f1_ab, 0, 0.6, 1, 0.5) +
  draw_plot(f1_c, 0, 0.39, 0.1666, 0.33) + draw_plot(f1_d, 0.1666, 0.39, 0.1666, 0.33) + draw_plot(f1_e, 0.3333, 0.39, 0.1666, 0.33) +
  draw_plot(f1_f, 0.5, 0.39, 0.5, 0.33) +
  draw_plot(f1_g1, 0, 0.205, 0.75, 0.185) + draw_plot(f1_g2, 0.73, 0.205, 0.25, 0.185) +
  draw_plot(f1_g1_legend, 0, 0.185, 0.75, 0.02) + draw_plot(f1_g2_legend, 0.735, 0.185, 0.25, 0.02) +
  draw_plot(f1_h1, 0, 0, 0.75, 0.185) + draw_plot(f1_h2, 0.73, 0, 0.25, 0.185) +
  annotate("text",x=c(0.4,0.4),y=c(0.32,0.12),label=c("63 nt\nt1/2 = 1.31 min +/- 2.6 sec","1,711 nt\nt1/2 = 2.44 min +/- 0.75 sec"),size=2.5) +
  draw_plot_label(c("A","B","C","D","E","F","G","H"), c(0, 0.6, 0, 0.1666, 0.3333, 0.5, 0, 0), c(1, 1, 0.72, 0.72, 0.72, 0.72, 0.39, 0.185))
dev.off()

##### FIGURE 2 - optimality #####
#### A - overall length dist ####
length.min <- combo.median.data.200txpts$length[10]
length.max <- combo.median.data.200txpts$length[nrow(combo.median.data.200txpts)-10]
f2_a <- ggplot(combo.median.data.200txpts, aes(x=length,y=half)) + geom_point(alpha=0.25, color="grey25") + 
  geom_point(data=subset(combo.median.data.200txpts, length > length.min & length < length.max),aes(x=length, y=fit_half),color=wes_palette("FantasticFox")[3]) + 
  scale_y_log10(lim=c(0.1,100),breaks=c(0.25,1,4,16,64,256,1024)) + scale_x_log10(limits=c(40,100000),breaks=c(40,60,100,1000,10000,100000),labels=comma) + labs(x="intron length (nt)",y="half-life (min)") + 
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),axis.text.y=element_text(size=8))
#### B - focus on 60-70 nt ####
f2_b_box <- ggplot(subset(combo.juncratio.data.parsed, intronlen <=100), aes(x=factor(lenbin5),y=fitvalue)) + geom_boxplot(notch=T,fill=wes_palette("FantasticFox")[3],outlier.color="lightgrey") + 
  scale_y_log10(limits=c(0.25,20),breaks=c(0.25,1,4,16,64,256,1024)) + labs(x="",y="half-life (min)") + theme(axis.text.x=element_blank(),axis.title.y=element_text(color=wes_palette("FantasticFox")[3]))
f2_b_bar <- ggplot(combo.lenbin.data[c(1:6),], aes(x=lenbin, y=counts)) + geom_bar(stat="identity",fill="goldenrod1",alpha=0.35) + labs(x="intron length (nt)", y="intron count") +  
  scale_y_log10(limits=c(1,10500),breaks=c(10,100,1000,10000),labels=comma) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),axis.text.y=element_text(color="darkgrey"),axis.title.y=element_text(color="goldenrod1",angle=270),axis.ticks.y=element_line(color="darkgrey"))
#### C - branchpoint distance ####

f2_c <- ggplot(bpscan_bar, aes(x=factor(bin),y=-mean)) + geom_bar(stat="identity", fill=wes_palette("FantasticFox")[3],width=0.75) + geom_errorbar(aes(ymin=-mean-se, ymax=-mean+se),width=0.5) + ylim(-40,0) + coord_flip() + 
  annotate("text",c(1,2,3,4,5,6),y=-5,label=rev(c("40-50nt","50-60nt","60-70nt","70-80nt","80-90nt","90-100nt")),color="white") +
  labs(y="distance from 3' splice site",x="intron length") + theme(axis.title.y=element_text(angle=270),axis.text.x=element_text(size=8),
                                                                   axis.text.y=element_text(size=0),axis.line.y=element_line(size=0),axis.ticks.y=element_line(size=0)) +
  background_grid(major="x")
#### D - regulatory mode ####
combo.juncratio.data.parsed$len_bin <- factor(combo.juncratio.data.parsed$len_bin, levels=c("20%","40%","60%","80%","100%"))
f2_d <- ggplot(subset(combo.juncratio.data.parsed, type!="SEcontaining"), aes(x=factor(len_bin),y=fitvalue, fill=factor(type))) + geom_boxplot(notch=T,outlier.color="lightgrey")+ 
  scale_y_log10(limits=c(0.25,10),breaks=c(0.5,1,2,4,8)) + scale_x_discrete(labels=c("40-60nt","61-65nt","66-77nt","78-284nt",">284nt")) + 
  scale_fill_manual(values=wes_palette("Darjeeling")[c(5,4,3)]) + labs(x="intron length (quintiles)",y="half-life (min)",fill=NULL) + 
  theme_cowplot() + background_grid(major="y",minor="y") + 
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),legend.direction="vertical",legend.position="right",legend.key.size=unit(5,"mm"),legend.text=element_text(size=10))
#### COMBINE ####
pdf("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/fig2.pdf",width=8,height=10.5, useDingbats = F)
ggdraw() + 
  draw_plot(f2_a, 0, 0.66, 0.5, 0.33) + 
  draw_plot(switch_axis_position(f2_b_bar, axis='y'), 0.0955, 0.33, 0.42, 0.33) + draw_plot(f2_b_box, 0.0, 0.364, 0.407, 0.296) + 
  draw_plot(switch_axis_position(f2_c, axis='y'), 0.5, 0.66, 0.5, 0.33) + 
  draw_plot(f2_d, 0.5, 0.33, 0.5, 0.33) +
  draw_plot_label(c("A","C","B","D"), c(0, 0.5, 0, 0.5), c(1, 1, 0.66, 0.66))
dev.off()  

##### FIGURE 3 - definition #####
#### A - schematic + circ plot ####
pdf("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/fig2_circplot.pdf",width=4,height=3.5)
pdf("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/fig2_circplot_40radial_10L_yellowpurple4.pdf",width=4,height=3.5)
ggplot(circ.data.format, aes(x=RIMEpolar, y=factor(circ))) + geom_tile(aes(fill=hl_bin, width=RIMEwidth)) + coord_polar(theta="x",start=0) +
  scale_fill_gradient(low="yellow",high="magenta4",na.value="black") + geom_segment(aes(x=0.122,xend=0.122,y=0,yend=10.5),color="gold",size=0.75) +
  labs(x="percentiles of RImE",y="percentiles of aggregated intron & exon length",fill="percentiles of intron half-life") +
  theme(legend.position=c(0.25,0.25),axis.text.x=element_blank())
dev.off()

# 40 radial bins, 10 linear bins
# added manually from illustrator
rimeplot_img <- readPNG("Figures/revisedfigures/circplot-02.png")
rimeplot_grob <- rasterGrob(rimeplot_img, interpolate=T)
f3_a <- qplot(1:10, 1:10, geom="blank") + 
  annotation_custom(rimeplot_grob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  labs(x="",y="") + theme_classic() + theme(axis.text=element_blank(), axis.ticks=element_blank())

binplot_img <- readPNG("Figures/revisedfigures/circplot-04.png")
binplot_grob <- rasterGrob(binplot_img, interpolate=T)
f3_b <- qplot(1:10, 1:10, geom="blank") + 
  annotation_custom(binplot_grob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  labs(x="",y="") + theme_classic() + theme(axis.text=element_blank(), axis.ticks=element_blank())

circplot_img <- readPNG("Figures/revisedfigures/circplot-03.png")
circplot_grob <- rasterGrob(circplot_img, interpolate=T)
f3_c <- qplot(1:10, 1:10, geom="blank") + 
  annotation_custom(circplot_grob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  labs(x="",y="") + theme_classic() + theme(axis.text=element_blank(), axis.ticks=element_blank())
#### D/E - RIME scatter ####

rine.min <- combo.median.data.200txpts.RINE$IEratio[10]
rine.max <- combo.median.data.200txpts.RINE$IEratio[nrow(combo.median.data.200txpts.RINE)-10]
f3_d <- ggplot(combo.median.data.200txpts.RINE, aes(x=IEratio,y=half)) + geom_point(alpha=0.15,size=0.5,shape=19,color="grey27") + 
  scale_color_gradient2(low="deeppink4",mid="snow2",high="dodgerblue4",guide=F) +
  geom_vline(xintercept=1,color="gold",linetype="dashed",alpha=0.95) + 
  geom_point(data=subset(combo.median.data.200txpts.RINE, IEratio > rine.min & IEratio < rine.max), aes(x=IEratio, y=fit_half),color="magenta4") + 
  scale_y_log10(lim=c(0.4,60),breaks=c(0.5,1,2,4,8,16,32,64)) + scale_x_log10(limits=c(0.01,100),breaks=c(0.03,0.1,0.3,1,3,10,30),labels=comma) + labs(x="RIME",y="half-life (min)") + 
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),axis.text.y=element_text(size=8))# + geom_vline(xintercept=c(0.70,1.3),color="red")

combo.juncratio.data.parsed$IEratio_bin <- factor(combo.juncratio.data.parsed$IEratio_bin, levels=c("introndef","confused","exondef"))
f3_d_sub <- ggplot(combo.juncratio.data.parsed, aes(x=factor(IEratio_bin),y=fitvalue)) + geom_boxplot(notch=T,aes(color=factor(IEratio_bin)),outlier.color="lightgrey",size=0.75) + 
  annotate("segment",x=c(1.2,2.2),xend=c(1.8,2.8),y=13,yend=13,size=0.25) + annotate("text",x=c(1.5,2.5),y=14,label="***") +
  scale_color_manual(values=c("deeppink4", "snow3","dodgerblue4"),guide=F) +
  scale_x_discrete(labels=c("< 0.75","0.75 - 1.33","> 1.33")) + scale_y_log10(limits=c(0.4,15),breaks=c(1,2,4,8,16)) +
  labs(x="RIME",y="half-life (min)") + 
  theme(panel.background=element_rect(fill="white"),axis.text.x=element_text(size=7,angle=45,hjust=1),axis.text.y=element_text(size=7),axis.title.x=element_text(size=7),axis.title.y=element_text(size=7))
#### F - ID/ED boxplots ####

combo.juncratio.data.parsed.bothdef$IEratio_bin <- factor(combo.juncratio.data.parsed.bothdef$IEratio_bin, levels=c("introndef","exondef"))
combo.juncratio.data.parsed.bothdef$len_bin_intron_def <- factor(combo.juncratio.data.parsed.bothdef$len_bin_intron_def, 
                                                                 levels=c("20%_ID","40%_ID","60%_ID","80%_ID","100%_ID","20%_ED","40%_ED","60%_ED","80%_ED","100%_ED"))
f3_f <- ggplot(combo.juncratio.data.parsed.bothdef, aes(x=factor(IEratio_bin),y=fitvalue,fill=factor(len_bin_intron_def))) + geom_boxplot(notch=T,outlier.color="lightgrey") + 
  scale_y_log10(limits=c(0.4,25),breaks=c(0.5,1,2,4,8,16)) + scale_x_discrete(labels=c("R < 0.75\n(n=20,777)","RIME > 1.33\n(n=3,842)")) +
  scale_fill_manual(values=c(brewer.pal(9,"PuRd")[5:9], brewer.pal(9,"Blues")[5:9]), 
                    labels=c("43-59 nt","60-63 nt","64-68 nt","69-99 nt","99-2,493 nt","57-353 nt","354-674 nt","675-1,252 nt","1,253-2,943 nt","> 2,944 nt"), guide=guide_legend(ncol=5,byrow=T)) +
  labs(x="quintiles of intron length",y="half-life (min)",fill="intron length\n(quintiles)") + 
  theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),axis.title=element_text(size=10),
        legend.position="bottom",legend.key.size=unit(3,"mm"),legend.title=element_text(size=6),legend.text=element_text(size=5))
#### G - SRE enrichment ####
intron.enrich <- subset(all.enrich_rmSS, def=="introndef" & region=="intron")
intron.kmers <- c("UAUUAU","UAUUUA","CUGCUG","UGCUGC")
f3_g_up <- ggplot(intron.enrich, aes(x=log2(enrichment), y=-log2(pval))) + 
  geom_point(alpha=0.5, color="grey27") + geom_point(data=subset(intron.enrich, BH_corrected_pval<10^-35 & abs(log2(enrichment))>0.5), color="magenta4", alpha=0.5) + 
  geom_text_repel(data=subset(intron.enrich, BH_corrected_pval<10^-35 & abs(log2(enrichment)) > 0.5 & kmer %in% intron.kmers), aes(label=kmer), size=2) + 
  facet_grid(def~region, labeller=labeller(def = c("introndef"="RIME < 0.75"))) + xlim(-1.5,1.5) + 
  theme(strip.background=element_rect(fill="deeppink4"),strip.text=element_text(color="white",size=10),
        axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),axis.title=element_text(size=10)) + background_grid(major="only_minor")

exon.enrich <- subset(all.enrich_rmSS, def=="exondef" & (region=="upexon" | region=="downexon"))
exon.kmers <- c("AACAAC","CAACAA","ACAACA","GCAGCA")
f3_g_down <- ggplot(exon.enrich, aes(x=log2(enrichment), y=-log2(pval))) + 
  geom_point(alpha=0.5, color="grey27") + geom_point(data=subset(exon.enrich, BH_corrected_pval<10^-50 & abs(log2(enrichment))>0.5), color="magenta4", alpha=0.5) + 
  geom_text_repel(data=subset(exon.enrich,  BH_corrected_pval<10^-50 & abs(log2(enrichment)) > 0.5 & kmer %in% exon.kmers), aes(label=kmer), size=2) + 
  facet_grid(def~region, labeller=labeller(def = c("exondef"="RIME > 1.33"), region=c("upexon"="upstream exon", "downexon"="downstream exon"))) + xlim(-1.5,1.5) + 
  theme(strip.background=element_rect(fill="dodgerblue4"),strip.text=element_text(color="white",size=10),
        axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),axis.title=element_text(size=10)) + background_grid(major="only_minor")

both.enrich <- rbind(intron.enrich, exon.enrich)
both.kmers <- c(intron.kmers, exon.kmers)
f3_g <- ggplot(both.enrich, aes(x=log2(enrichment), y=-log2(pval))) + 
  geom_point(alpha=0.5, color="grey27") + geom_point(data=subset(both.enrich, BH_corrected_pval<10^-35 & abs(log2(enrichment))>0.5), color="magenta4", alpha=0.5) + 
  geom_text_repel(data=subset(both.enrich,  BH_corrected_pval<10^-35 & abs(log2(enrichment)) > 0.5 & kmer %in% both.kmers), aes(label=kmer), size=2) + 
  facet_grid(~region, labeller=labeller(region=c("upexon"="upstream exon\n(in RIME > 1.33)", "intron"="intron\n(in RIME < 0.75)","downexon"="downstream exon\n(in RIME > 1.33)"))) + xlim(-1.5,1.5) + 
  theme(strip.background=element_rect(fill="dodgerblue4"),strip.text=element_text(color="white",size=8),
        axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),axis.title=element_text(size=10)) + background_grid(major="only_minor")
#### H - accuracy results ####

f3_h <- ggplot(lenmatched, aes(x=factor(type),color=factor(type),y=accuracy*100)) + geom_boxplot(notch=T,fill=NA,outlier.color="lightgrey",size=0.75) + ylim(0,15) + 
  annotate("segment",x=c(1.2,1.2),xend=c(1.8,2.8),y=c(14, 12),yend=c(14,12)) + annotate("text",x=c(1.5,1.5),y=c(14.5,12.5),label="***") + 
  scale_color_manual(values=c("deeppink4","grey57","dodgerblue4"),guide=F) + scale_x_discrete(labels=c("< 0.75","0.75 - 1.33","> 1.33")) +
  labs(y="% non-canonical junctions",x="RIME",color="RIME score") +
  theme(legend.position="bottom",axis.text.x=element_text(angle=45,hjust=1,size=8),axis.text.y=element_text(size=8),axis.title=element_text(size=10),legend.title=element_text(size=6),legend.text=element_text(size=5))
#### COMBINE ####
pdf("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/fig3.pdf",width=8,height=10.5, useDingbats = F)
ggdraw() + #
  draw_plot(f3_a, 0, 0.8, 0.4, 0.15) + draw_plot(f3_b, 0, 0.66, 0.4, 0.15) +
  draw_plot(f3_c, 0.4, 0.66, 0.6, 0.33) +
  draw_plot(f3_d, 0, 0.33, 0.35, 0.33) + draw_plot(f3_d_sub, 0.32, 0.33, 0.2, 0.33) +
  draw_plot(f3_f, 0.5, 0.33, 0.5, 0.33) +
  draw_plot(f3_g, 0, 0, 0.66, 0.33) + 
  draw_plot(f3_h, 0.66, 0, 0.33, 0.33) +
  draw_plot_label(c("A","B","C","D","E","F","G","H"), c(0, 0, 0.4, 0, 0.32, 0.50, 0, 0.66), c(1, 0.8, 1, 0.66, 0.66, 0.66, 0.33, 0.33))
dev.off()

##### FIGURE 4 - gene-wise #####
#### A - regressions ####

replace.names <- c("intron length","gene expression","3' ss strength","5' ss strength",
                   "A+U%","A+U% in 3' region","A+U% in 5' region","intron position",
                   "first intron length","first intron half-life","enhancer in first intron",
                   "upstream exon length","downstream exon length","enhancer in intron",
                   "","","","")

test.lm.data.ID$def <- "intron"
test.lm.data.ID$relimp_fix <- -test.lm.data.ID$relimp*100
test.lm.data.ID$names_full <- replace.names
test.lm.data.ED$def <- "exon"
test.lm.data.ED$relimp_fix <- test.lm.data.ED$relimp*100
test.lm.data.ED$names_full <- replace.names
test.lm.data.def <- rbind(test.lm.data.ID, test.lm.data.ED)

offset <- 12
f4_a <- ggplot(test.lm.data.def, aes(x=factor(names),color=factor(sign))) +
  geom_linerange(data=subset(test.lm.data.def, def=="intron"), aes(ymin=-offset, ymax=-offset+relimp_fix),size=4.5) +
  geom_linerange(data=subset(test.lm.data.def, def=="exon"), aes(ymin=offset, ymax=offset+relimp_fix),size=4.5) +
  geom_errorbar(data=subset(test.lm.data.def, def=="intron"), aes(ymin=-offset-(relimp_upper*100),ymax=-offset-(relimp_lower*100)),color="grey50",width=0.25) +
  geom_errorbar(data=subset(test.lm.data.def, def=="exon"), aes(ymin=offset+(relimp_lower*100),ymax=offset+(relimp_upper*100)),color="grey50",width=0.25) +
  geom_label(aes(x=factor(names), y=0, label=names_full), 
             inherit.aes=F, size=3.5, label.padding=unit(0.0, "lines"), label.size=0,label.r = unit(0.0, "lines"), fill="white", color="black") +
  scale_y_continuous(limits=c(-36.5-offset, 36.5+offset),breaks=c(seq(-35,0,5)+-offset, seq(0,35,5)+offset),
                     labels=as.character(c(rev(seq(0,35,5)),seq(0,35,5)))) +
  scale_color_manual(values=c(wes_palette("Darjeeling")[5], wes_palette("Darjeeling")[3]),labels=c("shorter half-life","longer half-life")) +
  labs(y="relative importance (%)",color="correlated with") + 
  theme(legend.position=c(0.1,0.5),axis.title.y=element_blank(),axis.text.y=element_blank(),
        axis.line.y=element_blank(),axis.ticks.y=element_blank()) + coord_flip() + background_grid(major="x")
#### B - gene-wise stdev of halflives ####
f4_b <- ggplot(combo.stddev.data.median, aes(x=factor(type),y=mean_stdev,fill=factor(type))) + geom_bar(stat="identity",width=0.5) + 
  geom_errorbar(aes(ymin=mean_stdev-sem, ymax=mean_stdev+sem),width=0.25) +
  annotate("segment",x=1,xend=2,y=3.7,yend=3.7,size=0.25,linetype="longdash") + ylim(0,4) +
  annotate("text",x=1.5,y=3.8,label="***") +
  scale_fill_manual(values=c(wes_palette("Royal1")[1], "grey85"),guide=F) + scale_x_discrete(labels=c("introns from\nsame gene","randomly sampled\nintrons")) +
  background_grid(major="y",minor="y") + labs(y="half-life SD (min)",x="",fill="") + theme(axis.text=element_text(size=7),axis.title.y=element_text(size=9))
#### C - gene-wise intron length ####
f4_c <- ggplot(combo.stddev.data.median, aes(x=factor(type),y=mean_intron_stddev,fill=factor(type))) + geom_bar(stat="identity",width=0.5) + 
  geom_errorbar(aes(ymin=mean_intron_stddev-sem_intron_stddev, ymax=mean_intron_stddev+sem_intron_stddev),width=0.25) +
  annotate("segment",x=1,xend=2,y=900,yend=900,size=0.25,linetype="longdash") + annotate("text",x=1.5,y=910,label="***") +
  scale_fill_manual(values=c(wes_palette("Royal1")[1],"grey85"),guide=F) + 
  scale_x_discrete(labels=c("same\ngene","random")) + ylim(0,950) +
  background_grid(major="y",minor="y") + labs(y="intron length SD (nt)",x="",fill="") + theme(axis.text=element_text(size=7),axis.title.y=element_text(size=9))
#### D - gene-wise stddev ####

f4_d <- ggplot(combo.stddev.data.def, aes(x=factor(def),y=stddev,fill=factor(def),alpha=factor(type),color=factor(type))) + geom_boxplot(notch=T,outlier.color="lightgrey",size=0.75) + 
  scale_fill_manual(values=c("deeppink4","snow3","dodgerblue4"),guide=F) + scale_x_discrete(labels=c("all intron\ndefined","mixed \ndefinition","all exon\ndefined")) +
  scale_alpha_manual(values=c(1,0.25),labels=c("introns from the same gene","randomly sampled introns")) + 
  scale_color_manual(values=c("black","grey60"),labels=c("introns from the same gene","randomly sampled introns")) +
  annotate("segment",x=c(0.87,0.87,1.87),xend=c(1.73,2.73,2.73),y=c(2.6,2.4,2.2),yend=c(2.6,2.4,2.2),size=0.35,linetype="longdash") + annotate("text",x=c(2.35,2.35,1.35),y=c(2.25,2.45,2.65),label="***") + 
  ylim(0,3) + labs(x="",y="half-life SD (min)",alpha="",color="") + 
  theme(legend.position="bottom",axis.title.x=element_blank(),axis.text=element_text(size=7),axis.title.y=element_text(size=9),legend.key.size=unit(3,"mm"),legend.text=element_text(size=7))
  
#### E - gene ontology ####
f4_e <- ggplot(siglevel.chosen, aes(x=factor(description),y=log2(enrich),fill=factor(type))) + geom_bar(stat="identity",position="dodge") + 
  scale_fill_manual(values=rev(c("deeppink4","grey75","grey50","dodgerblue4")),guide=guide_legend(reverse=T)) + labs(x="",y="log2(enrichment over background)",fill="") + coord_flip() + 
  theme(legend.position="bottom",axis.title=element_text(size=10),axis.text=element_text(size=7),legend.key.size=unit(3,"mm"),legend.text=element_text(size=7))
#### COMBINE ####
pdf("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/fig4.pdf", width=8, height=10.5, useDingbats = F)
ggdraw() +
  geom_rect(xmin=c(0.035,0.595),xmax=c(0.41,0.985),ymin=c(0.97,0.97),ymax=c(0.99,0.99),color="grey30",fill="grey30") +
  geom_text(x=c(0.2075, 0.7925),y=c(0.98, 0.98), label=c("intron definition","exon definition"),color="white") + 
  draw_plot(f4_a, 0, 0.66, 1, 0.315) +
  draw_plot(f4_b, 0, 0.33, 0.33, 0.33) + draw_plot(f4_c, 0.33, 0.33, 0.17, 0.33) + draw_plot(f4_d, 0.5, 0.33, 0.5, 0.33) +
  draw_plot(f4_e, 0, 0, 0.5, 0.33) +
  draw_plot_label(c("A","B","C","D","E"), c(0, 0, 0.33, 0.5, 0), c(1, 0.66, 0.66, 0.66, 0.33))  
dev.off()





##### SUPP FIGURE 1 - simulation schematic #####
simschem_img <- readPNG("Figures/revisedfigures/SimulationSchematic.png")
simschem_grob <- rasterGrob(simschem_img, interpolate=T)
fs1 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(simschem_grob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  labs(x="",y="") + theme_classic() + theme(axis.text=element_blank(),axis.ticks=element_blank())
#### COMBINE ####
pdf("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/suppfig1.pdf",width = 8,height=10.5, useDingbats = F)
ggdraw() +
  draw_plot(fs1, -0.05, 0.05, 1.1, 1.1) +
  draw_plot_label(c("A","B","C"), c(0.05, 0.05, 0.05), c(1, 0.66, 0.45))
dev.off()

##### SUPP FIGURE 2 - simulations #####
#### A - full spearman correlations ####
fs2_a <- ggplot(subset(full.sim.cors.data, cor_type=="spearman"), aes(x=factor(type),y=meancor, fill=factor(type),color=factor(type))) + 
  geom_bar(stat="identity",position="dodge",color=NA) + geom_errorbar(aes(ymin=meancor-secor, ymax=meancor+secor),width=0.5) + 
  scale_fill_manual(values=c(brewer.pal(9,"RdPu")[3], brewer.pal(9,"Oranges")[3], brewer.pal(9,"BuPu")[6]),guide=F) +
  scale_color_manual(values=c(brewer.pal(9,"RdPu")[5], brewer.pal(9,"Oranges")[5], brewer.pal(9,"BuPu")[8]),guide=F) +
  scale_x_discrete(labels=c("intron ratios",expression(paste(Psi," decay",sep="")),"junction ratios")) + scale_y_continuous(limits=c(0,1.15),breaks=c(0,0.25,0.5,0.75,1)) +
  labs(x="approach",y=expression(paste("correlation btwn est and simulated t"["1/2"]))) + background_grid(major="y",minor="none") + 
  theme(axis.text.x=element_text(size=8,angle=45,hjust=1),axis.text.y=element_text(size=8),axis.title.x=element_blank(),axis.title.y=element_text(size=10)) 
#### B - full % truth across ####
fs2_b <- ggplot(fullsims.error, aes(x=factor(type),y=abslogpererror_mean, fill=factor(type),color=factor(type))) + 
  geom_bar(stat="identity",position="dodge",color=NA) + geom_text(x=1,y=0.05,label="NA",color=brewer.pal(9,"RdPu")[3],size=3) +
  geom_errorbar(aes(ymin=abslogpererror_mean-abslogpererror_sd, ymax=abslogpererror_mean+abslogpererror_sd),width=0.5) + 
  scale_fill_manual(values=c(brewer.pal(9,"RdPu")[3], brewer.pal(9,"Oranges")[3], brewer.pal(9,"BuPu")[6]),guide=F) +
  scale_color_manual(values=c(brewer.pal(9,"RdPu")[5], brewer.pal(9,"Oranges")[5], brewer.pal(9,"BuPu")[8]),guide=F) +
  scale_x_discrete(labels=c("intron ratios",expression(paste(Psi," decay",sep="")),"junction ratios")) + #scale_y_continuous(limits=c(0,10),breaks=c(-50,0,50,100)) +
  labs(x="approach",y=expression(paste("| log2(est t" ["1/2"], " / sim t" ["1/2"], " ) |"))) +  background_grid(major="y",minor="none") +
  theme(axis.text.x=element_text(size=8,angle=45,hjust=1),axis.text.y=element_text(size=8),axis.title.x=element_blank(),axis.title.y=element_text(size=10))
#### C - full % relative truth ####
fs2_c <- ggplot(fullsims.relative.error, aes(x=factor(type),y=abslogpererror_mean, fill=factor(type),color=factor(type))) + 
  geom_bar(stat="identity",position="dodge",color=NA) + 
  geom_errorbar(aes(ymin=abslogpererror_mean-abslogpererror_se, ymax=abslogpererror_mean+abslogpererror_se),width=0.5) + 
  scale_fill_manual(values=c(brewer.pal(9,"RdPu")[3], brewer.pal(9,"Oranges")[3], brewer.pal(9,"BuPu")[6]),guide=F) +
  scale_color_manual(values=c(brewer.pal(9,"RdPu")[5], brewer.pal(9,"Oranges")[5], brewer.pal(9,"BuPu")[8]),guide=F) +
  scale_x_discrete(labels=c("intron ratios",expression(paste(Psi," decay",sep="")),"junction ratios")) + #scale_y_continuous(limits=c(0,10),breaks=c(-50,0,50,100)) +
  labs(x="approach",y=expression(paste("| log2( est t" ["1/2"]^"j", "/t" ["1/2"]^"k", " / sim t" ["1/2"]^"j", "/t" ["1/2"]^"k", ") |"))) +  background_grid(major="y",minor="none") +
  theme(axis.text.x=element_text(size=8,angle=45,hjust=1),axis.text.y=element_text(size=8),axis.title.x=element_blank(),axis.title.y=element_text(size=9))
#### D - full dist of all methods, exp sorted ####

fs2_d1 <- ggplot(subset(fullsims, type=="ratio"), aes(x=half_life,y=sim_hl,color=factor(expression_level))) + geom_point(size=1.5) + stat_smooth(size=2) + xlim(0,60) +
  scale_color_manual(values=c(brewer_pal(palette="YlOrRd")(9),"black"),guide=F) + geom_abline(color="yellow",linetype="longdash",size=1.5) + 
  labs(x=" ",y="estimated half-life (min)",color="expression level") + background_grid(major="xy",minor="none") + 
  theme(legend.position="bottom",axis.text=element_text(size=10),axis.title=element_text(size=12)) + ggtitle("intron ratio")
fs2_d2 <- ggplot(subset(fullsims, type=="psi"), aes(x=half_life,y=sim_hl,color=factor(expression_level))) + geom_point(size=1.5) + stat_smooth(size=2) + xlim(0,60) +
  scale_color_manual(values=c(brewer_pal(palette="YlOrRd")(9),"black"),guide=F) + geom_abline(color="yellow",linetype="longdash",size=1.5) + 
  labs(x="simulated half-life (min)",y="",color="expression level") + background_grid(major="xy",minor="none") + 
  theme(legend.position="bottom",axis.text=element_text(size=10),axis.title=element_text(size=12)) + ggtitle(expression(paste(Psi," decrease",sep="")))
fs2_d3 <- ggplot(subset(fullsims, type=="junc"), aes(x=half_life,y=sim_hl,color=factor(expression_level))) + geom_point(size=1.5) + stat_smooth(size=2) + xlim(0,60) +
  scale_color_manual(values=c(brewer_pal(palette="YlOrRd")(9),"black"),guide=F) + geom_abline(color="yellow",linetype="longdash",size=1.5) + 
  labs(x=" ",y="",color="expression level") + background_grid(major="xy",minor="none") + 
  theme(legend.position="bottom",axis.text=element_text(size=10),axis.title=element_text(size=12)) + ggtitle("junction dynamics")

legend.simexp <- ggplot(subset(fullsims, type=="psi"), aes(x=half_life,y=sim_hl,color=factor(expression_level))) + geom_point(size=2) + stat_smooth(size=2,fill="lightgrey") + 
  scale_color_manual(values=c(brewer_pal(palette="YlOrRd")(9),"black"),guide=guide_legend(nrow=1)) + 
  labs(x="simulated half-life (min)",y="",color="expression level") +theme(legend.position="bottom")
fs2_legend <- g.legend(legend.simexp)

fs2_d <- ggdraw() + 
#  geom_rect(xmin=0, xmax=0.33, ymin=0.91, ymax=0.99, color="white",fill="grey30",size=3) +
#  geom_rect(xmin=c(0.33, 0.66), xmax=c(0.66, 1), ymin=0.91, ymax=0.99, color="white", fill="grey30",size=3) +
#  geom_text(x=0.165,y=0.95,label=c("intron ratio"),color="white") +
#  geom_text(x=c(0.495, 0.825), y=0.95, label=c("relative~Psi","junction~ratio"),parse=T,color="white") + 
  draw_plot(fs2_d1, 0, 0.05, 0.33, 0.95) + draw_plot(fs2_d2, 0.33, 0.05, 0.33, 0.95) + draw_plot(fs2_d3, 0.66, 0.05, 0.33, 0.95) +
  draw_plot(fs2_legend, 0.25, 0, 0.5, 0.1)
#### NIX - full directional % truth ####
#fs1_e <- ggplot(fullsims.relative.error.hl20, aes(x=factor(type),y=logpererror_mean, fill=factor(type),color=factor(type))) + 
#  geom_bar(stat="identity",position="dodge",color=NA) + 
#  geom_errorbar(aes(ymin=logpererror_mean-logpererror_se, ymax=logpererror_mean+logpererror_se),width=0.5) + 
#  scale_fill_manual(values=c(brewer.pal(9,"RdPu")[3], brewer.pal(9,"Oranges")[3], brewer.pal(9,"BuPu")[6]),guide=F) +
#  scale_color_manual(values=c(brewer.pal(9,"RdPu")[5], brewer.pal(9,"Oranges")[5], brewer.pal(9,"BuPu")[8]),guide=F) +
#  scale_x_discrete(labels=c("intron ratios",expression(paste("relative ",Psi,sep="")),"junction ratios")) + #scale_y_continuous(limits=c(0,10),breaks=c(-50,0,50,100)) +
#  labs(x="approach",y="log2( relative estimate / relative truth )") +  background_grid(major="y",minor="none") +
#  theme(axis.text.x=element_text(size=8,angle=45,hjust=1),axis.text.y=element_text(size=8),axis.title.x=element_blank(),axis.title.y=element_text(size=10))
#### E - intron length divides ####
fullsims.cor$type <- factor(fullsims.cor$type, levels=c("ratio","psi","junc"))
fs2_e <- ggplot(subset(fullsims.cor), aes(x=factor(intron),y=relate,fill=factor(type))) + 
  geom_boxplot(data=subset(fullsims.cor, expression_level >= 5 & cor_type=="ME" & type=="psi"), color="ivory",outlier.color=NA) + 
  geom_boxplot(data=subset(fullsims.cor, expression_level >= 5 & cor_type=="ME" & type=="junc"), color="ivory",outlier.color=NA) + 
  scale_y_log10(limits=c(1,500),breaks=c(0.3,1,3,10,30,100,300),label=comma) + scale_x_discrete(breaks=c(40,100,500,1000,11000,21000,31000,41000)) +
  scale_fill_manual(values=c(brewer.pal(9,"BuPu")[8], brewer.pal(9,"Oranges")[5]),labels=rev(c("junction dynamics", expression(paste(Psi," decrease",sep=""))))) +
  labs(x="intron length (nt)",y="mean error (min)",fill="approach") + background_grid(major="y",minor="none") + 
  theme(legend.position=c(0.75,0.75),axis.text.y=element_text(size=8),axis.text.x=element_text(size=7,angle=45,hjust=1),axis.title=element_text(size=10),
        legend.text.align=0,legend.key.size=unit(3,"mm"),legend.title=element_text(size=10),legend.text=element_text(size=10))
  
#### F - D_dist divides ####
fs2_f <- ggplot(subset(fullsims.cor), aes(x=factor(D_dist*1000),y=relate,fill=factor(type))) + 
  geom_boxplot(data=subset(fullsims.cor, expression_level >= 5 & cor_type=="ME" & type=="psi"), color="ivory",outlier.color=NA,notch=T) + 
  geom_boxplot(data=subset(fullsims.cor, expression_level >= 5 & cor_type=="ME" & type=="junc"), color="ivory",outlier.color=NA,notch=T) + 
  scale_y_log10(limits=c(1,500),breaks=c(0.3,1,3,10,30,100,300),label=comma) + 
  scale_fill_manual(values=c(brewer.pal(9,"BuPu")[8], brewer.pal(9,"Oranges")[5]),labels=rev(c("junction dynamics", expression(paste(Psi," decrease",sep="")))),guide=guide_legend(nrow=1)) +
  labs(x="downstream distance (nt)",y="mean error (min)",fill="approach") + background_grid(major="y",minor="none") + 
  theme(legend.position=c(0.6,0.9),axis.text.y=element_text(size=8),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title=element_text(size=10),
        legend.text.align=0,legend.key.size=unit(3,"mm"),legend.title=element_text(size=10),legend.text=element_text(size=10))
#### COMBINE ####

pdf("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/suppfig2.pdf",width=8,height=10.5,useDingbats = F)
ggdraw() +
  draw_plot(fs2_a, 0, 0.66, 0.33, 0.33) + draw_plot(fs2_b, 0.33, 0.66, 0.33, 0.33) + draw_plot(fs2_c, 0.66, 0.66, 0.33, 0.33) +
  draw_plot(fs2_d, 0, 0.33, 1, 0.33) +
  draw_plot(fs2_e, 0, 0, 0.5, 0.33) + draw_plot(fs2_f, 0.5, 0, 0.5, 0.33) +
  draw_plot_label(c("A","B","C","D","E","F"), c(0, 0.33, 0.66, 0, 0, 0.5), c(1, 1, 1, 0.66, 0.33, 0.33))
dev.off()
  

##### SUPP FIGURE 3 - fit confidence #####
#### A - ratios by timepoint/replicate ####
fs3_a <- ggplot(juncratio.data, aes(x=factor(time),y=ratio,fill=factor(rep))) + geom_boxplot(notch=T) + scale_fill_manual(values=c("dodgerblue3","dodgerblue2","dodgerblue1")) +
  scale_y_log10(limits=c(0.05,7.5),breaks=c(0.1, 0.25, 0.75, 1, 2.5),labels=comma) + labs(x="labeling period",y="ratio of IE/EE\njunction reads",fill="replicates") +
  theme(legend.position="bottom",axis.text=element_text(size=8))
#### B - top1000 genes fits ####
fs3_b <- ggplot(coef.data.exp, aes(x=factor(type),y=coef, fill=factor(time))) + geom_boxplot(notch=T) + ylim(0,1) +
  scale_fill_manual(values=c(brewer.pal(11,"BrBG")[8:10], rev(brewer.pal(11,"BrBG")[2:4]), "darkgrey")) + scale_x_discrete(labels=c("within time","across time","across all")) +
  labs(x="comparison",y="coefficient of variation",fill="labeling periods") + theme(legend.position="bottom",legend.key.size=unit(3,"mm"),legend.title=element_text(size=10),axis.text=element_text(size=8))
#### C - confidence intervals for fits ####
fs3_c <- ggplot(combo.juncratio.data.parsed, aes(x=fitvalue,y=se)) + geom_point(alpha=0.15,size=0.75,color="grey25") + 
  scale_x_log10(limit=c(0.25,100),breaks=c(0.5,1,2,4,8,16,32,64)) + scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100),labels=comma) +
  labs(x="half-life (min)",y="bootstrap SE (min)") + theme(axis.text.x=element_text(angle=45,hjust=1,size=8),axis.text.y=element_text(size=8)) + background_grid(major="xy")
#### D - residual sum of squares ####
fs3_d <- ggplot(data.frame(rsq = rsq.zero), aes(x=log10(rsq))) + geom_histogram(fill="black",color="white") + xlim(-25,2) +
  labs(x="log10(residual sum of squares)",y="count") +
  theme(axis.text=element_text(size=8)) + background_grid(major="xy")
#### E - variability in half-lives with different transcription rates ####
# read line is median of 1500 values
fs3_e <- ggplot(sumsqfit.data.all, aes(x=factor(txnrate), y=halflife)) +  
  geom_boxplot(notch=T) + scale_y_log10(limits=c(0.1,100),breaks=c(0.5,2,8,1e2,1e3),label=comma) + 
  geom_hline(yintercept=median(combo.juncratio.data.parsed$fitvalue),color="red",linetype="dashed") +
  labs(x="transcription rate (nt/min)",y="half-life (min)") + theme(axis.text=element_text(size=8))
#### COMBINE ####
pdf("~/Desktop/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/suppfig3.pdf",width=8,height=10.5,useDingbats = F)
ggdraw() +
  draw_plot(fs3_a, 0, 0.66, 0.5, 0.33) +
  draw_plot(fs3_b, 0.5, 0.66, 0.5, 0.33) +
  draw_plot(fs3_c, 0, 0.33, 0.5, 0.33) +
  draw_plot(fs3_d, 0.5, 0.33, 0.5, 0.33) +
  draw_plot(fs3_e, 0, 0, 0.5, 0.33) +
#  geom_rect(xmin=0.25,xmax=0.75,ymin=0.33,ymax=0.66) +
#  geom_text(x=0.5,y=0.495,label="standard errors",color="white") + 
  draw_plot_label(c("A","B","C","D","E"), c(0, 0.5, 0, 0.5, 0), c(1, 1, 0.66, 0.66, 0.33))
dev.off()

##### SUPP FIGURE 4  - 60nt #####
#### A - intron length distribution ####
fs4_a <- ggplot(combo.juncratio.data.parsed, aes(x=intronlen)) + geom_histogram(color="white") + scale_x_log10(breaks=c(50,100,1000,10000,50000),labels=comma) + 
  scale_y_continuous(limits=c(0,15000),breaks=c(5000,10000,15000),labels=comma) + labs(x="intron length (nt)") + background_grid(major="xy") +
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))
#### B - splice sites for short introns ####
fs4_b <- ggplot(subset(miso.splicesite.data, len <=100), aes(x=factor(len.bin5),y=SS,fill=factor(type))) + geom_boxplot(notch=T) + 
  scale_fill_manual(values=c("darkorange2","deepskyblue"),labels=c("3'ss","5'ss"),guide=guide_legend(nrow=1)) + ylim(-3,15) + 
  labs(x="intron length (nt)",y="maxEnt score",fill="splice site") + theme_cowplot() + 
  background_grid(major="y",minor="y") + theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom")
#### C - matching splice sites to explain slow splicing of short introns ####
t.test(subset(ss.match.data, type=="40-50 nt")$half.lives, subset(ss.match.data, type=="60-70 nt\n(matched for both ss scores)")$half.lives)

fs4_c <- ggplot(ss.match.data, aes(x=factor(type),y=half.lives,fill=factor(type))) + geom_boxplot(notch=T,outlier.color="white") + scale_y_log10(limits=c(0.25,10),breaks=c(0.5,1,2,4)) + 
  #geom_segment(aes(x=1.25,xend=2.75,y=60,yend=60),linetype="dashed",size=0.5,color="darkgrey") + geom_text(aes(x=2,y=65,label="***"),size=3) + 
  #geom_segment(aes(x=1.25,xend=3.75,y=80,yend=80),linetype="dashed",size=0.5,color="darkgrey") + geom_text(aes(x=2.5,y=85,label="***"),size=3) + 
  scale_x_discrete(labels=c("40-50 nt","60-70 nt\n(3'ss)","60-70 nt\n(5'ss)","60-70 nt\n(both)")) + guides(fill=guide_legend(nrow=2)) + 
  scale_fill_manual(values=c("royalblue1","lightblue1","lightblue1","lightblue1"),breaks=c("40-50 nt","60-70 nt\n(matched for 3'ss score)"),labels=c("actual distribution","matched for splice site score")) + 
  labs(x="intron length",y="half-life (min)",fill=NULL) + background_grid(major="y",minor="y") + 
  theme(legend.position="bottom",axis.text.x=element_text(angle=45,hjust=1,size=8),legend.text=element_text(size=9))
#### D - ss distribution for all introns ####
fs4_d <- ggplot(miso.splicesite.data, aes(x=factor(len_bin),y=SS,fill=factor(type))) + geom_boxplot(notch=T) + 
  scale_fill_manual(values=c("darkorange2","deepskyblue"),labels=c("3'ss","5'ss"),guide=F) + ylim(0,15) + 
  scale_x_discrete(labels=c("40-60nt","61-65nt","66-77nt","78-284nt",">284nt")) +
  labs(x="intron length (quintiles)",y="maxEnt score",fill="splice site") + theme_cowplot() + 
  background_grid(major="y",minor="y") + theme(axis.text.x=element_text(angle=45,hjust=1))
#### NIX?? E - distribution of introns by regulation type ####
fs4_e <- ggplot(subset(combo.juncratio.data.parsed, type!="SEcontaining"), aes(x=intronlen, y=..density..,fill=factor(type))) + geom_histogram(position="dodge",color="white") + 
  scale_x_log10(limits=c(30,50000),breaks=c(50,100,1000,10000),labels=comma) + scale_fill_manual(values=wes_palette("Darjeeling")[c(5,4,3)]) + 
  labs(x="intron length",y="density",fill="intron type") + theme(legend.position=c(0.6,0.5))
#### COMBINE ####
pdf("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/suppfig4.pdf",width=8,height=10.5, useDingbats = F)
ggdraw() + 
  draw_plot(fs4_a, 0, 0.66, 1, 0.33) + 
  draw_plot(fs4_b, 0, 0.33, 0.5, 0.33) + draw_plot(fs4_c, 0.5, 0.33, 0.5, 0.33) + 
  draw_plot(fs4_d, 0, 0, 0.5, 0.33) + #draw_plot(fs4_e, 0.5, 0, 0.5, 0.33) +
  draw_plot_label(c("A","B","C","D"), c(0, 0, 0.5, 0), c(1, 0.66, 0.66, 0.33))
dev.off()

##### SUPP FIGURE 5 - RIME #####
#### NIX - up vs. down exon ####
#fs5_a <- ggplot(combo.juncratio.data.parsed, aes(x=upexon_len, y=downexon_len)) + geom_point(alpha=0.25, color="grey25") + geom_smooth(method="lm") +
#  scale_x_log10(limits=c(10,10000),breaks=c(10,100,1000,10000),labels=comma) + scale_y_log10(limits=c(10,10000),breaks=c(10,100,1000,10000),labels=comma) + 
#  annotate("text",x=20,y=10,label="R = 0.05",color="blue") +
#  labs(x="upstream exon length (nt)",y="downstream exon length (nt)") + background_grid(major="xy") + 
#  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))
#### A - binning scheme for aggregate I&E ####
fs5_a <- ggplot(circ.data, aes(x=intronlen, y=meanexon, color=factor(circ_bin))) + geom_point() + 
  scale_x_log10(limits=c(40,10000),breaks=c(10,100,1000,10000,100000),labels=comma) + scale_y_log10(limits=c(40,10000),breaks=c(10,100,1000,10000),labels=comma) + 
  scale_color_manual(values=c(brewer_pal(palette="YlGn")(9),"black"),labels=c("10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")) + 
  labs(x="intron length (nt)", y="mean exon length (nt)", color="joint deciles of intron \nand mean exon length") + background_grid(major="xy") +
  theme(legend.position="bottom",axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=8),legend.title=element_text(size=10),legend.key.size=unit(3,"mm"))
#### B - binning scheme for RIME ####
fs5_b <- ggplot(circ.data, aes(x=intronlen, y=meanexon, color=RIMEbin)) + geom_point(alpha=0.5) + 
  scale_x_log10(limits=c(40,10000),breaks=c(10,100,1000,10000,100000),labels=comma) + scale_y_log10(limits=c(40,10000),breaks=c(10,100,1000,10000),labels=comma) + 
  scale_color_gradient2(low="darkmagenta",mid="yellow",high="darkblue") + labs(x="intron length (nt)", y="mean exon length (nt)", color="percentiles of RIME\n(40 bins)") + background_grid(major="xy") +
  theme(legend.position="bottom",axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(angle=45,hjust=1,size=8),legend.title=element_text(size=10),legend.key.size=unit(3,"mm"))
#### C - RIME vs. intron/exon length ####
fs5_c <- ggplot(RIME.heatmap.data, aes(y=factor(RIME), x=factor(unit), fill=hl_bin)) + geom_tile() + 
  scale_fill_gradient(low="white",high="black",na.value="lightblue") + geom_hline(yintercept=32.5,color="yellow") +
  labs(x="percentiles of intron | exon length", y="percentiles of RIME",fill="bins of half-life") + facet_wrap(~type, nrow=1) + 
  theme(legend.position="bottom", axis.text.y=element_text(size=7),axis.text.x=element_text(size=7,angle=45,hjust=1),axis.title=element_text(size=10),
        legend.text=element_text(size=7,angle=45,hjust=1),legend.title=element_text(size=8),legend.key.size=unit(4,"mm"),
        strip.background=element_rect(fill="grey17"),strip.text=element_text(color="white",size=10))
#### D - highlighting stripes
stripeplot_img <- readPNG("~/Desktop/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/circplot-06.png")
stripeplot_grob <- rasterGrob(stripeplot_img, interpolate=T)
fs5_d <- qplot(1:10, 1:10, geom="blank") + 
  annotation_custom(stripeplot_grob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  labs(x="",y="") + theme_classic() + theme(axis.text=element_blank(), axis.ticks=element_blank())
#### E - intron length in bins ####
fs5_e <- ggplot(circ.intron.data, aes(x=factor(RIME),y=factor(circ),fill=hl_bin,label=round(mean_intron,0))) + geom_tile(alpha=0.5) + 
  scale_fill_gradient(low="yellow",high="magenta4",na.value="black") + 
  geom_text(data=subset(circ.intron.data, RIME < -0.8),size=1.5) + geom_text(data=subset(circ.intron.data, circ <= 0.8 & RIME <= -0.75),size=1.5) +
  geom_text(data=subset(circ.intron.data, circ <= 0.7 & RIME < -0.68),size=1.5) + geom_text(data=subset(circ.intron.data, circ <= 0.6 & RIME <= -0.625),size=1.5) +
  geom_text(data=subset(circ.intron.data, circ <= 0.5 & RIME < -0.56),size=1.5) + geom_text(data=subset(circ.intron.data, circ <= 0.4 & RIME <= -0.5),size=1.5) +
  geom_text(data=subset(circ.intron.data, circ <= 0.3 & RIME < -0.40),size=1.5) + geom_text(data=subset(circ.intron.data, circ <= 0.2),size=1.5) +
  labs(x="percentiles of RIME",y="percentiles of joint \nintron & mean exon length",fill="half-life") +
  theme(legend.position="bottom",axis.text.x=element_text(size=7,angle=45,hjust=1),axis.text.y=element_text(size=7),axis.title=element_text(size=9),
        legend.text=element_text(angle=45,hjust=1,size=8),legend.title=element_text(size=9),legend.key.size=unit(3,"mm"))
#### F - exon binning ####
combo.juncratio.data.parsed.bothdef$IEratio_bin <- factor(combo.juncratio.data.parsed.bothdef$IEratio_bin, levels=c("introndef","exondef"))
combo.juncratio.data.parsed.bothdef$len_bin_exon_def <- factor(combo.juncratio.data.parsed.bothdef$len_bin_exon_def, 
                                                               levels=c("20%_ID","40%_ID","60%_ID","80%_ID","100%_ID","20%_ED","40%_ED","60%_ED","80%_ED","100%_ED"))
fs5_f <- ggplot(combo.juncratio.data.parsed.bothdef, aes(x=factor(IEratio_bin),y=fitvalue,fill=factor(len_bin_exon_def))) + geom_boxplot(notch=T,outlier.color="lightgrey") + 
  scale_y_log10(limits=c(0.4,25),breaks=c(1,2,4,8,16)) + scale_x_discrete(labels=c("intron definition (n = 20,277)","exon definition (n = 3,842)")) +
  scale_fill_manual(values=c(brewer.pal(9,"PuRd")[5:9], brewer.pal(9,"Blues")[5:9]), 
                    labels=c("77-217 nt","218-342 nt","343-514 nt","515-799 nt","780-9,415 nt","34-155 nt","156-228 nt","289-334 nt","335-543 nt","544-4,300 nt"), guide=guide_legend(nrow=4,byrow=T)) +
  labs(x="",y="half-life (min)",fill="mean exon length\n(quintiles)") + 
  theme(legend.position="bottom",axis.title.x=element_blank(),axis.title.y=element_text(size=10),axis.text=element_text(size=8),
        legend.text=element_text(size=8),legend.title=element_text(size=10),legend.key.size=unit(3,"mm"))
#### COMBINE ####
pdf("~/Desktop/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/suppfig5.pdf",width=8,height=10.5, useDingbats = F)
ggdraw() + #
  draw_plot(fs5_a, 0, 0.66, 0.5, 0.33) + draw_plot(fs5_d, 0.55, 0.7, 0.45, 0.3) + 
  draw_plot(fs5_b, 0, 0.33, 0.5, 0.33) + draw_plot(fs5_e, 0.5, 0.33, 0.5, 0.33) + 
  draw_plot(fs5_c, 0, 0, 0.5, 0.33) + draw_plot(fs5_f, 0.5, 0, 0.5, 0.33) +
  draw_plot_label(c("A","D","B","E","C","F"), c(0, 0.5, 0, 0.5, 0, 0.5), c(1, 1, 0.66, 0.66, 0.33, 0.33))
dev.off()

#### FOR BINNING SCHEMATIC ####
pdf("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/binningscheme.pdf",width=8,height=8,useDingbats = F)
ggplot(circ.data, aes(x=intronlen, y=meanexon, color=factor(circ_bin))) + geom_point() + 
  scale_x_log10(limits=c(40,10000),breaks=c(10,100,1000,10000,100000),labels=comma) + scale_y_log10(limits=c(40,10000),breaks=c(10,100,1000,10000),labels=comma) + 
  scale_color_manual(values=c(rep(c("grey","black"),5)),labels=c("10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"),guide=F) + 
  labs(x="intron length (nt)", y="mean exon length (nt)", color="joint deciles of intron \nand mean exon length") + background_grid(major="xy") +
  theme(legend.position="bottom",axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=8),legend.title=element_text(size=10),legend.key.size=unit(3,"mm"))
ggplot(circ.data, aes(x=intronlen, y=meanexon, color=factor(RIMEbin))) + geom_point(alpha=0.5) + 
  scale_x_log10(limits=c(40,10000),breaks=c(10,100,1000,10000,100000),labels=comma) + scale_y_log10(limits=c(40,10000),breaks=c(10,100,1000,10000),labels=comma) + 
  scale_color_manual(values=c(rep(c("grey","black"),20)),guide=F) + labs(x="intron length (nt)", y="mean exon length (nt)", color="percentiles of RIME\n(40 bins)") + background_grid(major="xy") +
  theme(legend.position="bottom",axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(angle=45,hjust=1,size=8),legend.title=element_text(size=10),legend.key.size=unit(3,"mm"))
dev.off()


##### SUPP FIGURE 6 - full kmers #####
#### A - full kmer plots ####
all.enrich_rmSS$def <- factor(all.enrich_rmSS$def, levels=c("introndef","exondef"))
all.enrich_rmSS$region <- factor(all.enrich_rmSS$region, levels=c("upexon","intron","downexon"))

fs6_a <- ggplot(all.enrich_rmSS, aes(x=log2(enrichment), y=-log10(pval))) + geom_point(alpha=0.5) + 
  geom_text_repel(data=subset(all.enrich_rmSS, BH_corrected_pval<10^-35 & abs(log2(enrichment)) > 0.5), aes(label=kmer),size=2) +
  facet_grid(def~region, labeller=labeller(def=c("introndef"="RIME < 0.75","exondef"="RIME > 1.33"), region=c("upexon"="upstream exon", "intron"="intron","downexon"="downstream exon"))) +
  theme(strip.background=element_rect(fill="grey17"),strip.text=element_text(color="white",size=10),
        axis.text=element_text(size=8),axis.title=element_text(size=10)) + background_grid(major="xy")
#### B - kmer enrichments ####
fs6_b <- ggplot(upexon.enrich, aes(x=factor(quant),y=n.kmers/len, fill=factor(quant))) + geom_boxplot(notch=T) + ylim(0,0.05) +
  scale_fill_manual(values=c(brewer.pal(9,"PuRd")[5:9], "snow3", brewer.pal(9,"Blues")[5:9]),guide=F) + 
  scale_x_discrete(labels=c("43-59 nt","60-63 nt","64-68 nt","69-99 nt","99-2,493 nt","RIME ~ 1","57-353 nt","354-674 nt","675-1,252 nt","1,253-2,943 nt","> 2,944 nt")) + 
  labs(x="intron bin",y="density of kmers") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=6),axis.text.y=element_text(size=6),axis.title=element_text(size=8))
#### COMBINE ####
pdf("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/suppfig6.pdf",width=8,height=10.5,useDingbats = F)
ggdraw() +
  draw_plot(fs6_a, 0, 0.5, 1, 0.5) +
  draw_plot(fs6_b, 0.25, 0.2, 0.5, 0.3) +
  draw_plot_label(c("A","B"), c(0, 0.25), c(1, 0.5))
dev.off()

#### SUPP FIGURE 7 - regression coefficients ####
#### A - ID coefficients #### 
fs7_a <- ggplot(test.lm.data.ID, aes(x=factor(names),y=estimate)) + geom_hline(yintercept=0,linetype="dashed") +
  geom_errorbar(aes(ymin=estimate-err,ymax=estimate+err),width=0.5,color=wes_palette("Rushmore")[4]) + geom_point(aes(size=-log10(pval)),color=wes_palette("Royal1")[1]) + ylim(-0.35,0.35) +
  scale_x_discrete(labels=rev(c("intron length","intron position","enhancer in intron","","5'ss strength","3'ss strength","",
                                "upstream exon length","downstream exon length","","A+U%","A+U% in 3' region","A+U% in 5' region","",
                                "first intron length","first intron half-life","enhancer in first intron","gene expression"))) +
  labs(y="regression coefficient",x="",size="-log10(regression p-value)") + coord_flip() + background_grid(major="x") + ggtitle("Intron-defined introns") + 
  theme(legend.position="bottom",axis.title.y=element_blank(),axis.text.y=element_text(size=8),legend.text=element_text(size=8),legend.title=element_text(size=8))
#### B - ED coefficients ####
fs7_b <- ggplot(test.lm.data.ED, aes(x=factor(names),y=estimate)) + geom_hline(yintercept=0,linetype="dashed") +
  geom_errorbar(aes(ymin=estimate-err,ymax=estimate+err),width=0.5,color=wes_palette("Rushmore")[4]) + geom_point(aes(size=-log10(pval)),color=wes_palette("Royal1")[1]) + ylim(-0.35,0.35) +
  scale_x_discrete(labels=rev(c("intron length","intron position","enhancer in intron","","5'ss strength","3'ss strength","",
                                "upstream exon length","downstream exon length","","A+U%","A+U% in 3' region","A+U% in 5' region","",
                                "first intron length","first intron half-life","enhancer in first intron","gene expression"))) +
  labs(y="regression coefficient",x="",size="-log10(regression p-value)") + coord_flip() + background_grid(major="x") + ggtitle("Exon-defined introns") +
  theme(legend.position="bottom",axis.title.y=element_blank(),axis.text.y=element_text(size=8),legend.text=element_text(size=8),legend.title=element_text(size=8))
#### C - gene expression for ID/ED genes ####

fs7_c <- ggplot(combo.metadata.polyA.half.median, aes(x=factor(definition),y=TPM,color=factor(definition))) + geom_boxplot(notch=T,outlier.color = "lightgrey",size=0.75) + scale_y_log10(limits=c(5,5000)) + 
  annotate("segment",x=c(1.2, 1.2, 2.2),xend=c(2.8, 1.8, 2.8),y=c(500, 700, 700),yend=c(500, 700, 700)) + annotate("text",x=c(2,1.5,2.5),y=c(600,800,800),label=c("***","***","***")) +
  scale_color_manual(values=c("deeppink4","snow3","dodgerblue4"),guide=F) + 
  scale_x_discrete(labels=c("all intron\ndefined","mixed\ndefinition","exon\ndefined")) +
  labs(x="",y="gene expression (TPM)") + theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10),axis.text=element_text(size=7))
  
#### D - SD when varying txn rates ####
fs7_d <- ggplot(combo.stddev.data.median.txnrate, aes(x=factor(type),y=mean_stdev,fill=factor(type))) + geom_bar(stat="identity",width=0.5) + geom_errorbar(aes(ymin=mean_stdev-sem, ymax=mean_stdev+sem),width=0.25) +
  ylim(0,6.5) + annotate("segment",x=c(1,2),xend=c(3,3),y=c(6.2,5.2),yend=c(6.2, 5.2),size=0.25) +
  annotate("text",x=c(2.5,2.5),y=c(6.3,5.2),label="***") +
  scale_fill_manual(values=c(wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[4], wes_palette("Royal1")[1]),guide=F) + scale_x_discrete(labels=c("introns from\nsame gene","introns from same gene \nw/ diff txn rates","randomly sampled\nintrons")) +
  background_grid(major="y",minor="y") + labs(y="half-life SD (min)",fill=NA) + 
  theme(axis.text=element_text(size=7),axis.title.x=element_blank(),axis.title.y=element_text(size=10))

#### E - coef for same gene/random introns ####
setnum = 10
fs7_e <- ggplot(combo.metadata.polyA.half.mean) + stat_ecdf(data=subset(sampling.data.combo.half.mean, type=="all"), aes(x=ActualDist.stddev/ActualDist.splicing.rate,color=factor(both))) +
  stat_ecdf(data=subset(sampling.data.combo.half.mean, type=="remaining"), aes(x=FirstAlways.stddev/FirstAlways.splicing.rate,color=factor(both))) +
  stat_ecdf(aes(x=Std.Dev.all.introns/Splicing.rate.all.introns,color='blue'),size=1.5) + ylim(0,1) + 
  stat_ecdf(aes(x=Std.Dev.remaining.introns/Splicing.rate.remaining.introns,color='lightblue'),size=1.5) + 
  scale_color_manual(breaks=c("all.1","remaining.1","blue","lightblue"),values=c(rep(wes_palette("Moonrise1")[4],setnum),wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[4],rep(wes_palette("Royal1")[1],setnum)),guide=guide_legend(nrow=4),
                     labels=c("randomly sampled introns","randomly sampled introns (no first intron)","introns from same gene","introns from same gene (no first intron)")) + theme_bw() + 
  labs(x="SD/mean of half-lives",y="cumulative density",color="") + 
  theme(legend.position=c(0.35,0.85),legend.key=element_blank(),legend.text=element_text(size=7),legend.margin=unit(0.5,"mm"),legend.key.size=unit(1,"mm"),axis.text=element_text(size=8),axis.title=element_text(size=10)) +
  scale_x_log10(limits=c(0.05,2),breaks=c(0.1,0.5,1),labels=comma)
#### NIX? E - alt splicing enrichment for ID/ED genes ####
#fs7_e <- ggplot(combo.metadata.polyA.half.median, aes(x=factor(definition),y=skipped_exon_count/all_exon_count,color=factor(definition))) + geom_boxplot(notch=T,outlier.color="lightgrey",size=0.75) +
#  scale_color_manual(values=c("deeppink4","snow3","dodgerblue4"),guide=F) + scale_x_discrete(labels=c("all intron defined","mixed definition","exon defined")) +
#  ylim(0,0.75) + labs(x="",y="proportion alternative exons") + theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10),axis.text=element_text(size=8))
#### COMBINE ####
pdf("~/Desktop/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/suppfig7.pdf",width=8,height=10.5, useDingbats = F)
ggdraw() + 
  draw_plot(fs7_a, 0, 0.5, 0.5, 0.5) + draw_plot(fs7_b, 0.5, 0.5, 0.5, 0.5) + 
  draw_plot(fs7_c, 0, 0.25, 0.5, 0.25) + draw_plot(fs7_d, 0.5, 0.25, 0.5, 0.25) + 
  draw_plot(fs7_e, 0, 0, 0.5, 0.25) +
  draw_plot_label(c("A","B","C","D","E"), c(0, 0.5, 0, 0.5, 0), c(1, 1, 0.5, 0.5, 0.25))
dev.off()

#### SUPP FIGURE 8 - first intron ####
#### An - 1st intron longer/slower ####
fs8_a <- ggplot(pointrange.data, aes(x=len, y=half)) + geom_point(size=3,alpha=0.5) + geom_errorbar(aes(ymin=half-half_se, ymax=half+half_se),size=1.5,alpha=0.25) + 
  geom_errorbarh(aes(xmin=len-len_se, xmax=len+len_se),size=1.5,alpha=0.25) +
  geom_text_repel(data=subset(pointrange.data, intron_num <= 6),aes(label=intron_num),color="darkred",segment.color=NA,force=2,nudge_y=-0.5,nudge_x=25) + 
  ylim(1,6.5) + labs(x="mean intron length (nt)",y="mean half-life (min)",color="intron number") + theme_bw() + guides(color=guide_colorbar(barwidth=5,barheight=0.5)) +
  theme(legend.position=c(0.5,0.1),legend.direction="horizontal",legend.text=element_text(size=7),legend.title=element_text(size=7))
#### oB - enhancers in first introns ####
fs8_b <- ggplot(enhancerfirst, aes(x=factor(len_bin),y=enhancer/total*100,fill=factor(type))) + geom_bar(stat="identity",position="dodge",color="white") + 
  scale_fill_manual(values=c(wes_palette("Royal1")[1], wes_palette("FantasticFox")[3])) + 
  scale_x_discrete(labels=c("43-59 nt","60-63 nt","64-68 nt","69-99 nt","99-2,493 nt","RIME ~ 1","57-353 nt","354-674 nt","675-1,252 nt","1,253-2,943 nt","> 2,944 nt")) + 
  labs(x="intron length (quintiles)",y="% introns \nwith enhancer",fill="intron position") + theme_bw() +
  theme(legend.position="bottom",axis.text.x=element_text(angle=45,hjust=1))
####r C - enhancers in all introns ####
combo.juncratio.data.parsed$len_bin_intron_def <- factor(combo.juncratio.data.parsed$len_bin_intron_def, levels=c(paste0(seq(20,100,20),"%_ID"),"confused",paste0(seq(20,100,20),"%_ED")))
combo.juncratio.data.parsed$arnold_binary <- factor(combo.juncratio.data.parsed$arnold_binary, levels=c("no enhancer","enhancer"))

fs8_c <- ggplot(combo.juncratio.data.parsed, aes(x=factor(IEratio_bin),y=fitvalue,fill=factor(arnold_binary))) + geom_boxplot(notch=T,outlier.color=alpha(wes_palette("Royal1")[1],0.05)) + 
  scale_fill_manual(values=c(wes_palette("Royal1")[1],wes_palette("FantasticFox")[3])) + scale_y_log10(limits=c(0.5,16),breaks=c(0.5,1,2,4,8,16)) +
  scale_x_discrete(labels=c("intron-defined","RIME ~ 1","exon-defined")) + 
  annotate("segment",x=c(0.8, 2.8), xend=c(1.2, 3.2), y=c(14, 14), yend=c(14, 14), size=0.25) +
  annotate("text",x=c(1, 3),y=c(15, 15), label=c("**","***")) + theme_bw() +
  labs(x="intron length (quintiles)",y="half-life (min)",fill="") + theme(legend.position="bottom",axis.text.x=element_text(angle=45,hjust=1))

#ggplot(combo.juncratio.data.parsed, aes(x=factor(len_bin_intron_def),y=fitvalue,fill=factor(arnold_binary))) + geom_boxplot(notch=T,outlier.color=alpha(wes_palette("Royal1")[1],0.05)) + 
#  scale_fill_manual(values=c(wes_palette("Royal1")[1],wes_palette("FantasticFox")[3])) + scale_y_log10(limits=c(0.5,16),breaks=c(0.5,1,2,4,8,16)) +
#  scale_x_discrete(labels=c("43-59 nt","60-63 nt","64-68 nt","69-99 nt","99-2,493 nt","RIME ~ 1","57-353 nt","354-674 nt","675-1,252 nt","1,253-2,943 nt","> 2,944 nt")) + 
#  annotate("segment",x=c(0.8, 1.8, 2.8, 3.8, 4.8), xend=c(1.2, 2.2, 3.2, 4.2, 5.2), y=c(20, 20, 25, 45, 60), yend=c(20, 20, 25, 45, 60), size=0.25) +
#  annotate("text",x=c(1, 2, 3, 4, 5),y=c(21, 21, 27, 47, 62), label=c("**","*","**","***","***")) +
#  labs(x="intron length (quintiles)",y="half-life (min)",fill="") + theme(legend.position="bottom",axis.text.x=element_text(angle=45,hjust=1))
###t# D - 1st intron correlation ####

fs8_d <- ggplot(combo.first.thirds.data, aes(x=factor(bin),y=mean)) + geom_bar(aes(alpha=factor(bin)),stat="identity",fill=wes_palette("FantasticFox")[3],width=0.85) + 
  geom_errorbar(aes(ymin=mean-sderr,ymax=mean+sderr),color=wes_palette("Royal1")[1],width=0.5) + scale_y_continuous(limits=c(0,4),breaks=c(0,2,4,6)) +
  annotate("segment",x=c(1,2,3),xend=c(2,3,4),y=c(3.3,3.1,2.9),yend=c(3.3,3.1,2.9),size=0.25) +
  annotate("text",x=c(1.5,2.5,3.5),y=c(3.5,3.2,3),label=c("**","***","***")) +theme_bw() +
  scale_alpha_manual(values=c(0.625,0.75,0.875,1),guide=F) + labs(x="first intron length (nt)",y="median half-life (min)") + theme(axis.text.x=element_text(angle=45,hjust=1))
##n## E - condition on 60-70nt introns vs. 1st intron length ####
fs8_e <- ggplot(combo.metadata.polyA.half.5intron.data, aes(x=first.len, y=remaining.halflife)) + geom_smooth(method="lm",se=F,color=wes_palette("FantasticFox")[5],size=3,alpha=0.5) + geom_point(size=3) + 
  scale_x_log10() + scale_y_log10() + theme_bw() +
  #geom_text_repel(aes(label=genename)) + 
  labs(x="first intron length",y="median half life (min)")
#i### F - frequency of remaining introns with 1st intron length ####
fs8_f <- ggplot(subset(combo.metadata.polyA.half.median, Number.of.annotated.introns <= 10 & Number.of.detected.introns > 0), aes(x=Length.first.intron,color=factor(Number.of.annotated.introns))) + 
  stat_ecdf() + scale_x_log10(label=comma) + scale_color_brewer(palette = "RdBu",guide=guide_legend(ncol=4,byrow=T)) + 
  labs(x="first intron length (nt)",y="cumulative frequency",color="number of annotated introns") + background_grid(major="y",minor="y") +
  theme(legend.position=c(0.7,0.2),axis.text.x=element_text(size=7),legend.text=element_text(size=7),legend.title=element_text(size=6))
#### COMBINE ####
pdf("~/Dropbox (MIT)/Projects/Adelman/timecourse/Figures/revisedfigures/suppfig8.pdf",width=8,height=10.5,useDingbats = F)
ggdraw() + 
  draw_plot(fs8_a, 0, 0.66, 0.5, 0.33) + draw_plot(fs8_b, 0.5, 0.66, 0.5, 0.33) + 
  draw_plot(fs8_c, 0, 0.33, 0.5, 0.33) + draw_plot(fs8_d, 0.5, 0.33, 0.5, 0.33) + 
  draw_plot(fs8_e, 0, 0, 0.5, 0.33) + draw_plot(fs8_f, 0.5, 0, 0.5, 0.33) +
  draw_plot_label(c("A","B","C","D","E","F"), c(0, 0.5, 0, 0.5, 0, 0.5), c(1,1,0.66,0.66, 0.33,0.33))
dev.off()

 

#### SUPP TABLE 1 - intron data ####

supptable1 <- data.frame(intron = combo.juncratio.data.parsed$intron,
                         gene = combo.juncratio.data.parsed$gene,
                         TPM = combo.juncratio.data.parsed$TPM_total,
                         PSI = combo.juncratio.data.parsed$PSI_total,
                         intron_position = combo.juncratio.data.parsed$intronnum + 1,
                         intron_length = combo.juncratio.data.parsed$intronlen,
                         intron_type = combo.juncratio.data.parsed$type,
                         ss5_maxEnt = combo.juncratio.data.parsed$ss5,
                         ss3_maxEnt = combo.juncratio.data.parsed$ss3,
                         contains_enhancer = combo.juncratio.data.parsed$arnold_enhancers,
                         upexon_length = combo.juncratio.data.parsed$upexon_len,
                         downexon_length = combo.juncratio.data.parsed$downexon_len,
                         three_length = combo.juncratio.data.parsed$threelength,
                         ie_count_5m = combo.juncratio.data.parsed$ie_count_5,
                         ie_count_10m = combo.juncratio.data.parsed$ie_count_10,
                         ie_count_20m = combo.juncratio.data.parsed$ie_count_20,
                         ee_count_5m = combo.juncratio.data.parsed$ee_count_5,
                         ee_count_10m = combo.juncratio.data.parsed$ee_count_10,
                         ee_count_20m = combo.juncratio.data.parsed$ee_count_20,
                         halflife = combo.juncratio.data.parsed$fitvalue,
                         halflife_error = combo.juncratio.data.parsed$se,
                         accuracy = combo.juncratio.data.parsed$accuracy)

write.table(supptable1, "~/Dropbox (MIT)/Collaboration/DrosophilaSplicing_manuscript/SplicingKinetics_manuscript/tables_supp/SupplementaryTable1.txt",sep="\t",quote=F,row.names=F,col.names=T)

#### SUPP TABLE 2 - ID gene ontology ####

genes.introndef.enrich <- summary(enrichGO(gene=genes.introndef.tr$ENTREZID, organism="fly", ont="BP", universe=testallgenes_tr$ENTREZID, pAdjustMethod="BH", qvalueCutoff=0.1, readable=TRUE))
genes.introndef.enrich <- genes.introndef.enrich[c(2:nrow(genes.introndef.enrich)),c(1:7,9)]

write.table(genes.introndef.enrich, "~/Dropbox (MIT)/Collaboration/DrosophilaSplicing_manuscript/SplicingKinetics_manuscript/tables_supp/SupplementaryTable2.txt",sep="\t",quote=F,row.names=F,col.names=T)

#### SUPP TABLE 3 - ED gene ontology ####

genes.exondef.enrich <- summary(enrichGO(gene=genes.exondef.tr$ENTREZID, organism="fly", ont="BP", universe=testallgenes_tr$ENTREZID, pAdjustMethod="BH", qvalueCutoff=0.1, readable=TRUE))
genes.exondef.enrich <- genes.exondef.enrich[c(2:nrow(genes.exondef.enrich)), c(1:7,9)]

write.table(genes.exondef.enrich, "~/Dropbox (MIT)/Collaboration/DrosophilaSplicing_manuscript/SplicingKinetics_manuscript/tables_supp/SupplementaryTable3.txt",sep="\t",quote=F,row.names=F,col.names=T)





