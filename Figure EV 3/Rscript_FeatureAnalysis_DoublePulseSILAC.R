
# clear the workspace 
rm(list = ls())
#setting working directory
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

## LOAD LIBRARIES ##
library(plotrix)
library(grid)
library(gridExtra)
library(scales)
library(ggplot2)


## UNIVERSAL VARIABLES ##
mylevhmn=c(
"bioplex","piette_et_al","a2_bioID2","tcl","human")
## READ DATA ##
hmn_all = read.table("dat/human_dat.txt", header=T)
hmn_aas = read.table("dat/human_pme.aa", header=T)
hmn_pme = read.table("dat/human_pme.txt")

bioplex = read.table("src/bioplex.txt")
tcl = read.table("src/tcl.txt")
piette_et_al = read.table("src/piette_et_al.txt")
a2_bioID2 = read.table("src/a2_bioID2.txt")

hmn_all$abd_log <- log(hmn_all$abd,10)
hmn_all$rbp_cnt <- lapply(hmn_all$rbp, round)
hmn_all$pap_cnt <- ifelse(hmn_all$pap>=0.05,1,0)
hmn_all$plc_cnt <- ifelse(hmn_all$plc>0.0,1,0)
hmn_all$abd_log_cen <- hmn_all$abd_log-mean(na.omit(hmn_all$abd_log))
hmn_all$abd_log_nrm <- hmn_all$abd_log_cen/sd(na.omit(hmn_all$abd_log_cen))
hmn_all$csl_cen <- hmn_all$csl-mean(na.omit(hmn_all$csl))
hmn_all$csl_nrm <- hmn_all$csl_cen/sd(na.omit(hmn_all$csl_cen))
hmn_all$tsl     <- hmn_all$stc/hmn_all$len
hmn_all$sst     <- hmn_all$abd_log - hmn_all$csl_nrm

my.hmn      <- hmn_all[which(hmn_all$nam %in% hmn_pme$V1),]
my.hmn$tag <- "human"
my.aam      <- hmn_aas[which(hmn_aas$nam %in% hmn_pme $V1),]
my.aam$tag <- "human"
ref.hmn     <- my.hmn
ref.aam     <- my.aam

#######################
#######################

#
 my.tmp      <- ref.hmn[which(ref.hmn$nam %in% bioplex$V1),]
 my.tmp$tag <- "bioplex"
 my.hmn      <- rbind(my.hmn,my.tmp)
 my.tmp      <- ref.aam[which(ref.aam$nam %in% bioplex$V1),]
 my.tmp$tag <- "bioplex"
 my.aam      <- rbind(my.aam,my.tmp)
#
 my.tmp      <- ref.hmn[which(ref.hmn$nam %in% tcl$V1),]
 my.tmp$tag <- "tcl"
 my.hmn      <- rbind(my.hmn,my.tmp)
 my.tmp      <- ref.aam[which(ref.aam$nam %in% tcl$V1),]
 my.tmp$tag <- "tcl"
 my.aam      <- rbind(my.aam,my.tmp)
#
 my.tmp      <- ref.hmn[which(ref.hmn$nam %in% piette_et_al$V1),]
 my.tmp$tag <- "piette_et_al"
 my.hmn      <- rbind(my.hmn,my.tmp)
 my.tmp      <- ref.aam[which(ref.aam$nam %in% piette_et_al$V1),]
 my.tmp$tag <- "piette_et_al"
 my.aam      <- rbind(my.aam,my.tmp)
#
 my.tmp      <- ref.hmn[which(ref.hmn$nam %in% a2_bioID2$V1),]
 my.tmp$tag <- "a2_bioID2"
 my.hmn      <- rbind(my.hmn,my.tmp)
 my.tmp      <- ref.aam[which(ref.aam$nam %in% a2_bioID2$V1),]
 my.tmp$tag <- "a2_bioID2"
 my.aam      <- rbind(my.aam,my.tmp)
#
###########
# PLOT   
###########
pdf("plot_data.pdf")

my.hmn$tag <-factor(my.hmn$tag, levels=mylevhmn)
my.aam$tag <-factor(my.aam$tag, levels=mylevhmn)
## dis ## 
ggplot(my.hmn, aes(x=tag,y=dis)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="Disorder") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,100) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## len ## 
ggplot(my.hmn, aes(x=tag,y=len)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="Seq_Length") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,2000) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## run ## 
ggplot(my.hmn, aes(x=tag,y=run)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="n_IDRs") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,10) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## max ## 
ggplot(my.hmn, aes(x=tag,y=max)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="Max_IDR") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,600) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## chg ## 
ggplot(my.hmn, aes(x=tag,y=chg)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="%Charged_AA") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,60) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## net ## 
ggplot(my.hmn, aes(x=tag,y=net)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="Net_Charge") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(-50,50) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## abd_log ## 
ggplot(my.hmn, aes(x=tag,y=abd_log)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="Log_Abundance") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(-4,5) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## gvy ## 
ggplot(my.hmn, aes(x=tag,y=gvy)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="GRAVY") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(-3,2) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## pip ## 
ggplot(my.hmn, aes(x=tag,y=pip)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="PScore") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,8) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## iac ## 
ggplot(my.hmn, aes(x=tag,y=iac)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="Binary_Interactions") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,35) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## mfc ## 
ggplot(my.hmn, aes(x=tag,y=mfc)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="n_MoRFs") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## csl ## 
ggplot(my.hmn, aes(x=tag,y=csl)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="Camsol") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(-10,10) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## tsl ## 
ggplot(my.hmn, aes(x=tag,y=tsl)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="%hi_TANGO") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0.0,0.75) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## ssa ## 
ggplot(my.hmn, aes(x=tag,y=ssa)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="%Alpha") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,90) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## ssb ## 
ggplot(my.hmn, aes(x=tag,y=ssb)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="%Beta") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,90) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## ssc ## 
ggplot(my.hmn, aes(x=tag,y=ssc)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="%Coil") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,100) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## pol ## 
ggplot(my.hmn, aes(x=tag,y=pol)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="Polarizabilty") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(60,80) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## tmx ## 
ggplot(my.hmn, aes(x=tag,y=tmx)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="Tm") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(30,70) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## toh ## 
ggplot(my.hmn, aes(x=tag,y=toh)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="Turnover") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,100) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## sst ## 
ggplot(my.hmn, aes(x=tag,y=sst)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="SuperSaturation") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(-3,5) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## phs ## 
ggplot(my.hmn, aes(x=tag,y=phs)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="Annotated_PSites") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,30) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## mag ## 
ggplot(my.hmn, aes(x=tag,y=mag)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="MaGS") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(-3,5) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## rbp ## 
ggplot(my.hmn, aes(x=tag,y=rbp)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x="Set", y="RBPPred") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,1) + scale_fill_manual(values=c("firebrick","firebrick","firebrick","firebrick","gainsboro")) + theme(legend.position="none")
#
## pap_cnt ##
barplot(c(length(my.hmn[which(my.hmn$tag=="bioplex" & my.hmn$pap_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="bioplex"),]$nam),length(my.hmn[which(my.hmn$tag=="tcl" & my.hmn$pap_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="tcl"),]$nam),length(my.hmn[which(my.hmn$tag=="piette_et_al" & my.hmn$pap_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="piette_et_al"),]$nam),length(my.hmn[which(my.hmn$tag=="a2_bioID2" & my.hmn$pap_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="a2_bioID2"),]$nam),length(my.hmn[which(my.hmn$tag=="human" & my.hmn$pap_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="human"),]$nam)),ylab="pap_cnt Fraction",names=c("bioplex","tcl","piette_et_al","a2_bioID2","human"),col=c("firebrick","firebrick","firebrick","firebrick","gainsboro"))
## plc_cnt ##
barplot(c(length(my.hmn[which(my.hmn$tag=="bioplex" & my.hmn$plc_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="bioplex"),]$nam),length(my.hmn[which(my.hmn$tag=="tcl" & my.hmn$plc_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="tcl"),]$nam),length(my.hmn[which(my.hmn$tag=="piette_et_al" & my.hmn$plc_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="piette_et_al"),]$nam),length(my.hmn[which(my.hmn$tag=="a2_bioID2" & my.hmn$plc_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="a2_bioID2"),]$nam),length(my.hmn[which(my.hmn$tag=="human" & my.hmn$plc_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="human"),]$nam)),ylab="plc_cnt Fraction",names=c("bioplex","tcl","piette_et_al","a2_bioID2","human"),col=c("firebrick","firebrick","firebrick","firebrick","gainsboro"))
## rbp_cnt ##
barplot(c(length(my.hmn[which(my.hmn$tag=="bioplex" & my.hmn$rbp_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="bioplex"),]$nam),length(my.hmn[which(my.hmn$tag=="tcl" & my.hmn$rbp_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="tcl"),]$nam),length(my.hmn[which(my.hmn$tag=="piette_et_al" & my.hmn$rbp_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="piette_et_al"),]$nam),length(my.hmn[which(my.hmn$tag=="a2_bioID2" & my.hmn$rbp_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="a2_bioID2"),]$nam),length(my.hmn[which(my.hmn$tag=="human" & my.hmn$rbp_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="human"),]$nam)),ylab="rbp_cnt Fraction",names=c("bioplex","tcl","piette_et_al","a2_bioID2","human"),col=c("firebrick","firebrick","firebrick","firebrick","gainsboro"))
dev.off()
#

###########
# P-Values
###########
#
p_dis <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$dis, my.hmn[which(my.hmn$tag=='human'),]$dis)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$dis, my.hmn[which(my.hmn$tag=='human'),]$dis)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$dis, my.hmn[which(my.hmn$tag=='tcl'),]$dis)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$dis, my.hmn[which(my.hmn$tag=='tcl'),]$dis)$p.val)
p_len <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$len, my.hmn[which(my.hmn$tag=='human'),]$len)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$len, my.hmn[which(my.hmn$tag=='human'),]$len)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$len, my.hmn[which(my.hmn$tag=='tcl'),]$len)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$len, my.hmn[which(my.hmn$tag=='tcl'),]$len)$p.val)
p_run <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$run, my.hmn[which(my.hmn$tag=='human'),]$run)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$run, my.hmn[which(my.hmn$tag=='human'),]$run)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$run, my.hmn[which(my.hmn$tag=='tcl'),]$run)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$run, my.hmn[which(my.hmn$tag=='tcl'),]$run)$p.val)
p_max <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$max, my.hmn[which(my.hmn$tag=='human'),]$max)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$max, my.hmn[which(my.hmn$tag=='human'),]$max)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$max, my.hmn[which(my.hmn$tag=='tcl'),]$max)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$max, my.hmn[which(my.hmn$tag=='tcl'),]$max)$p.val)
p_chg <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$chg, my.hmn[which(my.hmn$tag=='human'),]$chg)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$chg, my.hmn[which(my.hmn$tag=='human'),]$chg)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$chg, my.hmn[which(my.hmn$tag=='tcl'),]$chg)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$chg, my.hmn[which(my.hmn$tag=='tcl'),]$chg)$p.val)
p_net <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$net, my.hmn[which(my.hmn$tag=='human'),]$net)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$net, my.hmn[which(my.hmn$tag=='human'),]$net)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$net, my.hmn[which(my.hmn$tag=='tcl'),]$net)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$net, my.hmn[which(my.hmn$tag=='tcl'),]$net)$p.val)
p_abd_log <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$abd_log, my.hmn[which(my.hmn$tag=='human'),]$abd_log)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$abd_log, my.hmn[which(my.hmn$tag=='human'),]$abd_log)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$abd_log, my.hmn[which(my.hmn$tag=='tcl'),]$abd_log)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$abd_log, my.hmn[which(my.hmn$tag=='tcl'),]$abd_log)$p.val)
p_gvy <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$gvy, my.hmn[which(my.hmn$tag=='human'),]$gvy)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$gvy, my.hmn[which(my.hmn$tag=='human'),]$gvy)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$gvy, my.hmn[which(my.hmn$tag=='tcl'),]$gvy)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$gvy, my.hmn[which(my.hmn$tag=='tcl'),]$gvy)$p.val)
p_pip <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$pip, my.hmn[which(my.hmn$tag=='human'),]$pip)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$pip, my.hmn[which(my.hmn$tag=='human'),]$pip)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$pip, my.hmn[which(my.hmn$tag=='tcl'),]$pip)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$pip, my.hmn[which(my.hmn$tag=='tcl'),]$pip)$p.val)
p_iac <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$iac, my.hmn[which(my.hmn$tag=='human'),]$iac)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$iac, my.hmn[which(my.hmn$tag=='human'),]$iac)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$iac, my.hmn[which(my.hmn$tag=='tcl'),]$iac)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$iac, my.hmn[which(my.hmn$tag=='tcl'),]$iac)$p.val)
p_mfc <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$mfc, my.hmn[which(my.hmn$tag=='human'),]$mfc)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$mfc, my.hmn[which(my.hmn$tag=='human'),]$mfc)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$mfc, my.hmn[which(my.hmn$tag=='tcl'),]$mfc)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$mfc, my.hmn[which(my.hmn$tag=='tcl'),]$mfc)$p.val)
p_csl <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$csl, my.hmn[which(my.hmn$tag=='human'),]$csl)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$csl, my.hmn[which(my.hmn$tag=='human'),]$csl)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$csl, my.hmn[which(my.hmn$tag=='tcl'),]$csl)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$csl, my.hmn[which(my.hmn$tag=='tcl'),]$csl)$p.val)
p_tsl <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$tsl, my.hmn[which(my.hmn$tag=='human'),]$tsl)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$tsl, my.hmn[which(my.hmn$tag=='human'),]$tsl)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$tsl, my.hmn[which(my.hmn$tag=='tcl'),]$tsl)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$tsl, my.hmn[which(my.hmn$tag=='tcl'),]$tsl)$p.val)
p_ssa <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$ssa, my.hmn[which(my.hmn$tag=='human'),]$ssa)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$ssa, my.hmn[which(my.hmn$tag=='human'),]$ssa)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$ssa, my.hmn[which(my.hmn$tag=='tcl'),]$ssa)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$ssa, my.hmn[which(my.hmn$tag=='tcl'),]$ssa)$p.val)
p_ssb <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$ssb, my.hmn[which(my.hmn$tag=='human'),]$ssb)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$ssb, my.hmn[which(my.hmn$tag=='human'),]$ssb)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$ssb, my.hmn[which(my.hmn$tag=='tcl'),]$ssb)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$ssb, my.hmn[which(my.hmn$tag=='tcl'),]$ssb)$p.val)
p_ssc <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$ssc, my.hmn[which(my.hmn$tag=='human'),]$ssc)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$ssc, my.hmn[which(my.hmn$tag=='human'),]$ssc)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$ssc, my.hmn[which(my.hmn$tag=='tcl'),]$ssc)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$ssc, my.hmn[which(my.hmn$tag=='tcl'),]$ssc)$p.val)
p_pol <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$pol, my.hmn[which(my.hmn$tag=='human'),]$pol)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$pol, my.hmn[which(my.hmn$tag=='human'),]$pol)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$pol, my.hmn[which(my.hmn$tag=='tcl'),]$pol)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$pol, my.hmn[which(my.hmn$tag=='tcl'),]$pol)$p.val)
p_tmx <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$tmx, my.hmn[which(my.hmn$tag=='human'),]$tmx)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$tmx, my.hmn[which(my.hmn$tag=='human'),]$tmx)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$tmx, my.hmn[which(my.hmn$tag=='tcl'),]$tmx)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$tmx, my.hmn[which(my.hmn$tag=='tcl'),]$tmx)$p.val)
p_toh <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$toh, my.hmn[which(my.hmn$tag=='human'),]$toh)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$toh, my.hmn[which(my.hmn$tag=='human'),]$toh)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$toh, my.hmn[which(my.hmn$tag=='tcl'),]$toh)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$toh, my.hmn[which(my.hmn$tag=='tcl'),]$toh)$p.val)
p_sst <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$sst, my.hmn[which(my.hmn$tag=='human'),]$sst)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$sst, my.hmn[which(my.hmn$tag=='human'),]$sst)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$sst, my.hmn[which(my.hmn$tag=='tcl'),]$sst)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$sst, my.hmn[which(my.hmn$tag=='tcl'),]$sst)$p.val)
p_phs <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$phs, my.hmn[which(my.hmn$tag=='human'),]$phs)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$phs, my.hmn[which(my.hmn$tag=='human'),]$phs)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$phs, my.hmn[which(my.hmn$tag=='tcl'),]$phs)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$phs, my.hmn[which(my.hmn$tag=='tcl'),]$phs)$p.val)
p_mag <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$mag, my.hmn[which(my.hmn$tag=='human'),]$mag)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$mag, my.hmn[which(my.hmn$tag=='human'),]$mag)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$mag, my.hmn[which(my.hmn$tag=='tcl'),]$mag)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$mag, my.hmn[which(my.hmn$tag=='tcl'),]$mag)$p.val)
p_rbp <- c(wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$rbp, my.hmn[which(my.hmn$tag=='human'),]$rbp)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$rbp, my.hmn[which(my.hmn$tag=='human'),]$rbp)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='piette_et_al'),]$rbp, my.hmn[which(my.hmn$tag=='tcl'),]$rbp)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='a2_bioID2'),]$rbp, my.hmn[which(my.hmn$tag=='tcl'),]$rbp)$p.val)
p_hmn <- rbind(p_dis,p_len,p_run,p_max,p_chg,p_net,p_abd_log,p_gvy,p_pip,p_iac,p_mfc,p_csl,p_tsl,p_ssa,p_ssb,p_ssc,p_pol,p_tmx,p_toh,p_sst,p_phs,p_mag,p_rbp)
#
piette_et_al_vs_pme <- p.adjust(p_hmn[,1],method="hochberg")
a2_bioID2_vs_pme <- p.adjust(p_hmn[,2],method="hochberg")
piette_et_al_vs_tcl <- p.adjust(p_hmn[,3],method="hochberg")
a2_bioID2_vs_tcl <- p.adjust(p_hmn[,4],method="hochberg")

props<-c("dis","len","run","max","chg","net","abd_log","gvy","pip","iac","mfc","csl","tsl","ssa","ssb","ssc","pol","tmx","toh","sst","phs","mag","rbp")
bioplex_vs_pme<-p.adjust(sapply(props,function(x)wilcox.test(my.hmn[which(my.hmn$tag=="bioplex"),x],my.hmn[which(my.hmn$tag=="human"),x])$p.value),method="hochberg")
bioplex_vs_tcl<-p.adjust(sapply(props,function(x)wilcox.test(my.hmn[which(my.hmn$tag=="bioplex"),x],my.hmn[which(my.hmn$tag=="tcl"),x])$p.value),method="hochberg")
corr_pval<-cbind(piette_et_al_vs_pme,a2_bioID2_vs_pme,bioplex_vs_pme,piette_et_al_vs_tcl,a2_bioID2_vs_tcl,bioplex_vs_tcl)
 
ptab<-signif(corr_pval,digits=2)

pdf("plot_pval.pdf",width=16,height=10)
grid.table(ptab)
dev.off()

#

