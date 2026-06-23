
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
"tcl_duplicate","ns_lc","other_lc","tcl","human")
## READ DATA ##
hmn_all = read.table("dat/human_dat.txt", header=T)
hmn_aas = read.table("dat/human_pme.aa", header=T)
hmn_pme = read.table("dat/human_pme.txt")

tcl_duplicate = read.table("src/tcl_duplicate.txt")
tcl = read.table("src/tcl.txt")
ns_lc = read.table("src/ns_lc.txt")
other_lc = read.table("src/other_lc.txt")

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
 my.tmp      <- ref.hmn[which(ref.hmn$nam %in% tcl_duplicate$V1),]
 my.tmp$tag <- "tcl_duplicate"
 my.hmn      <- rbind(my.hmn,my.tmp)
 my.tmp      <- ref.aam[which(ref.aam$nam %in% tcl_duplicate$V1),]
 my.tmp$tag <- "tcl_duplicate"
 my.aam      <- rbind(my.aam,my.tmp)
#
 my.tmp      <- ref.hmn[which(ref.hmn$nam %in% tcl$V1),]
 my.tmp$tag <- "tcl"
 my.hmn      <- rbind(my.hmn,my.tmp)
 my.tmp      <- ref.aam[which(ref.aam$nam %in% tcl$V1),]
 my.tmp$tag <- "tcl"
 my.aam      <- rbind(my.aam,my.tmp)
#
 my.tmp      <- ref.hmn[which(ref.hmn$nam %in% ns_lc$V1),]
 my.tmp$tag <- "ns_lc"
 my.hmn      <- rbind(my.hmn,my.tmp)
 my.tmp      <- ref.aam[which(ref.aam$nam %in% ns_lc$V1),]
 my.tmp$tag <- "ns_lc"
 my.aam      <- rbind(my.aam,my.tmp)
#
 my.tmp      <- ref.hmn[which(ref.hmn$nam %in% other_lc$V1),]
 my.tmp$tag <- "other_lc"
 my.hmn      <- rbind(my.hmn,my.tmp)
 my.tmp      <- ref.aam[which(ref.aam$nam %in% other_lc$V1),]
 my.tmp$tag <- "other_lc"
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
barplot(c(length(my.hmn[which(my.hmn$tag=="tcl_duplicate" & my.hmn$pap_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="tcl_duplicate"),]$nam),length(my.hmn[which(my.hmn$tag=="tcl" & my.hmn$pap_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="tcl"),]$nam),length(my.hmn[which(my.hmn$tag=="ns_lc" & my.hmn$pap_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="ns_lc"),]$nam),length(my.hmn[which(my.hmn$tag=="other_lc" & my.hmn$pap_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="other_lc"),]$nam),length(my.hmn[which(my.hmn$tag=="human" & my.hmn$pap_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="human"),]$nam)),ylab="pap_cnt Fraction",names=c("tcl_duplicate","tcl","ns_lc","other_lc","human"),col=c("firebrick","firebrick","firebrick","firebrick","gainsboro"))
## plc_cnt ##
barplot(c(length(my.hmn[which(my.hmn$tag=="tcl_duplicate" & my.hmn$plc_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="tcl_duplicate"),]$nam),length(my.hmn[which(my.hmn$tag=="tcl" & my.hmn$plc_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="tcl"),]$nam),length(my.hmn[which(my.hmn$tag=="ns_lc" & my.hmn$plc_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="ns_lc"),]$nam),length(my.hmn[which(my.hmn$tag=="other_lc" & my.hmn$plc_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="other_lc"),]$nam),length(my.hmn[which(my.hmn$tag=="human" & my.hmn$plc_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="human"),]$nam)),ylab="plc_cnt Fraction",names=c("tcl_duplicate","tcl","ns_lc","other_lc","human"),col=c("firebrick","firebrick","firebrick","firebrick","gainsboro"))
## rbp_cnt ##
barplot(c(length(my.hmn[which(my.hmn$tag=="tcl_duplicate" & my.hmn$rbp_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="tcl_duplicate"),]$nam),length(my.hmn[which(my.hmn$tag=="tcl" & my.hmn$rbp_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="tcl"),]$nam),length(my.hmn[which(my.hmn$tag=="ns_lc" & my.hmn$rbp_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="ns_lc"),]$nam),length(my.hmn[which(my.hmn$tag=="other_lc" & my.hmn$rbp_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="other_lc"),]$nam),length(my.hmn[which(my.hmn$tag=="human" & my.hmn$rbp_cnt==1),]$nam) / length(my.hmn[which(my.hmn$tag=="human"),]$nam)),ylab="rbp_cnt Fraction",names=c("tcl_duplicate","tcl","ns_lc","other_lc","human"),col=c("firebrick","firebrick","firebrick","firebrick","gainsboro"))
dev.off()
#
###############
# AA Plot
###############
pdf("plot_aa.pdf")
## R ##
ggplot(my.aam, aes(x=tag,y=R)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="R") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## H ##
ggplot(my.aam, aes(x=tag,y=H)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="H") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## K ##
ggplot(my.aam, aes(x=tag,y=K)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="K") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## D ##
ggplot(my.aam, aes(x=tag,y=D)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="D") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## E ##
ggplot(my.aam, aes(x=tag,y=E)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="E") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## S ##
ggplot(my.aam, aes(x=tag,y=S)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="S") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## T ##
ggplot(my.aam, aes(x=tag,y=T)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="T") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## N ##
ggplot(my.aam, aes(x=tag,y=N)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="N") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## Q ##
ggplot(my.aam, aes(x=tag,y=Q)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="Q") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## C ##
ggplot(my.aam, aes(x=tag,y=C)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="C") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## U ##
ggplot(my.aam, aes(x=tag,y=U)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="U") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## G ##
ggplot(my.aam, aes(x=tag,y=G)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="G") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## P ##
ggplot(my.aam, aes(x=tag,y=P)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="P") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## A ##
ggplot(my.aam, aes(x=tag,y=A)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="A") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## V ##
ggplot(my.aam, aes(x=tag,y=V)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="V") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## I ##
ggplot(my.aam, aes(x=tag,y=I)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="I") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## L ##
ggplot(my.aam, aes(x=tag,y=L)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="L") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## M ##
ggplot(my.aam, aes(x=tag,y=M)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="M") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## F ##
ggplot(my.aam, aes(x=tag,y=F)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="F") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## Y ##
ggplot(my.aam, aes(x=tag,y=Y)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="Y") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
## W ##
ggplot(my.aam, aes(x=tag,y=W)) + geom_violin(trim=FALSE, fill=c("firebrick2"),adjust=1.5) + labs(x="Set", y="W") + geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#
dev.off()
#
###########
# P-Values
###########
#
p_dis <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$dis, my.hmn[which(my.hmn$tag=='human'),]$dis)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$dis, my.hmn[which(my.hmn$tag=='human'),]$dis)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$dis, my.hmn[which(my.hmn$tag=='tcl'),]$dis)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$dis, my.hmn[which(my.hmn$tag=='tcl'),]$dis)$p.val)
p_len <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$len, my.hmn[which(my.hmn$tag=='human'),]$len)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$len, my.hmn[which(my.hmn$tag=='human'),]$len)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$len, my.hmn[which(my.hmn$tag=='tcl'),]$len)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$len, my.hmn[which(my.hmn$tag=='tcl'),]$len)$p.val)
p_run <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$run, my.hmn[which(my.hmn$tag=='human'),]$run)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$run, my.hmn[which(my.hmn$tag=='human'),]$run)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$run, my.hmn[which(my.hmn$tag=='tcl'),]$run)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$run, my.hmn[which(my.hmn$tag=='tcl'),]$run)$p.val)
p_max <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$max, my.hmn[which(my.hmn$tag=='human'),]$max)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$max, my.hmn[which(my.hmn$tag=='human'),]$max)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$max, my.hmn[which(my.hmn$tag=='tcl'),]$max)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$max, my.hmn[which(my.hmn$tag=='tcl'),]$max)$p.val)
p_chg <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$chg, my.hmn[which(my.hmn$tag=='human'),]$chg)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$chg, my.hmn[which(my.hmn$tag=='human'),]$chg)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$chg, my.hmn[which(my.hmn$tag=='tcl'),]$chg)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$chg, my.hmn[which(my.hmn$tag=='tcl'),]$chg)$p.val)
p_net <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$net, my.hmn[which(my.hmn$tag=='human'),]$net)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$net, my.hmn[which(my.hmn$tag=='human'),]$net)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$net, my.hmn[which(my.hmn$tag=='tcl'),]$net)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$net, my.hmn[which(my.hmn$tag=='tcl'),]$net)$p.val)
p_abd_log <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$abd_log, my.hmn[which(my.hmn$tag=='human'),]$abd_log)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$abd_log, my.hmn[which(my.hmn$tag=='human'),]$abd_log)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$abd_log, my.hmn[which(my.hmn$tag=='tcl'),]$abd_log)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$abd_log, my.hmn[which(my.hmn$tag=='tcl'),]$abd_log)$p.val)
p_gvy <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$gvy, my.hmn[which(my.hmn$tag=='human'),]$gvy)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$gvy, my.hmn[which(my.hmn$tag=='human'),]$gvy)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$gvy, my.hmn[which(my.hmn$tag=='tcl'),]$gvy)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$gvy, my.hmn[which(my.hmn$tag=='tcl'),]$gvy)$p.val)
p_pip <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$pip, my.hmn[which(my.hmn$tag=='human'),]$pip)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$pip, my.hmn[which(my.hmn$tag=='human'),]$pip)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$pip, my.hmn[which(my.hmn$tag=='tcl'),]$pip)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$pip, my.hmn[which(my.hmn$tag=='tcl'),]$pip)$p.val)
p_iac <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$iac, my.hmn[which(my.hmn$tag=='human'),]$iac)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$iac, my.hmn[which(my.hmn$tag=='human'),]$iac)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$iac, my.hmn[which(my.hmn$tag=='tcl'),]$iac)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$iac, my.hmn[which(my.hmn$tag=='tcl'),]$iac)$p.val)
p_mfc <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$mfc, my.hmn[which(my.hmn$tag=='human'),]$mfc)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$mfc, my.hmn[which(my.hmn$tag=='human'),]$mfc)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$mfc, my.hmn[which(my.hmn$tag=='tcl'),]$mfc)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$mfc, my.hmn[which(my.hmn$tag=='tcl'),]$mfc)$p.val)
p_csl <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$csl, my.hmn[which(my.hmn$tag=='human'),]$csl)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$csl, my.hmn[which(my.hmn$tag=='human'),]$csl)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$csl, my.hmn[which(my.hmn$tag=='tcl'),]$csl)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$csl, my.hmn[which(my.hmn$tag=='tcl'),]$csl)$p.val)
p_tsl <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$tsl, my.hmn[which(my.hmn$tag=='human'),]$tsl)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$tsl, my.hmn[which(my.hmn$tag=='human'),]$tsl)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$tsl, my.hmn[which(my.hmn$tag=='tcl'),]$tsl)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$tsl, my.hmn[which(my.hmn$tag=='tcl'),]$tsl)$p.val)
p_ssa <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$ssa, my.hmn[which(my.hmn$tag=='human'),]$ssa)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$ssa, my.hmn[which(my.hmn$tag=='human'),]$ssa)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$ssa, my.hmn[which(my.hmn$tag=='tcl'),]$ssa)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$ssa, my.hmn[which(my.hmn$tag=='tcl'),]$ssa)$p.val)
p_ssb <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$ssb, my.hmn[which(my.hmn$tag=='human'),]$ssb)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$ssb, my.hmn[which(my.hmn$tag=='human'),]$ssb)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$ssb, my.hmn[which(my.hmn$tag=='tcl'),]$ssb)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$ssb, my.hmn[which(my.hmn$tag=='tcl'),]$ssb)$p.val)
p_ssc <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$ssc, my.hmn[which(my.hmn$tag=='human'),]$ssc)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$ssc, my.hmn[which(my.hmn$tag=='human'),]$ssc)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$ssc, my.hmn[which(my.hmn$tag=='tcl'),]$ssc)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$ssc, my.hmn[which(my.hmn$tag=='tcl'),]$ssc)$p.val)
p_pol <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$pol, my.hmn[which(my.hmn$tag=='human'),]$pol)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$pol, my.hmn[which(my.hmn$tag=='human'),]$pol)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$pol, my.hmn[which(my.hmn$tag=='tcl'),]$pol)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$pol, my.hmn[which(my.hmn$tag=='tcl'),]$pol)$p.val)
p_tmx <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$tmx, my.hmn[which(my.hmn$tag=='human'),]$tmx)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$tmx, my.hmn[which(my.hmn$tag=='human'),]$tmx)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$tmx, my.hmn[which(my.hmn$tag=='tcl'),]$tmx)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$tmx, my.hmn[which(my.hmn$tag=='tcl'),]$tmx)$p.val)
p_toh <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$toh, my.hmn[which(my.hmn$tag=='human'),]$toh)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$toh, my.hmn[which(my.hmn$tag=='human'),]$toh)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$toh, my.hmn[which(my.hmn$tag=='tcl'),]$toh)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$toh, my.hmn[which(my.hmn$tag=='tcl'),]$toh)$p.val)
p_sst <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$sst, my.hmn[which(my.hmn$tag=='human'),]$sst)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$sst, my.hmn[which(my.hmn$tag=='human'),]$sst)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$sst, my.hmn[which(my.hmn$tag=='tcl'),]$sst)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$sst, my.hmn[which(my.hmn$tag=='tcl'),]$sst)$p.val)
p_phs <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$phs, my.hmn[which(my.hmn$tag=='human'),]$phs)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$phs, my.hmn[which(my.hmn$tag=='human'),]$phs)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$phs, my.hmn[which(my.hmn$tag=='tcl'),]$phs)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$phs, my.hmn[which(my.hmn$tag=='tcl'),]$phs)$p.val)
p_mag <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$mag, my.hmn[which(my.hmn$tag=='human'),]$mag)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$mag, my.hmn[which(my.hmn$tag=='human'),]$mag)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$mag, my.hmn[which(my.hmn$tag=='tcl'),]$mag)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$mag, my.hmn[which(my.hmn$tag=='tcl'),]$mag)$p.val)
p_rbp <- c(wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$rbp, my.hmn[which(my.hmn$tag=='human'),]$rbp)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$rbp, my.hmn[which(my.hmn$tag=='human'),]$rbp)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='ns_lc'),]$rbp, my.hmn[which(my.hmn$tag=='tcl'),]$rbp)$p.val,wilcox.test(my.hmn[which(my.hmn$tag=='other_lc'),]$rbp, my.hmn[which(my.hmn$tag=='tcl'),]$rbp)$p.val)
p_hmn <- rbind(p_dis,p_len,p_run,p_max,p_chg,p_net,p_abd_log,p_gvy,p_pip,p_iac,p_mfc,p_csl,p_tsl,p_ssa,p_ssb,p_ssc,p_pol,p_tmx,p_toh,p_sst,p_phs,p_mag,p_rbp)
#
ns_lc_vs_pme <- p.adjust(p_hmn[,1],method="hochberg")
other_lc_vs_pme <- p.adjust(p_hmn[,2],method="hochberg")
ns_lc_vs_tcl <- p.adjust(p_hmn[,3],method="hochberg")
other_lc_vs_tcl <- p.adjust(p_hmn[,4],method="hochberg")
corr_pval<-cbind(ns_lc_vs_pme,other_lc_vs_pme,ns_lc_vs_tcl,other_lc_vs_tcl)
 
ptab<-signif(corr_pval,digits=2)

pdf("plot_pval.pdf")
grid.table(ptab)
dev.off()
#

