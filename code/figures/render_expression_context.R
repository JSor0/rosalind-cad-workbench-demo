#!/usr/bin/env Rscript
# Figure 3: source-reviewed individual measurements with direct mean labels.
# Design specification recorded before plotting in design/FIGURE_3_DESIGN.md.
# Uses only preserved measurements and means; no glyphs or anatomy are drawn.
# Usage: Rscript scripts/render_expression_context.R [package_directory] [output_directory]
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(png))
a <- commandArgs(trailingOnly=TRUE)
filearg <- grep("^--file=", commandArgs(), value=TRUE)
default_pkg <- normalizePath(file.path(dirname(sub("^--file=", "", filearg[1])), ".."))
pkg <- if(length(a)>=1) normalizePath(a[1]) else default_pkg
out <- if(length(a)>=2) normalizePath(a[2]) else pkg
# Display-only whitespace crops; source PNGs remain byte-for-byte unchanged.
# Bounds are zero-based half-open [x0,y0,x1,y1]. No anatomy is edited.
active_png <- readPNG(file.path(pkg,"assets/preserved_v7/approved_icons/tardigrade_active.png"))
tun_png <- readPNG(file.path(pkg,"assets/preserved_v7/approved_icons/tardigrade_tun.png"))
active_icon <- active_png[371:940,111:1155,,drop=FALSE]
tun_icon <- tun_png[316:950,251:1050,,drop=FALSE]
active_icon[,,4] <- active_icon[,,4]*.60 # approved compositional opacity only
source_values <- file.path(pkg,"source_data/figure3_expression_values.tsv")
source_summary <- file.path(pkg,"source_data/figure3_expression_summary.tsv")
d <- read.delim(source_values, stringsAsFactors=FALSE, check.names=FALSE)
s <- read.delim(source_summary, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(nrow(d)==6L,nrow(s)==1L,s$gene_id=="bHd16413",
  identical(d$sample_id,paste0("GSM247250",1:6)),
  identical(d$condition,c(rep("active",3),rep("tun",3))),
  all(d$individuals==10000L),all(d$gene_id=="bHd16413"),
  isTRUE(all.equal(d$tpm,c(1.565613,.997161,1.597606,7.91215,11.60757,10.53063))),
  abs(s$active_mean_tpm-1.38679333333)<1e-10,
  abs(s$tun_mean_tpm-10.0167833333)<1e-10,
  abs(s$tun_over_active_mean_ratio-7.22298203529)<1e-10)
# Read recorded means; do not recompute aggregation or tests.
means <- c(s$active_mean_tpm,s$tun_mean_tpm)
ink <- "#1D3238"; active <- "#D6A472"; tun <- "#B86725"; gray <- "#C8CED0"
muted <- "#687579"; gridgray <- "#E9ECED"; orange <- "#C8792C"
state_cols <- c(active,tun)
txt <- function(label,x,y,size=15,col=ink,font=1,just="left",rot=0) {
 grid.text(label,x=x,y=y,just=just,rot=rot,
  gp=gpar(fontsize=size,col=col,fontface=font,fontfamily="Helvetica"))
}
line <- function(x,y,col=gray,lwd=1,...) grid.lines(x=x,y=y,gp=gpar(col=col,lwd=lwd),...)
ypos <- function(y) .245+y/13*.520
paint <- function() {
 grid.newpage()
 txt("The X5S1-associated locus is expressed in public tardigrade RNA-seq",.046,.952,24,font=2)
 txt("Yoshida et al. 2017. Gene Expression Omnibus, GEO GSE94295",.046,.906,16,col=muted)
 # Compact physiological preparation, registered above the two groups.
 grid.raster(active_icon,x=.240,y=.804,width=.065,height=.065*1.6*570/1045,interpolate=TRUE)
 grid.raster(tun_icon,x=.510,y=.804,width=.045,height=.045*1.6*635/800,interpolate=TRUE)
 txt("Tun preparation: 48 h at 85% relative humidity",.375,.862,15,col=muted,just="centre")
 # Solid restrained process arrow, separated from both icons (Wong, PoV p15 Fig. 2c-e).
 grid.lines(x=c(.290,.462),y=c(.800,.800),
   gp=gpar(col=muted,fill=muted,lwd=1.7),
   arrow=arrow(length=unit(2.8,"mm"),type="closed"))
 # Shared linear TPM coordinates; sparse guides are secondary to observations.
 for(v in c(0,3,6,9,12)) {
  line(c(.124,.691),rep(ypos(v),2),if(v==0) gray else gridgray,if(v==0) 1 else .6)
  txt(as.character(v),.106,ypos(v),17,col=muted,just="right")
 }
 txt("Locus abundance (TPM)",.061,.516,18,just="centre",rot=90)
 centers <- c(.240,.510)
 for(i in 1:2) line(centers[i]+c(-.040,.040),rep(ypos(means[i]),2),state_cols[i],2.5)
 xs <- rep(centers,each=3)+rep(c(-.025,0,.025),2)
 grid.circle(x=unit(xs,"npc"),y=unit(ypos(d$tpm),"npc"),r=unit(2.4,"mm"),
             gp=gpar(fill=rep(state_cols,each=3),col="white",lwd=.8))
 # Direct labels bind each mean to its observed group.
 txt(paste0(sprintf("%.3f",means[1])," TPM"),.302,ypos(means[1]),18,col=ink,font=2)
 txt(paste0(sprintf("%.3f",means[2])," TPM"),.572,ypos(means[2]),18,col=ink,font=2)
 txt("Hydrated active",centers[1],.207,18,col=ink,font=2,just="centre")
 txt("Anhydrobiotic tun",centers[2],.207,18,col=ink,font=2,just="centre")
 txt("3 libraries",centers[1],.176,15,col=muted,just="centre")
 txt("3 libraries",centers[2],.176,15,col=muted,just="centre")
 # Compact supporting conclusion, no significance or causal claim.
 txt(paste0(sprintf("%.2f",s$tun_over_active_mean_ratio),"\u00d7"),.753,.780,24,col=orange,font=2)
 txt("higher mean abundance",.753,.729,17,font=2)
 txt("in tun",.753,.697,17,font=2)
 txt("Descriptive only;",.753,.647,15,col=muted)
 txt("no new significance test",.753,.617,15,col=muted)
 # Open aligned identity rows; product assignments remain unresolved.
 txt("X5S1 protein",.753,.505,16,col=orange,font=2)
 txt("UniProt A0A1W0X5S1",.753,.473,15)
 txt("Associated genomic locus",.753,.424,15,font=2)
 txt("NCBI locus tag BV898_03327",.753,.392,14)
 txt("Historical RNA-seq gene ID",.753,.343,15,font=2)
 txt("bHd16413",.753,.311,15)
 species_label <- textGrob("H. exemplaris",x=.753,y=.260,just="left",
     gp=gpar(fontsize=15,col=muted,fontface=3,fontfamily="Helvetica"))
 grid.draw(species_label)
 txt(" Z151",.753+convertWidth(grobWidth(species_label),"npc",valueOnly=TRUE),.260,15,col=muted)
 txt(".1/.2 product assignments unresolved",.753,.225,13,col=muted)
 txt("Each dot: one library pooled from 10,000 adults. Orange lines: means.",.124,.124,16,col=muted)
 txt("Deposited Kallisto abundance estimates",.124,.091,15,col=muted)
 txt("What does this protein actually do?",.046,.037,20,font=2)
}
f_pdf <- file.path(out,"FIGURE_3_EXPRESSION_CONTEXT.pdf")
f_png <- file.path(out,"FIGURE_3_EXPRESSION_CONTEXT.png")
stopifnot(!file.exists(f_pdf),!file.exists(f_png))
quartz(type="pdf",file=f_pdf,width=16,height=10,family="Helvetica",bg="white",
 title="The X5S1-associated locus is expressed in public tardigrade RNA-seq")
paint();dev.off()
png(f_png,width=3200,height=2000,res=200,type="quartz",bg="white",family="Helvetica")
paint();dev.off()
cat(R.version.string,"\n")
cat("Wrote ",f_pdf,"\n",f_png,"\n",sep="")
