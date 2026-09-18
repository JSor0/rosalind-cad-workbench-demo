suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(png))
suppressPackageStartupMessages(library(jsonlite))
a<-commandArgs(trailingOnly=TRUE)
pkg<-normalizePath(if(length(a)>0)a[1] else '.')
out<-if(length(a)>1)normalizePath(a[2]) else pkg
mol<-file.path(pkg,'assets/hero')
s<-read.delim(file.path(pkg,'source_data/FIGURE1_SEARCH_DATA.tsv'),check.names=FALSE)
fam<-fromJSON(file.path(pkg,'source_data/FIGURE1_FAMILY_DATA.json'))
stopifnot(nrow(s)==2,s$returned_count[1]==307,s$provider_rank[1]==290,s$returned_count[2]==74,s$provider_rank[2]==1,fam$count==22,abs(s$evalue[1]-8.044e-7)<1e-15,abs(s$evalue[2]-3.069e-6)<1e-15,s$reported_identity[1]==.151,s$significant_count[2]==1,s$query_length[1]==197,s$significant_count[1]==298,s$reported_identity[2]==.148)
trim_white<-function(im){
 fg<-pmin(im[,,1],im[,,2],im[,,3])<.995
 if(dim(im)[3]==4)fg<-fg & im[,,4]>.01
 yy<-which(rowSums(fg)>0);xx<-which(colSums(fg)>0)
 b<-c(max(1,min(yy)-25),min(nrow(im),max(yy)+25),max(1,min(xx)-25),min(ncol(im),max(xx)+25))
 list(image=im[b[1]:b[2],b[3]:b[4],,drop=FALSE],bounds=b)
}
mouse_raw<-readPNG(file.path(mol,'mouse_c3.png'));x5_raw<-readPNG(file.path(mol,'x5s1_full.png'))
mt<-trim_white(mouse_raw);xt<-trim_white(x5_raw);mouse<-mt$image;x5<-xt$image
write.table(data.frame(asset=c('mouse_c3.png','x5s1_full.png'),row_first=c(mt$bounds[1],xt$bounds[1]),row_last=c(mt$bounds[2],xt$bounds[2]),column_first=c(mt$bounds[3],xt$bounds[3]),column_last=c(mt$bounds[4],xt$bounds[4])),file.path(out,'THUMBNAIL_CROP_BOUNDS.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
ink<-'#22363C';muted<-'#627277';teal<-'#247F82';orange<-'#C8792C';rule<-'#D6DEDF'
txt<-function(label,x,y,size=16,col=ink,font=1,just='left')grid.text(label,x=x,y=y,just=just,gp=gpar(fontsize=size,col=col,fontface=font,fontfamily='Helvetica'))
line<-function(x,y,col=rule,lwd=.7)grid.lines(x=x,y=y,gp=gpar(col=col,lwd=lwd))
arr<-function(x1,x2,y){grid.lines(x=c(x1,x2),y=c(y,y),arrow=arrow(length=unit(.110,'inches'),type='closed'),gp=gpar(col=ink,fill=ink,lwd=1.7))}
ras<-function(im,x,y,w=.25,h=.32){ar<-dim(im)[2]/dim(im)[1];wi<-min(w*16,h*10*ar);hi<-wi/ar;grid.raster(im,x=x,y=y,width=unit(wi,'inches'),height=unit(hi,'inches'),interpolate=TRUE)}
# Approved identity artwork: whitespace crop and proportional display only.
icon_dir<-file.path(pkg,'assets/preserved_v7/approved_icons')
mouse_icon<-readPNG(file.path(icon_dir,'mouse_teal.png'))[371:960,91:1165,,drop=FALSE]
x5_icon<-readPNG(file.path(icon_dir,'tardigrade_active.png'))[371:940,111:1155,,drop=FALSE]
draw<-function(){
 grid.newpage();grid.rect(gp=gpar(fill='white',col=NA))
 txt('Rosalind follows a divergent protein back to CAD',.045,.954,26,font=2)
 txt('Two Foldseek searches, followed by a 22-structure FoldMason alignment',.045,.908,17,col=muted)
 txt('Rosalind Workbench with Tamarind Bio',.045,.853,14.5,font=2)
 # Repeated noun and status rows keep the experiment/prediction identities stable.
 ras(mouse_icon,.075,.790,.056,.067)
 txt('Mouse CAD C3',.115,.790,21,col=teal,font=2)
 txt('Experimental query',.155,.756,14.5,col=muted,just='center')
 txt('PDB 1V0D, chain A residues 132-328',.155,.724,13.3,col=muted,just='center')
 ras(x5_icon,.472,.790,.053,.061)
 txt('X5S1',.508,.790,23,col=orange,font=2)
 txt('Full-length predicted monomer',.500,.756,14.5,col=muted,just='center')
 txt('UniProt A0A1W0X5S1',.500,.724,13.3,col=muted,just='center')
 ras(mouse_icon,.781,.790,.056,.067)
 txt('Mouse CAD',.821,.790,21,col=teal,font=2)
 txt('Experimental structure',.845,.756,14.5,col=muted,just='center')
 txt('PDB 1V0D, C3 shown',.845,.724,13.3,col=muted,just='center')
 # Structure thumbnails use the original archive, not AI or schematic substitutes.
 ras(mouse,.155,.535,.245,.345)
 ras(x5,.500,.535,.245,.345)
 ras(mouse,.845,.535,.245,.345)
 # Arrows denote completed database searches, not protein conversion.
 txt('Foldseek',.329,.670,15.8,font=2,just='center')
 txt('Predicted structures',.329,.639,13.3,col=muted,just='center')
 txt('AlphaFold-UniProt50',.329,.610,12.5,col=muted,just='center')
 arr(.281,.376,.548)
 txt('Foldseek',.674,.670,15.8,font=2,just='center')
 txt('Experimental structures',.674,.639,13.3,col=muted,just='center')
 txt('Provider PDB database',.674,.610,12.5,col=muted,just='center')
 arr(.626,.721,.548)
 # Directly owned result labels establish returned rank and significance.
 txt('197-residue C3 query',.155,.339,15.2,col=muted,just='center')
 txt(expression(plain('298 hits at ') * italic(E) <= 0.001),.155,.307,13.3,col=muted,just='center')
 txt('Two significant hits below 20% identity',.155,.282,12.5,col=muted,just='center')
 txt('Rank 290 / 307',.500,.349,19,col=orange,font=2,just='center')
 txt(expression(italic(E)==8.044%*%10^{-7}),.500,.313,16,just='center')
 txt('15.1% reported alignment identity',.500,.282,13.3,col=muted,just='center')
 txt('Rank 1 / 74',.845,.349,19,col=teal,font=2,just='center')
 txt(expression(italic(E)==3.069%*%10^{-6}),.845,.313,16,just='center')
 txt('14.8% reported alignment identity',.845,.282,13.3,col=muted,just='center')
 txt(expression(plain('Only returned hit at ') * italic(E) <= 0.001),.845,.257,12.0,col=muted,just='center')
 line(c(.045,.955),c(.241,.241),lwd=.6)
 # A new family-scale native operation, with finite input and explicit outputs.
 txt('Then compare the family',.045,.209,19,font=2)
 txt('22 existing predicted structures',.045,.161,17,font=2)
 txt(expression(plain('14 species; all 6 ') * italic('Hypsibius') * plain(' candidates')),.045,.125,13.7,col=muted)
 txt('FoldMason',.490,.170,18,font=2,just='center')
 txt('Completed Tamarind job',.490,.137,12.7,col=muted,just='center')
 arr(.331,.402,.147);arr(.578,.644,.147)
 txt('Amino-acid + 3Di alignments',.674,.181,15.2,font=2)
 txt('Interactive report, guide tree, settings and logs',.674,.148,12.5,col=muted)
 txt('Guide tree is not a resolved phylogeny.',.674,.115,12.5,col=muted)
 line(c(.045,.955),c(.080,.080),lwd=.6)
 txt('Native compute: 2 Foldseek jobs + 1 FoldMason job. Existing structures; no new predictions.',.045,.052,11.7,font=2)
 txt('Custom Mol* visualization; thumbnails independently sized. Exact inputs and database details in provenance.',.045,.027,11.3,col=muted)
}
pdf<-file.path(out,'FIGURE_1_WORKBENCH_DISCOVERY.pdf')
png<-file.path(out,'FIGURE_1_WORKBENCH_DISCOVERY.png')
quartz(type='pdf',file=pdf,width=16,height=10,family='Helvetica',bg='white');draw();dev.off()
png(png,width=3200,height=2000,res=200,type='quartz',bg='white');draw();dev.off()
