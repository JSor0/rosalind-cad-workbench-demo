#!/usr/bin/env Rscript
# Existing evidence only. Grouping correction requested by the user.
library(grid)
library(png)
a<-commandArgs(trailingOnly=TRUE)
pkg<-normalizePath(if(length(a))a[1] else '.')
out<-if(length(a)>1)normalizePath(a[2]) else pkg
mol<-file.path(pkg,'assets/hero')
d<-read.delim(file.path(pkg,'source_data/figure2_alignment_barcode.tsv'),stringsAsFactors=FALSE)
stopifnot(nrow(d)==211,sum(d$class=='identical')==32,sum(d$class=='different')==139,sum(d$class=='gap')==40)
ink<-'#1D3238';teal<-'#247F82';orange<-'#C8792C';muted<-'#5C6C72';gray<-'#CDD3D5';pale<-'#E9EDEE';red<-'#B35D57';wc<-'#59646A'
tx<-function(s,x,y,z=15,col=ink,face=1,j='left')grid.text(s,x,y,just=j,gp=gpar(fontfamily='Helvetica',fontsize=z,col=col,fontface=face))
ln<-function(x,y,col=gray,lwd=.7)grid.lines(x,y,gp=gpar(col=col,lwd=lwd))
im<-function(name,x,y,w,h,bounds){
 p<-readPNG(file.path(mol,name));p<-p[(bounds[2]+1):bounds[4],(bounds[1]+1):bounds[3],,drop=FALSE]
 ar<-dim(p)[2]/dim(p)[1]
 if(w*16/(h*11)>ar){ww<-h*11*ar/16;hh<-h}else{ww<-w;hh<-w*16/ar/11}
 grid.raster(p,x=x,y=y,width=ww,height=hh,interpolate=TRUE)
}
icon_dir<-file.path(pkg,'assets/preserved_v7/approved_icons')
icons<-list('Mouse CAD'=readPNG(file.path(icon_dir,'mouse_teal.png'))[371:960,91:1165,,drop=FALSE], 'X5S1'=readPNG(file.path(icon_dir,'tardigrade_active.png'))[371:940,111:1155,,drop=FALSE])
arch<-function(name,y,L,lo,hi,col){
 a<-icons[[name]]; ar<-dim(a)[2]/dim(a)[1]
 grid.raster(a,x=.062,y=y,width=.038,height=.038*16/11/ar,interpolate=TRUE)
 f<-function(n).162+n/423*.598
 tx(name,.090,y,14,col,2)
 grid.rect(x=(f(0)+f(L))/2,y=y,width=f(L)-f(0),height=.026,gp=gpar(fill=pale,col=NA))
 grid.rect(x=(f(lo-1)+f(hi))/2,y=y,width=f(hi)-f(lo-1),height=.026,gp=gpar(fill=col,col=NA))
 tx(paste0(if(name=='X5S1')'C3-like' else 'C3',': ',lo,'-',hi),(f(lo-1)+f(hi))/2,y,11.5,'white',2,'center')
 tx(paste0(L,' aa'),f(L)+.010,y,12,muted)
}
pair<-function(mouse,x5,y,difference=FALSE){
 tx(mouse,.702,y,14,if(difference)red else teal,if(difference)2 else 1)
 ln(c(.753,.784),c(y,y),if(difference)'#D7B1AE' else gray,.65)
 tx(x5,.800,y,14,if(difference)red else orange,if(difference)2 else 1)
}
draw<-function(){
 grid.newpage()
 tx('A CAD-like core in a remodeled protein',.045,.960,24,ink,2)
 tx('Mouse CAD and X5S1 share a fold despite sparse sequence identity.',.045,.918,15.5,muted)
 # Full-length context belongs to the protein identities before the core detail.
 arch('Mouse CAD',.861,344,132,328,teal)
 arch('X5S1',.820,423,204,403,orange)
 tx('Same amino-acid scale',.822,.861,11.5,muted)
 tx('Core spans highlighted',.822,.820,11.5,muted)
 tx('No canonical CIDE-N/C1 domain is reported for X5S1; C2-like elements remain possible.',.162,.781,11.5,muted)
 # One legend owns one shared picture. The middle gutter stays empty.
 tx('a  Fixed core superposition',.047,.737,17,ink,2)
 tx('Mouse CAD',.047,.706,13,teal,2)
 tx('experimental 1V0D, A132-328',.145,.706,12,muted)
 tx('X5S1',.047,.681,13,orange,2)
 tx('predicted monomer, A204-403',.145,.681,12,muted)
 im('hero_overlay.png',.338,.416,.583,.494,c(100,255,1670,1365))
 # The role of the right table is explicit and independent of the WC24 panel.
 tx('b  Mouse-X5S1 residue pairs',.691,.737,16,ink,2)
 tx('Fixed mapping; predicted geometry uncertain',.691,.707,11.5,muted)
 tx('Mouse',.702,.669,12.5,teal,2);tx('X5S1',.800,.669,12.5,orange,2)
 tx('Mouse Zn-site positions',.702,.637,12,muted)
 for(i in 1:4)pair(c('C229','C238','H242','C307')[i],c('C296','C311','H316','C392')[i],c(.609,.583,.557,.531)[i])
 tx('Mouse catalytic positions',.702,.490,12,muted)
 for(i in 1:3)pair(c('D262','H263','H308')[i],c('D337','H338','H393')[i],c(.462,.436,.410)[i])
 tx('Mapped nonconservative difference',.702,.370,11.5,red)
 pair('N299','K384',.341,TRUE)
 # WC24 is a separate comparison, not the owner of the table above.
 tx('c  Same-species comparator',.691,.283,16,ink,2)
 tx('WC24',.702,.254,13,wc,2)
 tx(expression(italic('H. exemplaris')),.766,.254,12,muted)
 im('wc24_overlay.png',.818,.153,.248,.180,c(111,490,1655,1354))
 tx('Mouse CAD',.702,.046,11.5,teal,2)
 tx('experimental',.795,.046,11.5,muted)
 tx('WC24 C3',.702,.021,11.5,wc,2)
 tx('predicted',.795,.021,11.5,muted)
 # Pairwise identity bar stays adjacent to its core view and preserves all gaps.
 tx('Sparse identity in the returned Foldseek alignment',.047,.155,14,ink,2)
 bx<-.047;bw<-.582;by<-.123;bh<-.019
 cls<-c(identical=ink,different='#CCD3D5',gap='white')
 grid.rect(x=bx+(seq_len(nrow(d))-.5)/nrow(d)*bw,y=by,width=bw/nrow(d),height=bh,gp=gpar(fill=unname(cls[d$class]),col=NA))
 grid.rect(x=bx+bw/2,y=by,width=bw,height=bh,gp=gpar(fill=NA,col=gray,lwd=.6))
 tx('1',bx,.103,9.5,muted);tx('211 columns',bx+bw,.103,9.5,muted,j='right')
 for(i in 1:3){xx<-.047+(i-1)*.174;grid.rect(x=xx+.004,y=.079,width=.008,height=.011,gp=gpar(fill=cls[i],col=gray,lwd=.5));tx(c('Identical pair','Different pair','Gap column')[i],xx+.016,.079,11,muted)}
 tx('15.1% reported identity. Forward search alignment; preserved columns and gaps.',.047,.052,11,muted)
 tx('Custom Mol* view; fixed coordinates. Residue pairs follow the saved pairwise mapping.',.047,.020,10.5,muted)
}
stopifnot(!file.exists(file.path(out,'FIGURE_2_CAD_X5S1_HERO.pdf')))
quartz(file=file.path(out,'FIGURE_2_CAD_X5S1_HERO.pdf'),type='pdf',width=16,height=11,family='Helvetica',pointsize=15,bg='white');draw();dev.off()
png(file.path(out,'FIGURE_2_CAD_X5S1_HERO.png'),width=3200,height=2200,res=200,type='quartz',bg='white');draw();dev.off()

