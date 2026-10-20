script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);root<-normalizePath(file.path(dirname(script),'../..'))
out<-file.path(root,'research/release');qa<-file.path(out,'qa');dir.create(qa,showWarnings=FALSE,recursive=TRUE)
manifest<-jsonlite::fromJSON(file.path(out,'figures_manifest.json'),simplifyVector=FALSE)
ids<-c(sort(grep('^F',names(manifest$figures),value=TRUE)),sort(grep('^S',names(manifest$figures),value=TRUE)))
sha<-function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1]
pdfs<-file.path(out,'figures',paste0(ids,'.pdf'))
stopifnot(system2('pdfunite',shQuote(c(pdfs,file.path(out,'FIGURE_ATLAS.pdf'))))==0)
records<-list();images<-list()
for(id in ids){
  prefix<-file.path(qa,id)
  stopifnot(system2('pdftoppm',c('-scale-to','1450','-singlefile','-png',shQuote(file.path(out,'figures',paste0(id,'.pdf'))),shQuote(prefix)))==0)
  im<-png::readPNG(paste0(prefix,'.png'));images[[id]]<-im
  gray<-im;light<-.2126*im[,,1]+.7152*im[,,2]+.0722*im[,,3]
  for(channel in 1:3)gray[,,channel]<-light
  png::writePNG(gray,paste0(prefix,'_gray.png'))
  rec<-list(color_preview_sha256=sha(paste0(prefix,'.png')),grayscale_sha256=sha(paste0(prefix,'_gray.png')))
  if(id %in% c('F03_single_cell_atlas','F04_expression_evidence','F06_clinical_associations','F07_functional_context')){
    hex<-grDevices::rgb(as.vector(im[,,1]),as.vector(im[,,2]),as.vector(im[,,3]));palette<-unique(hex)
    converted<-grDevices::col2rgb(colorspace::deutan(palette,severity=1))/255;index<-match(hex,palette);cvd<-im
    for(channel in 1:3)cvd[,,channel]<-matrix(converted[channel,index],nrow=dim(im)[1],ncol=dim(im)[2])
    png::writePNG(cvd,paste0(prefix,'_deutan.png'));rec$deutan_sha256<-sha(paste0(prefix,'_deutan.png'))
  }
  records[[id]]<-rec
}
png(file.path(out,'OVERVIEW.png'),width=2400,height=1740,res=220,type='cairo',bg='#F3F5F7')
grid::grid.newpage();grid::grid.text('Psoriasis | figures and evidence',x=.025,y=.97,just='left',gp=grid::gpar(fontfamily='Liberation Sans',fontface='bold',fontsize=19,col='#243746'))
short<-c('Genetic reconstruction','Shared-outcome sensitivity','Single-cell atlas','Expression evidence','Regional colocalization','Clinical associations','Functional context','All registered loci')
for(i in seq_along(ids)){
 col<-(i-1)%%4;row<-(i-1)%/%4;cx<-(col+.5)/4;cy<-.715-row*.443
 grid::grid.rect(x=cx,y=cy,width=.239,height=.408,gp=grid::gpar(fill='white',col='#D8E0E5',lwd=.6))
 im<-images[[ids[i]]];scale<-min(550/dim(im)[2],655/dim(im)[1]);w<-dim(im)[2]*scale/2400;h<-dim(im)[1]*scale/1740
 grid::grid.raster(im,x=cx,y=cy+.004,width=w,height=h,interpolate=TRUE)
 grid::grid.text(paste0(if(i<8)paste0('Figure ',i) else 'Supplement 1','  |  ',short[i]),x=cx,y=cy-.218,gp=grid::gpar(fontfamily='Liberation Sans',fontsize=8.2,col='#243746'))
}
dev.off()
capture.output(sessionInfo(),file=file.path(qa,'sessionInfo.txt'))
record<-list(status='review_artifacts_generated',figures=length(ids),script_sha256=sha(script),
             inputs=list(figures_manifest_sha256=sha(file.path(out,'figures_manifest.json'))),previews=records,
             outputs=list(atlas_sha256=sha(file.path(out,'FIGURE_ATLAS.pdf')),overview_sha256=sha(file.path(out,'OVERVIEW.png'))),
             scope='Rendering and color simulation only; human visual review and statistical checks are recorded separately')
jsonlite::write_json(record,file.path(qa,'render_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
