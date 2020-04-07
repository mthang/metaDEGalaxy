library('getopt')
library('data.table')
suppressPackageStartupMessages(library('phyloseq'))
suppressPackageStartupMessages(library('DESeq2'))

Sys.setenv("DISPLAY"=":1")

options(warn= -1)
option_specification = matrix(c(
   'infile','i',2,'character',
   'metafile','m',2,'character',
       'biom','b',2,'character',
     'obsfile','o',2,'character',
     'norm','n',2,'logical',
     'xcolumn','x',2,'numeric',
     'lcolumn','l',2,'numeric',
     'outdir','d',2,'character',
   'htmlfile','h',2,'character'
),byrow=TRUE,ncol=4);


options <- getopt(option_specification)
options(bitmapType="cairo")

matrix.format<-function(x) {
    m<-as.matrix(x[,-1])
    rownames(m)<-x[,1]
    m
}


gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}


tax_col_norm <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_col_extra <- c("None","Kingdom","Phylum","Class","Order","Family","Genus","Species") 

tax_col_norm_otu <- c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_col_extra_otu <- c("OTUID","None","Kingdom","Phylum","Class","Order","Family","Genus","Species")

if (!is.null(options$outdir)) {
  # Create the directory
  dir.create(options$outdir,FALSE)
}

is.biom<-options$biom


pdffile <- gsub("[ ]+", "", paste(options$outdir,"/pdffile.pdf"))
pngfile_net <- gsub("[ ]+", "", paste(options$outdir,"/net.png"))
htmlfile <- gsub("[ ]+", "", paste(options$htmlfile))

if(is.biom=="set_biom"){

galaxy_biom <- import_biom(options$infile)
galaxy_map <- import_qiime_sample_data(options$metafile)

number.of.tax.rank<-length(colnames(tax_table(galaxy_biom)))

	if(number.of.tax.rank == 7){
	    colnames(tax_table(galaxy_biom)) <- tax_col_norm
	} else {
	    colnames(tax_table(galaxy_biom)) <- tax_col_extra
	}

physeq_galaxy <- merge_phyloseq(galaxy_biom,galaxy_map)

} else {

count.table<-read.table(options$infile,header=T,sep="\t",comment.char="",stringsAsFactors = F)
meta.table<-read.table(options$metafile,header=T,sep="\t",comment.char="",stringsAsFactors = F)
tax.table<-read.table(options$obsfile,header=T,sep="\t",comment.char="",stringsAsFactors = F)


colnames(count.table)<-gsub("^X","",colnames(count.table))
colnames(meta.table)<-gsub("^X.","",colnames(meta.table))
colnames(tax.table)<-gsub("^X.","",colnames(tax.table))

count.table.formatted<-matrix.format(count.table)
OTU<-otu_table(count.table.formatted,taxa_are_rows = TRUE)

tax.table.new<-as.data.frame(cbind(tax.table[,1],t(as.data.table(strsplit(tax.table[,2],";")))))


if(length(colnames(tax.table.new)) != length(tax_col_extra_otu)){
colnames(tax.table.new)<-tax_col_norm_otu
}else{
colnames(tax.table.new)<-tax_col_extra_otu
}

tax.table.formatted<-matrix.format(tax.table.new)

TAX<-tax_table(tax.table.formatted)

physeq_galaxy <- phyloseq(OTU, TAX)


galaxy_map<-meta.table

rownames(galaxy_map)<-meta.table[,1]

sampledata<-sample_data(as.data.frame(galaxy_map,row.names=sample_names(galaxy_map),stringsAsFactos=F))

sample_data(physeq_galaxy)<-sampledata

}


x.selectedColumn<-colnames(galaxy_map)[options$xcolumn]
l.selectedColumn<-colnames(galaxy_map)[options$lcolumn]

### normalisation
if(is.null(options$norm) || options$norm =="false"){
suppressMessages(raw.count.deseq2.obj<-phyloseq_to_deseq2(physeq_galaxy,as.formula(paste('~',x.selectedColumn,sep=""))))
geoMeans = apply(counts(raw.count.deseq2.obj), 1, gm_mean)

deseq.obj = estimateSizeFactors(raw.count.deseq2.obj, geoMeans = geoMeans)
deseq.obj.norm<-otu_table(as.matrix(counts(deseq.obj,normalized=T)),taxa_are_rows=TRUE)

otu_table(physeq_galaxy)<-deseq.obj.norm
}


# Produce PDF file

pdf(pdffile); 
plot_net(physeq_galaxy,point_label=x.selectedColumn,color=l.selectedColumn)
garbage<-dev.off();

#Cairo(pngfile_net, type="png", bg="white",pointsize=12,dpi=100,units="in",width=6,height=6)
png(pngfile_net,units="in",width=6,height=6,pointsize=12,res=100,bg="white")
plot_net(physeq_galaxy,point_label=x.selectedColumn,color=l.selectedColumn)
garbage<-dev.off() 

# Produce the HTML file
htmlfile_handle <- file(htmlfile)
html_output = c('<html><body>',
	        '<table align="center>',
		'<tr>',
		'<td valign="middle" style="vertical-align:middle;">',
                '<a href="pdffile.pdf"><img src="net.png"/></a>',
		'</td>',
		'</tr>',
		'</table>',
                '</html></body>');
writeLines(html_output, htmlfile_handle);
close(htmlfile_handle);
