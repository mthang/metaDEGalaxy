library('getopt')
suppressPackageStartupMessages(library('phyloseq'))
Sys.setenv("DISPLAY"=":1")

options(warn=-1)
option_specification = matrix(c(
   'biomfile','b',2,'character',
   'metafile','m',2,'character',
     'xcolumn','x',2,'numeric',
     'lcolumn','l',2,'numeric',
     'taxonomy','t',2,'character',
     'outdir','o',2,'character',
   'htmlfile','h',2,'character'
),byrow=TRUE,ncol=4);


options <- getopt(option_specification);
options(bitmapType="cairo")
 

if (!is.null(options$outdir)) {
  # Create the directory
  dir.create(options$outdir,FALSE)
}


galaxy_biom <- import_biom(options$biomfile)
galaxy_map <- import_qiime_sample_data(options$metafile)
tax_col_norm <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_col_extra <- c("None","Kingdom","Phylum","Class","Order","Family","Genus","Species")

number.of.tax.rank<-length(colnames(tax_table(galaxy_biom)))

if( number.of.tax.rank == 7){
colnames(tax_table(galaxy_biom)) <- tax_col_norm
}else{
colnames(tax_table(galaxy_biom)) <- tax_col_extra
}


AIP_galaxy <- merge_phyloseq(galaxy_biom,galaxy_map)


x.selectedColumn<-colnames(galaxy_map)[options$xcolumn]
l.selectedColumn<-colnames(galaxy_map)[options$lcolumn]

tax.selected<-options$taxonomy

pdffile <- gsub("[ ]+", "", paste(options$outdir,"/kingdom.pdf"))
pngfile_kingdom <- gsub("[ ]+", "", paste(options$outdir,"/kingdom.png"))
htmlfile <- gsub("[ ]+", "", paste(options$htmlfile))


# Produce PDF file
pdf(pdffile);
plot_bar(AIP_galaxy,x=x.selectedColumn,facet_grid = paste('~', tax.selected),fill=l.selectedColumn)
garbage<-dev.off();

#png('kingdom.png')
bitmap(pngfile_kingdom,"png16m")
plot_bar(AIP_galaxy,x=x.selectedColumn,facet_grid = paste('~', tax.selected),fill=l.selectedColumn)
garbage<-dev.off()

# Produce the HTML file
htmlfile_handle <- file(htmlfile)
html_output = c('<html><body>',
	            '<table align="center">',
		        '<tr>',
		        '<td valign="middle" style="vertical-align:middle;">',
                '<a href="kingdom.pdf"><img src="kingdom.png"/></a>',
		        '</td>',
		        '</tr>',
		        '</table>',
                '</html></body>');
writeLines(html_output, htmlfile_handle);
close(htmlfile_handle);
