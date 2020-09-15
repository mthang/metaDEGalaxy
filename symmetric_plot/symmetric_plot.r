library('getopt')
library('tidyr')
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('phyloseq'))
suppressPackageStartupMessages(library('DESeq2'))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('data.table'))
Sys.setenv("DISPLAY"=":1")


options(warn= -1)
option_specification = matrix(c(
  'input.data','i',2,'character',
   'meta.data','m',2,'character',
   'obs.data','t',2,'character',
   'record','r',2,'numeric',
   'taxrank','x',2,'character',
	'norm','n',2,'logical',
    'n.column','c',2,'numeric',
     'g.group','g',2,'character',
      'outdir','o',2,'character',
    'htmlfile','h',2,'character'
),byrow=TRUE,ncol=4);


options <- getopt(option_specification);
options(bitmapType="cairo")
 

if (!is.null(options$outdir)) {
  # Create the directory
  dir.create(options$outdir,FALSE)
}


input.table<-read.table(options$input.data,header=T,sep="\t",stringsAsFactors = F)
metadata.table<-read.table(options$meta.data,header=T,sep="\t",stringsAsFactors = F,comment.char="")
obs.table<-read.table(options$obs.data,header=F,sep="\t",stringsAsFactors = F,comment.char="")

colnames(obs.table)<-c("OTUID","taxonomy")


tax_col <- c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_col_extra <- c("OTUID","None","Kingdom","Phylum","Class","Order","Family","Genus","Species")

### remove the leading #sign in column name if the column name begins with number
colnames(input.table)<-gsub("^X","",colnames(input.table))
colnames(metadata.table)<-gsub("^X.","",colnames(metadata.table))


column.name<-colnames(metadata.table)[options$n.column]
#in.column<-options$g.column
in.group<-options$g.group
nrecord<-options$record
ranking<-options$taxrank

### create data frame for group
group.df<-data.frame(group=unlist(strsplit(in.group,",")),stringsAsFactors = F)

### get the number of group
number.of.group<-dim(group.df)[1]

if(number.of.group != 2){
  print(paste("Number of group for comparision is",number.of.group,sep=""))
  quit("yes")
}



group1<-metadata.table[which(metadata.table[,column.name] %in% group.df$group[1]),]$SampleID
group2<-metadata.table[which(metadata.table[,column.name] %in% group.df$group[2]),]$SampleID

sample2group.map <- data.frame(sample=colnames(input.table[,c(group1,group2)]),
                          groups=c(rep(group.df$group[1],length(group1)),
                                   rep(group.df$group[2],length(group2))))


if(options$norm =="false"){
### raw count table
    count.table<-input.table
    rownames(count.table)<-count.table[,1]
    count.table<-count.table[,-1]

    suppressMessages(raw.count.deseq.obj<-DESeqDataSetFromMatrix(countData = count.table,colData=metadata.table, as.formula(paste('~',column.name,sep=""))))

    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
      }

    geoMeans = apply(counts(raw.count.deseq.obj), 1, gm_mean)
    deseq_obj = estimateSizeFactors(raw.count.deseq.obj, geoMeans = geoMeans)
    deseq_obj_norm<-counts(deseq_obj,normalized=T)
} else {
    deseq_obj_norm<-input.table
    rownames(deseq_obj_norm)<-deseq_obj_norm[,1]
    deseq_obj_norm<-deseq_obj_norm[,-1]
}


#first.50.otu.id<-rownames(counts(raw.count.deseq.obj))[1:50]
first.50.otu.id<-rownames(deseq_obj_norm)[1:nrecord]


#filtered.data<-counts(raw.count.deseq.obj)[first.50.otu.id,c(group1,group2)]
filtered.data<-deseq_obj_norm[first.50.otu.id,c(group1,group2)]
filtered.data<-as.data.frame(cbind(OTUID=rownames(filtered.data),filtered.data),stringsAsFactors=F)

nc <- match('OTUID',colnames(filtered.data))
filtered.data[,-nc]<-sapply(filtered.data[,-nc],as.integer)
stopifnot(min(range(filtered.data[,-nc]))>=0)

long <- gather(filtered.data,sample,expr,-OTUID)
suppressMessages((long <- left_join(long, sample2group.map)))
long$expr[long$groups == group.df$group[1]] <- long$expr*-1

sorted.OTU <- rev(sort(unique(long$OTUID)))
long$OTUID <- factor(long$OTUID, levels=sorted.OTU)
long$sample  <- as.factor(long$sample)


tax.table.new<-as.data.frame(cbind(obs.table[,1],t(as.data.table(strsplit(obs.table[,2],";")))))

if(length(colnames(tax.table.new)) != length(tax_col_extra))
{
colnames(tax.table.new)<-tax_col
} else {
colnames(tax.table.new)<-tax_col_extra
}

long<-cbind(long,tax.table.new[match(long$OTUID,tax.table.new$OTUID),-1])

comparison<-paste(column.name,paste(group.df$group[1],group.df$group[2],sep="-"),sep=" ")


p<-ggplot(long,aes(x=reorder(OTUID,expr),y=expr, fill=groups)) +
  geom_bar(stat='identity') + theme_bw() + 
  xlab("OTU ID") + 
  labs(title=comparison) + facet_grid( as.formula(paste(ranking,"~ .",sep="")),scales = "free", space = "free" ) + theme(strip.text.y = element_text(angle = 0)) +
  coord_flip()
q<-ggplot_build(p)
q$layout$panel_ranges[[1]]$x.labels <- gsub("-","",q$layout$panel_ranges[[1]]$x.labels)
#Reassemble the plot using ggplot_gtable()
q<-ggplot_gtable(q)




pdffile <- gsub("[ ]+", "", paste(options$outdir,"/symmetric.pdf"))
pngfile_symmetric <- gsub("[ ]+", "", paste(options$outdir,"/symmetric.png"))
htmlfile <- gsub("[ ]+", "", paste(options$htmlfile))


# Produce PDF file
pdf(pdffile);
plot(q)
garbage<-dev.off();

#png('richness.png')
bitmap(pngfile_symmetric,"png16m",height=10,width=10,res=100)
plot(q)
garbage<-dev.off()

# Produce the HTML file
 htmlfile_handle <- file(htmlfile)
 html_output = c('<html><body>',
 	        	 '<table align="center">',
 		         '<tr>',
 		         '<td valign="middle" style="vertical-align:middle;">',
                 '<a href="pdffile.pdf"><img src="symmetric.png"/></a>',
 		         '</td>',
 		         '</tr>',
 		         '</table>',
                 '</html></body>');
 writeLines(html_output, htmlfile_handle);
 close(htmlfile_handle);
