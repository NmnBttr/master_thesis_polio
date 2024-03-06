# clear workspace
rm(list = ls())

# set working directory
setwd("path/to/heatmap_tabular")

# library install and import

list.of.packages <- c("pheatmap", "RColorBrewer", "dplyr", "fs", "tools", "stringr", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

# options 
name_plot <- "Variant_frequency"               ## name your plot
frequency <- 0.1                   ## adjust the variant frequency
display_as <- T                     ## display also AS changes 
percent <- F                        ## should variants be shown as percent?
rounded <- F                        ## should rounded values be displayed in the map
multiplier_width <- 0.5            ## pdf weight multiplier - adjust if needed
multiplier_height <- 0.6            ## pdf height multiplier - adjust if needed
color_gene_annotation <- c("Set3")  ## adjust color of gene annotation ("Set2" or "Paired")
date <- F                           ## do the file names contain a date (format: dd.mm.yyyy) and you want to sort
number <- T                         ## do the file names contain a number and you want to sort
clustering <- F                     ## should the samples be clustered?
clustering_method <- "ward.D2"      ## what clustering method -> see ?hclust for further information
number_of_clusters <- 3             ## do you assume a particular amount of clusters?



# data preparation

## import data --- pattern = "SiB159_*""*.tabular"
files<-list.files(path = ".",pattern = "*.tabular", recursive = F ,full.names = TRUE)

## only in case it contains a date:
if(date) {
  files<-files[order(as.Date(regmatches(files, regexpr("\\d{2}[[:punct:]]\\d{2}[[:punct:]]\\d{4}", files)), format="%d.%m.%y"))]
}

## create data.frame for first dataset
variants<-read.table(files[1], header = T, sep = "\t",stringsAsFactors=FALSE, 
                     colClasses = c("character"), fill = T)
final<- paste0(variants$POS, variants$ALT)
final<- data.frame(cbind(final, variants$AF))
Sample_number<-path_file(files[1])
Sample_number = file_path_sans_ext(Sample_number)
colnames(final) <- c("Mutation", Sample_number)


## loop over all other datasets and join by Mutation
for(i in 2:length(files)) {
  variants<-read.table(files[i], header = T, sep = "\t",stringsAsFactors=FALSE, 
                       colClasses =c("character"), fill = T, skipNul = TRUE)
  varpos<- paste0(variants$POS, variants$ALT)
  varpos<- cbind(varpos, variants$AF)
  Sample_number<-path_file(files[i])
  Sample_number = file_path_sans_ext(Sample_number)
  colnames(varpos) <- c("Mutation", Sample_number)
  final <- full_join(final, varpos, by="Mutation", copy=T)
}

## sort
final<-final[str_order(final$Mutation, numeric = T),]

## set row names
# rownames(diasyhoras2) <- make.names(diasyhoras[,1], unique = TRUE)
row.names(final)<-final$Mutation
final<-final[,-1]


## order samples after number in name
if(number){final<-final[,str_order(colnames(final), numeric = T)]}

## adjust the variant frequency:
#final[final < frequency] <- NA

final <- final[rowSums(is.na(final)) !=ncol(final), ]
final <- t(final)
final[is.na(final)] <- 0
class(final) <- "numeric"

# convert values to percent
if(percent){final<-final*100}

# show rounded values in pheatmap
if(rounded){final_rounded <- round(final, 0)} else {final_rounded <- F}

# add annotations for Polio 1. For Polio 2 and 3 adjust protein positions

## readout annotations
ann_final <- read.table(files[1], header = T, sep = "\t", fill = T) 
# ann_final <- data.frame(paste0(ann_final$POS, ann_final$ALT, if(display_as){paste0(" ", "(", ann_final$EFF....AA, ")")}), ann_final$EFF....EFFECT, ann_final$EFF....GENE)
# colnames(ann_final) <- c("position","effect", "gene")
ann_pos <- data.frame(ann_final$POS)
colnames(ann_pos) <- c("position")
protein_pos = list()
for (i in 1:nrow(ann_pos)) {
  if(ann_pos$position[i] >= 1 & ann_pos$position[i] < 743) {
    protein_pos = append(protein_pos, "5'-UTR")
  } else if (ann_pos$position[i] >= 743 & ann_pos$position[i] <= 949) {
    protein_pos = append(protein_pos,"VP4")
  } else if (ann_pos$position[i] >= 950 & ann_pos$position[i] <= 1765) {
    protein_pos = append(protein_pos,"VP2")
  } else if (ann_pos$position[i] >= 1766 & ann_pos$position[i] <= 2479) {
    protein_pos = append(protein_pos,"VP3")
  } else if (ann_pos$position[i] >= 2480 & ann_pos$position[i] <= 3385) {
    protein_pos = append(protein_pos,"VP1")
  } else if (ann_pos$position[i] >= 3386 & ann_pos$position[i] <= 3832) {
    protein_pos = append(protein_pos,"2A")
  } else if (ann_pos$position[i] >= 3833 & ann_pos$position[i] <= 4123) {
    protein_pos = append(protein_pos,"2B")
  } else if (ann_pos$position[i] >= 4124 & ann_pos$position[i] <= 5110) {
    protein_pos = append(protein_pos,"2C")
  } else if (ann_pos$position[i] >= 5111 & ann_pos$position[i] <= 5371) {
    protein_pos = append(protein_pos,"3A")
  } else if (ann_pos$position[i] >= 5372 & ann_pos$position[i] <= 5437) {
    protein_pos = append(protein_pos,"3B")
  } else if (ann_pos$position[i] >= 5438 & ann_pos$position[i] <= 5986) {
    protein_pos = append(protein_pos,"3C")
  } else if (ann_pos$position[i] >= 5987 & ann_pos$position[i] <= 7369) {
    protein_pos = append(protein_pos,"3D")
  } else {
    protein_pos = append(protein_pos,"undefined region")
  }
}
ann_final <- data.frame(paste0(ann_final$POS, ann_final$ALT),unlist(protein_pos))
colnames(ann_final) <- c("position", "gene")

for(i in 2:length(files)) {
  ann <- read.table(files[i], header = T, sep = "\t", fill = T)
  ann_pos <- data.frame(ann$POS)
  colnames(ann_pos) <- c("position")
  protein_pos = list()
  for (i in 1:nrow(ann_pos)) {
    if (ann_pos$position[i] >= 1 & ann_pos$position[i] < 743) {
      protein_pos = append(protein_pos, "5'-UTR")
    } else if (ann_pos$position[i] >= 743 & ann_pos$position[i] <= 949) {
      protein_pos = append(protein_pos,"VP4")
    } else if (ann_pos$position[i] >= 950 & ann_pos$position[i] <= 1765) {
      protein_pos = append(protein_pos,"VP2")
    } else if (ann_pos$position[i] >= 1766 & ann_pos$position[i] <= 2479) {
      protein_pos = append(protein_pos,"VP3")
    } else if (ann_pos$position[i] >= 2480 & ann_pos$position[i] <= 3385) {
      protein_pos = append(protein_pos,"VP1")
    } else if (ann_pos$position[i] >= 3386 & ann_pos$position[i] <= 3832) {
      protein_pos = append(protein_pos,"2A")
    } else if (ann_pos$position[i] >= 3833 & ann_pos$position[i] <= 4123) {
      protein_pos = append(protein_pos,"2B")
    } else if (ann_pos$position[i] >= 4124 & ann_pos$position[i] <= 5110) {
      protein_pos = append(protein_pos,"2C")
    } else if (ann_pos$position[i] >= 5111 & ann_pos$position[i] <= 5371) {
      protein_pos = append(protein_pos,"3A")
    } else if (ann_pos$position[i] >= 5372 & ann_pos$position[i] <= 5437) {
      protein_pos = append(protein_pos,"3B")
    } else if (ann_pos$position[i] >= 5438 & ann_pos$position[i] <= 5986) {
      protein_pos = append(protein_pos,"3C")
    } else if (ann_pos$position[i] >= 5987 & ann_pos$position[i] <= 7369) {
      protein_pos = append(protein_pos,"3D")
    } else {
      protein_pos = append(protein_pos,"undefined region")
    }
  }
  ann <- data.frame(paste0(ann$POS, ann$ALT), unlist(protein_pos))
  colnames(ann) <- c("position", "gene")
  ann_final <- rbind(ann_final, ann)
}

ann_final<-data.frame(unique(ann_final))

## apply frequency filter
ann_final <- ann_final[ann_final$position %in% colnames(final),]
## sort
ann_final <- ann_final[str_order(ann_final$position, numeric = T),]
## reshape for annotation_col in pheatmap
ann_final <- data.frame(ann_final[,-1, drop = FALSE], row.names = ann_final$position)

# automatically determine gaps for the heatmap
gap_vector <- c()

for (i in 2:length(ann_final$gene)){
  if (ann_final$gene[i]!=ann_final$gene[i-1]){
    gap_vector <- c(gap_vector, i-1)
  }
}

# colormanagement

## colormangment heatmap
my_colors <- colorRampPalette(c("grey93","brown","black"))

## colormanagement annotations (genes)
count <- length(unique(ann_final$gene))
gene_color <- c(brewer.pal(color_gene_annotation, n=count))
names(gene_color) = unique(ann_final$gene)

color_list <- list(gene_color = gene_color)
names(color_list) <- c("gene")


# visualize heatmap

## The pdf scales with the number of variants and samples
# pdf(paste0(name_plot, "(", frequency, ")", ".pdf"), width = multiplier_width*ncol(final), height = multiplier_height*nrow(final))

pheatmap(final, 
         color = my_colors(100),
         #cellwidth = 20,
         #cellheight = 20,
         clustering_method = clustering_method,
         cluster_rows = clustering, 
         cluster_cols = F,
         cutree_rows = number_of_clusters,
         annotation_col = ann_final,
         gaps_col = gap_vector,
         display_numbers = final_rounded
         )

dev.off()

## file checker - checks if files have the expected column count of 12
for (i in 1:length(files)) {
  x <- fread(files[i], header = T, sep = "\t", fill = T)
  if (ncol(x) == 6){
    print(paste(files[i], "...............", "ok"))
  } else {
    print(paste(files[i], "...............", "error")) 
    
  }
  
}
