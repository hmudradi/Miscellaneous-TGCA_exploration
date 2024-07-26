if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("maftools")

 browseVignettes("maftools")
 # Set working directory
 setwd("C:\\Users\\Harini Mudradi\\Downloads\\vanallen-assessment")
 

 # Load the maftools package
 library(maftools)
 ###################################################QUESTION 1###########################################################################
 # Getting a list of MAF files using Sys.glob
 maf_files <- Sys.glob("mafs/Patient*.maf")
 clin_data <- read.table(file="sample-information.tsv", sep="\t", header=TRUE) 
 maf <- merge_mafs(maf_files, verbose = TRUE, clinicalData=clin_data)
 
 
 #Removing empty values from Protein_Change
 maf@data <- maf@data[!is.na(maf@data$Protein_Change) & maf@data$Protein_Change != "", ]
 ####################################################QUESTION 2#############################################################
 # Subsetting for mutations that are not of the Variant Classification "Silent"
 nonsyn_maf <- maf@data[maf@data$Variant_Classification != "Silent", ] # This returns the same num of rows as the data subset since MAF separated silent mutations when merging them together.
 
 #taking a look at the maf files
 plotmafSummary(
   maf,
   rmOutlier = TRUE,
   dashboard = TRUE,
   titvRaw = TRUE,
   log_scale = FALSE,
   addStat = NULL,
   showBarcodes = FALSE,
   fs = 1,
   textSize = 0.8,
   color = NULL,
   titleSize = c(1, 0.8),
   titvColor = NULL,
   top = 15
 )
 ######################################QUESTION 3####################################################################################
# Counting the occurrences of each unique combination of Hugo_Symbol and Protein_Change
mutation_counts <- table(nonsyn_maf$Hugo_Symbol, nonsyn_maf$Protein_Change)
# Converting the table to a data frame
mutation_counts_df <- as.data.frame(mutation_counts)
# Renaming the columns
colnames(mutation_counts_df) <- c('Hugo_Symbol', 'Protein_Change', 'Count')
# Sorting the data frame by Count in descending order
sorted_mutations_df <- mutation_counts_df[order(-mutation_counts_df$Count), ]
# Selecting the top 15 most common mutations(when considering just Gene Name and protein change)
top_mutations_df <- sorted_mutations_df[1:min(15, nrow(sorted_mutations_df)), ]
#print(top_mutations_df)

#After finding most common mutations, taking a look at the top frequency mutation will give us a lot more information about the type of mutation and the gene with the highest mutation.
oncoplot(maf = maf, top = 15, fontSize = 1, genesToIgnore = NULL, sampleOrder = maf@variants.per.sample$Tumor_Sample_Barcode)

#maf.sig = oncodrive(maf = maf, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')

#somaticInteractions(maf = maf, top = 15, pvalue = c(0.05,0.1))


#########################QUESTION 4###############################################################################################
# Clinical enrichment to look for the Enriched gene/genes
response.ce <- clinicalEnrichment(maf=maf, clinicalFeature= "Response" )

#This function uses Fisher's exact test to find genes that are differentially enriched in a group of interest compared to the rest of the sample.
# It calculates the odds ratio to measure the likelihood of observing mutants in the group of interest versus wild-type.

# Significant associations p-value < 0.05
res.sig <- response.ce$groupwise_comparision[p_value < 0.05]
write.csv(res.sig, "enriched_sig.csv", row.names=F)

#################################################QUESTION 5#######################################################################

# Load the required library
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)
library(dplyr)

# Convert 'enrichment_results' to a data frame
enrichment_results <- as.data.frame(res.sig)
# Creating a number of mutants in responders column
res.sig$num_mutated <- ifelse(res.sig$Group1 == "Responder",
                              as.numeric(gsub(" of 25", "", res.sig$n_mutated_group1)),
                              0)

#Scatterplot
ggplot(res.sig, aes(x = num_mutated, y = p_value, color = Group1)) +
  geom_point(position = position_jitter(width = 0.3, height = 0), size = 3, alpha = 0.7) +
  labs(x = "Number of Mutated Patients(out of 25) in Group 1",
       y = "P-value",
       title = "Enrichment Analysis of Mutated Genes") +
  theme_minimal() +
  scale_color_manual(values = c("Non-Responder" = "blue", "Responder" = "red")) +
  geom_text(aes(label = Hugo_Symbol), hjust = -0.1, vjust = 0.5, size = 3)
  
# Since bettering the parameters also didn't help in bettering the scatterplot, plotEnrichmentResults is used to perform the same function with more visual appeal.
# Setting the PDF file name and size
pdf("enrichment_plot.pdf", width = 12, height = 8)

#  plotEnrichmentResults to visualize CLinical ENrichment results. Here we set the parameters pvalue less than equal to 0.05 and odds ratio threshold as 1.
# When the odds ratio is greater than 1 and the p-value is less than 0.05, it indicates that there is a statistically significant association between the two groups(wild and mut) related to each other.
plotEnrichmentResults(
  response.ce,
  pVal = 0.05,
  ORthr = 1,
  featureLvls = NULL,
  cols = NULL,
  annoFontSize = 0.8,
  geneFontSize = 0.8,
  legendFontSize = 0.8,
  showTitle = TRUE
)

# Closing the PDF device
dev.off()

###################################QUESTION 6#########################################################################################3

#pfamDomains show us the enriched proteins for our maf files
pfamprotein = pfamDomains(maf=maf, AACol = 'Protein_Change', top = 15) #The top protein and the right protein in the plot is highly enriched. COG504 and 7tm_ respectively

#all the 7 genes listed here were the enriched genes which was visualized through plotEnrichmentResults function
subset_maf <- subsetMaf(maf = maf, genes = c('ERCC2' , 'AKAP9', 'HERC1', 'HECTD1', 'MACF1', 'MROH2B', 'KMT2C'))
tmb(maf = subset_maf) #tmb is the mutation burden function calculatated for nonsynonymous mutations per megabase

#mafSurvival(maf = maf, genes = 'ERCC2', isTCGA = FALSE, No time given)



# Calculating the sum of values who are responders in num_mutated to find out total number of mutants
total_mutants <- sum(res.sig$num_mutated[res.sig$Group1 == "Responder"]) #46

#Calculating total number of wild types
res.sig$num_total <- as.numeric(gsub(".* of ([0-9]+)$", "\\1", res.sig$n_mutated_group1))
total <- sum(res.sig$num_total[res.sig$Group1 == "Responder"]) #175
total_wildtype = total - total_mutants #129

# Creating a 2x2 table of mutations for ERCC2 gene(most enriched gene). Mutant is 9 , total number of patients is 25 and wildtype patients are 16
mutations_table <- matrix(c(9, 16, 16, 9), nrow = 2, byrow = TRUE)

# Perform Fisher's exact test
fisher_result <- fisher.test(mutations_table)

# Print the result
print(fisher_result)

#The results of Fisher's exact test indicate that the p-value is 0.08874, which is greater than 0.05. 
#This suggests that there is no statistically significant difference in the number of mutations between the mutant and wild-type groups for the enriched gene. 
#The odds ratio estimate of 0.3241952 indicates that the odds of having a mutation in the mutant group are lower than in the wild-type group, but the confidence interval includes 1, indicating that the difference is not statistically significant at the 0.05 level.


