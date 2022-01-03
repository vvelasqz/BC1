#start by getting data and compiling A,B,C samples into one
#sample descriptions
# Sample 1: Uninfected, beads, no antibody (Negative control)
# Sample 3: BC1, beads,  antibody
# Sample 4: Uninfected, beads, antibody
# Sample 5: Infected, beads, antibody
# Sample 7: Uninfected, beads, BC1, antibody
# Sample 8: Infected, beads, no antibody (Negative control)
# Sample 9: BC1, beads, no antibody (Negative control)

library(UniprotR)
file_list<- list.files("redata/")
samples<-list()
names<- substr(file_list, 6, 6)
remove_spp<- c("HUMAN", "PIG", "BOVIN", "REF_HEVBR", "SHEEP")

for(i in 1:length(file_list)){
  sample1<- read.table(paste0("redata/",file_list[[i]]), fill=T, sep="\t", row.names = NULL)
  colnames(sample1)[1:10]<- colnames(sample1)[2:11]
  sample1<- sample1[sample1$Unused>0 & !grepl(paste0(remove_spp, collapse = "|"), sample1$Accession),1:10]
  sample1$Uniprot_name<- sapply(sample1$Accession, FUN=function(x)strsplit(x, "\\|")[[1]][2])
  sample1[sample1$Accession=="BC1", "Uniprot_name"]<-"P14981"
  
  if(!names[[i]] %in% names(samples)){
    samples[[names[[i]]]]<- sample1
  }else{
    samples[[names[[i]]]]<- rbind(samples[[names[[i]]]], sample1)
  }
}

saveRDS(samples, "samples.RDS")
sample1<- samples[[7]]
samples<- readRDS("samples.RDS")
library(tidyverse)
samples_filtered<-list()
bait_numbers<- c("Uninfected_noAB", "Bc1_AB", "Uninfected_AB", "Infected_AB", "Uninfected_BC1_AB", "Infected_noAB", "Bc1_noAB")
names(bait_numbers)<- names(samples)
for(i in names(samples)){
  sample1<- samples[[i]] %>% group_by(Uniprot_name) %>% summarise(Counts=sum(Peptides.95..))
  sample1$length<- GetSeqLength(sample1$Uniprot_name)$Length
  sample1$bait<- bait_numbers[i]
  samples_filtered[[bait_numbers[i]]]<- sample1
  print(i)
}

saveRDS(samples_filtered, "samples_filtered.RDS")

#make saint files
#Bait file: IP name, bait name, test and negative control purifications (T = test, C = control).  
#Prey file: prey (protein) name, prey protein length, and prey gene name. 
#Interaction file: IP name, bait name, prey name, Counts. prey name as the first col of the prey file.  
#Interactions with zero counts must be removed from the file.
#No col names

samples_filtered<- readRDS("samples_filtered.RDS")
samples_filtered_all<- bind_rows(samples_filtered)
bait_file<- data.frame(IP=paste0(c("Uninfected_noAB", "Bc1_AB", "Uninfected_AB", "Infected_AB", "Uninfected_BC1_AB", "Infected_noAB"), rep(c("-1","-2"), each=6)),
                       bait=rep(c("Uninfected_noAB", "Bc1_AB", "Uninfected_AB", "Infected_AB", "Uninfected_BC1_AB", "Infected_noAB"), 2),
                       type=rep(c("C", "T", "C", "T", "T", "C"), 2))

write.table(bait_file, file="SAINTexpress_files/bait_file.dat",quote = FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

prey_file<- data.frame(unique(samples_filtered_all[,c(1,3)]))
write.table(prey_file, file="SAINTexpress_files/prey_file.dat",quote = FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

interaction_file<- data.frame(IP=paste0(samples_filtered_all$bait, rep(c("-1","-2"), each=203)),
                              bait=rep(samples_filtered_all$bait, 2),
                              prey=rep(samples_filtered_all$Uniprot_name, 2),
                              counts=rep(samples_filtered_all$Counts, 2))

write.table(interaction_file, file="SAINTexpress_files/interaction_file.dat",quote = FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

#Running SAINT
# copy files to /SAINTexpress_v3.6.3/Precompiled_binaries/Linux64
# dos2unix interaction_file.dat
# dos2unix prey_file.dat
# dos2unix bait_file.dat
# ./SAINTexpress-spc interaction_file.dat prey_file.dat bait_file.dat

#after running SAINT I will process the files to report
library(UniprotR)
saint_results<- read.table("SAINTexpress_files/list.txt", sep = '\t',header = TRUE, quote = "\"", fill=TRUE)
a<- GetNamesTaxa(saint_results$PreyGene)
saint_results<- cbind(saint_results, a)
saint_results$Protein.names<-str_replace_all(saint_results$Protein.names, pattern = ",", replacement = ";")
saint_results_summary<- saint_results[,c(1,2,6,13,15,16,18,19,22,26)]
write.csv(saint_results, file="SAINTexpress_files/saint_results.csv",quote = FALSE, row.names = FALSE)
write.csv(saint_results_summary, file="SAINTexpress_files/saint_results_summary.csv",quote = FALSE, row.names = FALSE)


##########################################################################################
#make Y2H-SCORES files for 3 baits, uninfected, infected and BC1, noAB will be non-selected sample
library(tidyverse)

save_raw_counts<- function(sample_tables, bait, prey_list, num_reps= 3, path, mul_counts=T){
  table_s<- sample_tables[[paste0(bait,"_AB")]][,c("Uniprot_name", "Counts")]
  table_n<- sample_tables[[paste0(bait,"_noAB")]][,c("Uniprot_name", "Counts")]
  table<- data.frame(prey=prey_list, s_col=0, ns_col=0)
  table[match(table_s$Uniprot_name, table$prey), "s_col"]<- table_s$Counts
  table[match(table_n$Uniprot_name, table$prey), "ns_col"]<- table_n$Counts
  
  raw_counts<-matrix(c(rep(table[,"s_col"], num_reps), rep(table[,"ns_col"], num_reps)), nrow = length(prey_list), ncol = 2*num_reps)
  raw_counts<- round(raw_counts)
  for(col in 1:ncol(raw_counts)){
    if(mul_counts){
      raw_counts[,col]<- 1000*raw_counts[,col] + round(sample(1:num_reps, 1))
    }
  }
  rownames(raw_counts)<- table$prey
  write.table(raw_counts, paste0(path,strsplit(path, "/")[[1]][length(strsplit(path, "/")[[1]])],"_salmon_counts.matrix"), row.names = T, col.names = F, sep = "\t", quote = F)
  
}

make_input_arguments_file<- function(draft_table, bait_name, num_reps){
  draft_table[2,"argument_value"]<- paste(paste0("/", bait_name, "R", 1:num_reps, "S.fastq"), collapse = ";")
  draft_table[4,"argument_value"]<- paste(paste0("/", bait_name, "R", 1:num_reps, "N.fastq"), collapse = ";")
  draft_table[7,"argument_value"]<- paste0("formated_files_for_scores/", bait_name)
  write.csv(draft_table, paste0("Y2H-SCORES_files/formated_files_for_scores/", bait_name, "_input_arguments.csv"), row.names = F, quote = F)
}

draft_table<- read.csv("Y2H_draft_input_arguments.csv")
dir_path<- "/Y2H-SCORES_files/formated_files_for_scores/"
bait_list<- c("Uninfected_BC1", "Bc1", "Infected", "Uninfected")
samples_filtered<- readRDS("samples_filtered.RDS")
samples_filtered[["Uninfected_BC1_noAB"]]<- samples_filtered[["Uninfected_noAB"]]
samples_filtered[["BC1_noAB"]]<- samples_filtered[["Uninfected_noAB"]]
samples_filtered_all<- bind_rows(samples_filtered)
prey_list<- unique(samples_filtered_all$Uniprot_name)

fofn<- list()
for(bait in bait_list){
  bait_folder<- paste0("Y2H-SCORES_files/formated_files_for_scores/", bait, "/")
  dir.create(bait_folder)
  save_raw_counts(sample_tables=samples_filtered, bait=bait, prey_list, num_reps= 3, path= bait_folder)
  make_input_arguments_file(draft_table, bait_name=bait, num_reps=3)
  fofn[[bait]]<- paste0("formated_files_for_scores/", bait, "_input_arguments.csv")
}
fofn_table<- t(bind_rows(fofn)) 
write.table(fofn_table, "Y2H-SCORES_files/fofn_for_compute_scores_busisiwe.txt", 
            col.names = F, row.names = F, quote = F)

#run Y2H-SCORES
#Rscript run_scores.R --fofn "fofn_for_compute_scores_busisiwe.txt" --out_dir "output_busisiwe/" --spec_p_val 1 --spec_fold_change 0 --enrich_p_val 1 --enrich_fold_change 0 --normalized T

##after running Y2H-SCORES I will process the files to report
library(UniprotR)
library(stringr)
library(tidyverse)
y2h_scores_results<- read.csv("Y2H-SCORES_files/output_busisiwe/Total_scores.csv",header = TRUE)
b<- GetNamesTaxa(y2h_scores_results$prey)
y2h_scores_results<- cbind(y2h_scores_results, b)
y2h_scores_results$bait<- paste0(y2h_scores_results$bait,"_AB")
y2h_scores_results$Protein.names<-str_replace_all(y2h_scores_results$Protein.names, pattern = ",", replacement = ";")
y2h_scores_results$Borda_quantile<- sapply(y2h_scores_results$Borda_scores, FUN=function(x)1-(sum(y2h_scores_results$Borda_scores<x)/length(y2h_scores_results$Borda_scores)))
y2h_scores_results_summary<- y2h_scores_results[,c(2,1,3:6,19,7,8,11,15)]
y2h_scores_results_summary_filtered<- y2h_scores_results_summary %>% group_by(prey) %>% slice(which.max(Borda_scores))
y2h_scores_results_summary_filtered<- y2h_scores_results_summary_filtered[y2h_scores_results_summary_filtered$bait!="Uninfected_AB",]
y2h_scores_results_summary_filtered<- y2h_scores_results_summary_filtered[y2h_scores_results_summary_filtered$bait!="Uninfected_AB",]

write.csv(y2h_scores_results, file="Y2H-SCORES_files/y2h_scores_results.csv",quote = FALSE, row.names = FALSE)
write.csv(y2h_scores_results_summary, file="Y2H-SCORES_files/y2h_scores_results_summary.csv",quote = FALSE, row.names = FALSE)
write.csv(y2h_scores_results_summary_filtered, file="Y2H-SCORES_files/y2h_scores_results_summary_filtered.csv",quote = FALSE, row.names = FALSE)


#check results
saint_results_summary<- read.csv(file="SAINTexpress_files/saint_results_summary.csv", row.names = NULL)
unique(saint_results_summary$Bait)
hist(saint_results_summary$SaintScore, breaks = length(unique(saint_results_summary$SaintScore)))
sum(saint_results_summary$SaintScore>0.9)
unique(saint_results_summary[saint_results_summary$SaintScore>0.9, "Prey"])
sum(saint_results_summary$SaintScore>0.85)
unique(saint_results_summary[saint_results_summary$SaintScore>0.88, "Prey"])
saint_results_summary_filtered<- saint_results_summary[saint_results_summary$SaintScore>0.88, ]
#take top 30 comparisons in "Bc1_AB"            "Infected_AB"       "Uninfected_BC1_AB" 

y2h_scores_results_summary<- read.csv(file="Y2H-SCORES_files/y2h_scores_results_summary.csv",row.names = NULL)
y2h_scores_results_summary_filtered2<- y2h_scores_results_summary[y2h_scores_results_summary$bait %in% c("Bc1_AB","Infected_AB","Uninfected_BC1_AB"),]
y2h_scores_results_summary_filtered2$rank<- rank(-y2h_scores_results_summary_filtered2$Borda_scores, ties.method = "random")
y2h_scores_results_summary_uninfected_AB<- y2h_scores_results_summary[y2h_scores_results_summary$bait=="Uninfected_AB",]

y2h_scores_results_summary_filtered3<- y2h_scores_results_summary_filtered2[y2h_scores_results_summary_filtered2$rank<=30,]
y2h_scores_results_summary_filtered3$prey %in% y2h_scores_results_summary_uninfected_AB[y2h_scores_results_summary_uninfected_AB$Borda_scores> min(y2h_scores_results_summary_filtered3$Borda_scores),"Prey"]
unique(y2h_scores_results_summary_filtered3$prey)

write.csv(y2h_scores_results_summary_filtered3, "y2h_scores_results_summary_filteredV2.csv", row.names = F)
write.csv(saint_results_summary_filtered, "saint_results_summary_filteredV2.csv", row.names = F)

#compare 
sum(unique(saint_results_summary_filtered$Prey) %in% unique(y2h_scores_results_summary_filtered3$prey))
#21 unique matches
#format and save file
BC1_interactors<- data.frame(Arabidopsis.homologue=c(unique(saint_results_summary_filtered$Gene.names...ordered.locus..),
                                                     unique(y2h_scores_results_summary_filtered3$Gene.names...ordered.locus..)),
                             source= c(rep("SAINTexpress", length(unique(saint_results_summary_filtered$Gene.names...ordered.locus..))),
                                       rep("Y2H-SCORES", length(unique(y2h_scores_results_summary_filtered3$Gene.names...ordered.locus..)))))

BC1_interactors<- BC1_interactors[!is.na(BC1_interactors$Arabidopsis.homologue),]
BC1_interactors$Arabidopsis.homologue<- substr(BC1_interactors$Arabidopsis.homologue, 1, 9)
BC1_interactors$source<- ifelse(duplicated(BC1_interactors$Arabidopsis.homologue), "SAINTexpress;Y2H-SCORES",BC1_interactors$source)
BC1_interactors<- BC1_interactors[-which(BC1_interactors[BC1_interactors$source=="SAINTexpress", "Arabidopsis.homologue"] %in% 
                                           BC1_interactors[BC1_interactors$source=="SAINTexpress;Y2H-SCORES", "Arabidopsis.homologue"]),]

write.csv(BC1_interactors, "filtered_interactors.csv", row.names = F)



