library (dplyr)
library (ggplot2)

avg_cov <- read.csv ("Processed_data/avg_coverage_pass_fail.csv")
avg_cov_pass <- filter (avg_cov,Pass_coverage == "Pass" )

snv_merged <- read.csv ("Processed_data/all_variants_filtered_merged")
snv_merged <- filter (snv_merged, sample %in% avg_cov_pass$sample)

snv_merged<- unite (snv_merged, sample_mutation, c (sample, mutation), sep= "-", remove =F)
snv_merged<- mutate (snv_merged, mutation_type=ifelse(is.na(REF_AA), "non_coding", ifelse(REF_AA==ALT_AA, "Syn", "Non")))

#mask 
problem_sites <- read.table ("ncov_references/problematic_sites_v7.txt", header=T)
snv_merged<- filter (snv_merged, !POS %in% problem_sites$POS)
write.csv (snv_merged, "Results/snv_filtered_merged.csv", row.names=F, quote=F)


## filter by high quality sequenced families
sample_pairs <- read.csv ("Processed_data/Nextclade_hhid.csv", colClasses = c(household = "character"))


snv_pairs_merged <- filter (snv_merged, sample %in% sample_pairs$sample)
snv_pairs_merged <- left_join (snv_pairs_merged, sample_pairs, by="sample")

snv_freq_merged <- ggplot (snv_pairs_merged, aes (ALT_FREQ))+geom_histogram(binwidth =0.05, boundary=0)+ theme_cowplot(12)+ xlab ("iSNV frequency")+ ylab ("# of iSNV")
ggsave ("Results/Plots/SNV_feequency_merged_reads.pdf")


snv_sample_merged <- count (snv_pairs_merged, sample)
snv_sample2_merged<- count (snv_sample_merged, n)
snv_sample2_merged[nrow(snv_sample2_merged) + 1,] = c("0","39")
snv_sample2_merged$nn <- as.numeric (snv_sample2_merged$nn)
snv_sample2_merged$n <- as.numeric (snv_sample2_merged$n)


snv_per_sample_plot_merged <- ggplot (snv_sample2_merged, aes (n, nn))+ geom_col () + theme_cowplot(12)+ xlab ("# of iSNV per sample")+ylab ("# of samples")+ scale_x_continuous(breaks = seq(0, 60, by = 5))
+ theme (axis.line = element_line(colour = 'black', size = .25), axis.ticks = element_line(colour = "black", size = .25))

plot_grid(snv_per_sample_plot_merged, snv_freq_merged, labels= c('A', 'B'), ncol = 2, label_size = 12, rel_widths=c(1,1))

## get positions for consensus allele
snv_positions_merged <- distinct (snv_pairs_merged, POS)
write.csv (snv_positions_merged,"Results/snv_positions_merged.csv", row.names=F, quote=F)

## add in consensus allele - run python script first

consensus_allele_merged <- read.table ("Results/consensus_allele_forsnv_sites_merged.txt", header=T)
consensus_allele_merged <- right_join ( consensus_allele_merged, sample_pairs, by = "sample")
consensus_allele_merged <-  consensus_allele_merged %>% rename (consensus_allele =Allele)
consensus_allele_merged <- unite (consensus_allele_merged, household_position, c (Household, POS), sep= "_", remove=F )
consensus_allele_merged <- filter (consensus_allele_merged, household_position%in% snv_pairs_merged$household_position)

snv_pairs_merged<- unite (snv_pairs_merged, household_position, c (Household, POS), sep= "_", remove=F )
snv_pairs_sample_consensus_merged <- left_join (consensus_allele_merged, snv_pairs_merged, by = c ("sample", "household_position", "POS", "Household", "study_id") )


trans_pairs <- read.csv ("Processed_data/Transmission_pairs.csv", colClasses = c(household = "character"))
trans_pairs <- rename (trans_pairs, Household = household)


snv_pairs_sample_consensus_merged <- left_join (snv_pairs_sample_consensus_merged,trans_pairs, by = c("sample", "study_id", "Household") )
snv_pairs_sample_consensus_merged <- select (snv_pairs_sample_consensus_merged, sample, consensus_allele, household_position, POS, study_id, Household, REF, ALT, mutation, mutation_type, ALT_FREQ, pair_id, Transmission_indiv, VOC.x)


# split into donor and reciepients
Donor_individual_merged <- filter (snv_pairs_sample_consensus_merged, Transmission_indiv == "A")
Recepient_individual_merged<- filter (snv_pairs_sample_consensus_merged, Transmission_indiv == "B")
Donor_individual_merged <- rename (Donor_individual_merged,Donor_sample=sample, Donor_consensus_allele = consensus_allele, Donor_study_id=study_id, Donor_REF=REF, Donor_ALT =ALT, Donor_mutation=mutation, Donor_ALT_FREQ =ALT_FREQ , Donor_mutation_type=mutation_type)

Recepient_individual_merged <- rename (Recepient_individual_merged, Recepient_sample=sample, Recepient_consensus_allele = consensus_allele, Recepient_study_id=study_id,  Recepient_REF=REF, Recepient_ALT =ALT, Recepient_mutation=mutation, Recepient_ALT_FREQ =ALT_FREQ , Recepient_mutation_type=mutation_type)

## make wide data frame with donor and recepient on one line
snv_pairs_wide_merged <- left_join (Donor_individual_merged,Recepient_individual_merged, by = c ("household_position", "POS", "pair_id", "Household") )
snv_pairs_wide_merged<- filter (snv_pairs_wide_merged, !is.na( Donor_REF) | !is.na (Recepient_REF))


snv_pairs_wide_merged  <- mutate (snv_pairs_wide_merged, Donor_freq_fix = ifelse (!is.na (Donor_ALT_FREQ), Donor_ALT_FREQ, ifelse (Donor_consensus_allele == Recepient_ALT, "1", "0")))
snv_pairs_wide_merged  <- mutate (snv_pairs_wide_merged, Recepient_freq_fix = ifelse (!is.na (Recepient_ALT_FREQ), Recepient_ALT_FREQ, ifelse (Recepient_consensus_allele == Donor_ALT, "1", "0")))

snv_pairs_wide_merged$Donor_freq_fix <- as.numeric (snv_pairs_wide_merged$Donor_freq_fix)
snv_pairs_wide_merged$Recepient_freq_fix <- as.numeric (snv_pairs_wide_merged$Recepient_freq_fix)
snv_trans_freq_pairs_merged <-ggplot (snv_pairs_wide_merged, aes (Donor_freq_fix, Recepient_freq_fix ))+ geom_point ()+ theme_bw()+ xlab ("Frequency in Donor")+ ylab ("Frequency in Recipient")

Donor_snv_merged <- filter (snv_pairs_wide_merged, !is.na (Donor_mutation ))


## prep input file for bottleneck analysis 

Donor_snv_bottleneck_merged <- rename (Donor_snv_merged, pos=POS, freq1 = Donor_freq_fix, freq2 = Recepient_freq_fix)
Donor_snv_bottleneck_merged <- mutate (Donor_snv_bottleneck_merged)

Donor_snv_bottleneck_merged$freq2 <- as.numeric (Donor_snv_bottleneck_merged$freq2)
Donor_snv_bottleneck2_merged <- Donor_snv_bottleneck_merged %>% mutate (freq1 = 1-freq1, freq2 = 1-freq2, Donor_ALT= Donor_REF)

Donor_snv_bottleneck_merged <- bind_rows (Donor_snv_bottleneck_merged, Donor_snv_bottleneck2_merged)

### for betabinomial aproximate

Donor_snv_bottleneck_found_merged <- filter (Donor_snv_bottleneck_merged, freq2>0)
Donor_snv_bottleneck_found_merged <- unite (Donor_snv_bottleneck_found_merged,  pairID_mutation, c (pair_id, Donor_mutation), sep= "_", remove=F)
Donor_snv_bottleneck_found_merged <- distinct (Donor_snv_bottleneck_found_merged,pairID_mutation,.keep_all=T )
Donor_snv_bottleneck_found_merged <- rename ( Donor_snv_bottleneck_found_merged,VOC = VOC.x.x )


## by transmission pair
Donor_snv_bottleneck_found_merged$pair_id <- as.factor (Donor_snv_bottleneck_found_merged$pair_id)
for(i in levels(Donor_snv_bottleneck_found_merged$pair_id)){
	trans_freq_pair <- Donor_snv_bottleneck_found_merged %>% filter ( pair_id == i)%>% select (freq1, freq2 )
	filename <- paste ("Results/Betabinomial/snv_betabinomial_merged",i, sep= "_")
	write.table ( trans_freq_pair, filename, sep= "\t", col.names=FALSE, row.names=FALSE)
}