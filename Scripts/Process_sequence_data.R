library (dplyr)
library (ggplot2)
library (stringr)
library (tidyr)


snv <- read.table ("Processed_data/all_variants_filtered", header=T)

### average coverage- get samples with 500x coverage in both replicates
avg_cov <- read.table ("Processed_data/AvgCoverage.all", header=T)
avg_cov <- separate (avg_cov, sample, c ("sample", "replicate"), "_") # split sample column
avg_cov <- spread (avg_cov, replicate, mean)
avg_cov<- rename (avg_cov, rep_2= 3, rep_1= 2)
avg_cov <- mutate (avg_cov, Pass_coverage = ifelse (rep_1 >= 500 & rep_2 >=500, "Pass", "Fail"))
write.csv (avg_cov, "Processed_data/avg_coverage_pass_fail.csv", quote=F, row.names=F)


## Filter SNV by coverage 
avg_cov_pass <- filter (avg_cov,Pass_coverage == "Pass" )
snv <- filter (snv, sample %in% avg_cov_pass$sample)


### get SNV type
snv <- mutate (snv, avg_freq = (ALT_FREQ_1+ALT_FREQ_2)/2)
snv<- unite (snv, sample_mutation, c (sample, mutation), sep= "-", remove =F)
snv<- mutate (snv , mutation_type=ifelse(is.na(REF_AA), "non_coding", ifelse(REF_AA==ALT_AA, "Syn", "Non")))

#mask 
problem_sites <- read.table ("ncov_references/problematic_sites_v7.txt", header=T)
snv<- filter (snv, !POS %in% problem_sites$POS)

write.csv (snv, "Results/snv_filtered.csv", row.names=F, quote=F)

### 

## add clades in 
sample_pairs <- read.csv ("Processed_data/Nextclade_hhid.csv", colClasses = c(household = "character"))

#Keep SNV from appropriate individuals
snv_pairs <- filter (snv, sample %in% sample_pairs$sample)
snv_pairs <- left_join (snv_pairs, sample_pairs, by="sample")
write.csv (snv_pairs, "Results/SNV_from_households_with_suffecient_sequencing.csv", row.names=F, quote=F)


### transmission pair set up 

trans_pairs <- read.csv ("Processed_data/Transmission_pairs.csv", colClasses = c(household = "character"))
trans_pairs <- rename (trans_pairs, Household = household)

snv_pairs<- unite (snv_pairs, household_position, c (Household, POS), sep= "_", remove=F )

## get positions for consensus allele
snv_positions_merged <- distinct (snv_pairs_merged, POS)
write.csv (snv_positions_merged,"Results/snv_positions.csv", row.names=F, quote=F)


## add consensus alleles - Run Python script first

consensus_allele <- read.table ("Results/consensus_allele_forsnv_sites.txt", header=T)
consensus_allele <- right_join ( consensus_allele, sample_pairs, by = "sample")
consensus_allele <-  consensus_allele %>% rename (consensus_allele =Allele)
consensus_allele <- unite (consensus_allele, household_position, c (Household, POS), sep= "_", remove=F )
consensus_allele <- filter (consensus_allele, household_position%in% snv_pairs$household_position)


snv_pairs_sample_consensus <- left_join (consensus_allele, snv_pairs, by = c ("sample", "household_position", "POS", "Household", "study_id") )


snv_pairs_sample_consensus <- left_join (snv_pairs_sample_consensus,trans_pairs, by = c("sample", "study_id", "Household") )
snv_pairs_sample_consensus <- select (snv_pairs_sample_consensus, sample, consensus_allele, household_position, POS, study_id, Household, REF, ALT, mutation, mutation_type, avg_freq, pair_id, Transmission_indiv, VOC.x)


# split into donor and reciepients
Donor_individual <- filter (snv_pairs_sample_consensus, Transmission_indiv == "A")
Recepient_individual<- filter (snv_pairs_sample_consensus, Transmission_indiv == "B")
Donor_individual <- rename (Donor_individual,Donor_sample=sample, Donor_consensus_allele = consensus_allele, Donor_study_id=study_id, Donor_REF=REF, Donor_ALT =ALT, Donor_mutation=mutation, Donor_avg_freq =avg_freq , Donor_mutation_type=mutation_type)

Recepient_individual <- rename (Recepient_individual, Recepient_sample=sample, Recepient_consensus_allele = consensus_allele, Recepient_study_id=study_id,  Recepient_REF=REF, Recepient_ALT =ALT, Recepient_mutation=mutation, Recepient_avg_freq =avg_freq , Recepient_mutation_type=mutation_type)

## make wide data frame with donor and recepient on one line
snv_pairs_wide <- left_join (Donor_individual,Recepient_individual, by = c ("household_position", "POS", "pair_id", "Household") )
snv_pairs_wide <- filter (snv_pairs_wide, !is.na( Donor_REF) | !is.na (Recepient_REF))


snv_pairs_wide  <- mutate (snv_pairs_wide, Donor_freq_fix = ifelse (!is.na (Donor_avg_freq), Donor_avg_freq, ifelse (Donor_consensus_allele == Recepient_ALT, "1", "0")))
snv_pairs_wide  <- mutate (snv_pairs_wide, Recepient_freq_fix = ifelse (!is.na (Recepient_avg_freq), Recepient_avg_freq, ifelse (Recepient_consensus_allele == Donor_ALT, "1", "0")))


snv_pairs_wide$Donor_freq_fix <- as.numeric (snv_pairs_wide$Donor_freq_fix)
snv_pairs_wide$Recepient_freq_fix <- as.numeric (snv_pairs_wide$Recepient_freq_fix)
write.csv (snv_pairs_wide, "Results/snv_pairs_wide.csv", quote=F, row.names=F)


Donor_snv <- filter (snv_pairs_wide, !is.na (Donor_mutation ))

## prep input file for bottleneck analysis 

Donor_snv_bottleneck <- rename (Donor_snv, pos=POS, freq1 = Donor_freq_fix, freq2 = Recepient_freq_fix)

Donor_snv_bottleneck2 <- Donor_snv_bottleneck 
Donor_snv_bottleneck$freq2 <- as.numeric (Donor_snv_bottleneck$freq2)
Donor_snv_bottleneck2 <- Donor_snv_bottleneck %>% mutate (freq1 = 1-freq1, freq2 = 1-freq2, Donor_ALT= Donor_REF)

Donor_snv_bottleneck <- bind_rows (Donor_snv_bottleneck, Donor_snv_bottleneck2)

# only keep snv that are foun in the donor
Donor_snv_bottleneck_found <- filter (Donor_snv_bottleneck, freq2>0)
Donor_snv_bottleneck_found <- unite (Donor_snv_bottleneck_found,  pairID_mutation, c (pair_id, Donor_mutation), sep= "_", remove=F)
Donor_snv_bottleneck_found <- distinct (Donor_snv_bottleneck_found,pairID_mutation,.keep_all=T )
Donor_snv_bottleneck_found <- rename ( Donor_snv_bottleneck_found,VOC = VOC.x.x )

# write file for each individual
Donor_snv_bottleneck_found$pair_id <- as.factor (Donor_snv_bottleneck_found$pair_id)
for(i in levels(Donor_snv_bottleneck_found$pair_id)){
	trans_freq_pair <- Donor_snv_bottleneck_found %>% filter ( pair_id == i)%>% select (freq1, freq2 )
	filename <- paste ("Results/BetaBinomial_bottleneck/snv_betabinmial",i, sep= "_")
	write.table ( trans_freq_pair, filename, sep= "\t", col.names=FALSE, row.names=FALSE)
}


