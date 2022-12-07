library (dplyr)
library (ggplot2)
library (cowplot)
library (scico)
library (patchwork)

### Sequnce overview data ###

sample_pairs <- read.csv ("Results/Nextclade_hhid.csv", colClasses=c("Household" = "character"))
snv_pairs <- read.csv ("Results/SNV_from_households_with_suffecient_sequencing.csv", colClasses=c("Household" = "character"))
snv_pairs_wide <- read.csv ("Results/snv_pairs_wide.csv", colClasses=c("Household" = "character"))

#  of individuals per household by clade
sample_pairs_count <- count (sample_pairs, Household, VOC )
sample_pairs_count2 <- count (sample_pairs_count,n, VOC )
sample_pairs_count2$VOC <- factor(sample_pairs_count2$VOC , levels=c("Non-VOC", "Alpha","Gamma", "Delta", "Omicron"))
scico_colors <- scico(5, alpha = NULL, begin = 0, end = .85, direction = -1, palette = "batlow")
indiv_house_plot <- ggplot (sample_pairs_count2, aes (n, nn, fill= VOC))+ geom_bar(position="stack", stat="identity", alpha=.7)+ xlab ("Individuals per household")+ ylab ("Households")+theme_cowplot(12)+scale_fill_manual(values=scico_colors)+ labs(fill='Variant')



# Num of SNV per individual
snv_sample <- count (snv_pairs, sample)
snv_sample2<- count (snv_sample, n)
snv_sample2[nrow(snv_sample2) + 1,] = c("0","56")
snv_sample2$nn <- as.numeric (snv_sample2$nn)
snv_sample2$n <- as.numeric (snv_sample2$n)

snv_per_sample_plot <- ggplot (snv_sample2, aes (n, nn))+ geom_col () + theme_cowplot(12)+ xlab ("iSNV per sample")+ylab ("Number of samples")+ scale_x_continuous(breaks = seq(0, 10, by = 2))
#+ theme (axis.line = element_line(colour = 'black', size = .25), axis.ticks = element_line(colour = "black", size = .25))

## Frequency of SNV
snv_per_freq <- ggplot (snv_pairs, aes (avg_freq))+geom_histogram(binwidth =0.05, boundary=0)+ theme_cowplot(12)+ xlab ("iSNV frequency")+ ylab ("Number of iSNV")

# frequency of samples
snv_trans_freq_pairs <-ggplot (snv_pairs_wide, aes (Donor_freq_fix, Recepient_freq_fix ))+ geom_point (alpha =.6, size =1.3)+ theme_cowplot(12)+ xlab ("Frequency in Donor")+ ylab ("Frequency in Recipient")


## plot Figures
plot_grid ( snv_per_sample_plot,snv_per_freq,
  labels = c('A', 'B')
)

ggasave ("Results/Plots/Figure_2.pdf")

plot_grid(
  indiv_house_plot,snv_trans_freq_pairs,
  labels = c('A', 'B', )
)

ggasave ("Results/Plots/Figure_3.pdf")


## Supplemental Figure 1 

# Genome wide coverage plot

 Coverage <- read.csv("Coverage.all", stringsAsFactors = FALSE)
 
 ###### Functions #####
 #coverage
 
 slide <- function(cov.df, setup.df)
{
  cov = rep(NA, nrow(setup.df))
  for(i in 1:nrow(setup.df))
  {
    s = setup.df$starts[i]
    e = setup.df$ends[i]
    subset(cov.df, pos >= s & pos < e, select = c(cov)) -> position
    mean(position$cov) -> cov[i]
  }
  out <- data.frame(mean =cov, pos = setup.df$pos)
  out$ID = unique(cov.df$ID)
  return(out)
}

cov_plot <- function(cov.df, title)
{
  cov.df  %>% summarize(first = min(pos), last = max(pos)) %>% plyr::adply(1,function(x) data.frame(starts = seq(x$first,x$last,by=400))) %>% mutate(ends = ifelse(starts + 400 < last, starts + 400, last)) -> setup
  setup %>% select(starts, ends) -> setup_means
  setup$pos <- apply(setup_means, 1, function(x) mean(x))
  plyr::ddply(cov.df, ~ID, slide, setup) -> cov.slid.df
  
  
  cov.plot <- ggplot(cov.slid.df, mapping = aes(x = as.factor(pos), y = mean)) + geom_boxplot(fill="white")
  cov.plot <- cov.plot + ggtitle(title) + ylab("Read depth")  + xlab(" Genome Position")
  cov.plot <- cov.plot + theme(axis.title.y = element_text(vjust=1.2))
  cov.plot <- cov.plot + theme(legend.position = "none") + theme_classic()
  cov.plot <- cov.plot + theme(text = element_text(size = 15), axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 10))
  return(cov.plot)
}

 # ======================== Plot ===================


cov_plot(Coverage, title = "") -> coverage.plot.all


# replicate iSNV frequency
snv <- read.table ("Processed_data/all_variants_filtered", header=T)

snv_rep_plot <- ggplot (snv, aes (ALT_FREQ_1, ALT_FREQ_2))+ geom_point()+ xlab ("iSNV frequency in replicate 1")+ ylab ("iSNV frequency in replicate 2")+theme_bw()
plot_insert <- ggplot (snv, aes (ALT_FREQ_1, ALT_FREQ_2))+ geom_point(size=0.3)+ xlab ("")+ ylab ("")+theme_bw()+xlim(0,.1)+ylim (0,0.1)

snv_rep_plot + inset_element(plot_insert, 0.01, .6, 0.4, .99)
