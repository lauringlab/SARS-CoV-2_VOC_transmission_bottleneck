## first set of figures - timing of infection and sample collection 
sym_onset <- read.csv ("Processed_data/sym_onset.csv", colClasses=c("Household" = "character"))

sym_onset$VOC <- factor(sym_onset$VOC , levels=c("Non-VOC", "Alpha","Gamma", "Delta", "Lambda", "Omicron"))


# time from index symptom onset to collection time
sym_onset_index <- filter (sym_onset, index ==1)
sym_onset_index.summary <- sym_onset_index %>% group_by (VOC) %>% summarize (time_from_index_sym_to_collection = median (time_from_index_sym_to_collection))

sym_onset_index$VOC <- factor(sym_onset_index$VOC , levels=c("Non-VOC", "Alpha", "Delta",  "Omicron"))

index_plot <- ggplot (sym_onset_index, aes (VOC, time_from_index_sym_to_collection, fill = study, color=study))+geom_crossbar(data=sym_onset_index.summary, aes(ymin = time_from_index_sym_to_collection, ymax = time_from_index_sym_to_collection, fill = NULL, color=NULL),size=.3, width = .5)+  geom_dotplot (stackdir="center",binaxis= "y", position=position_dodge(width = 0.25))+theme_cowplot (12) + ylab ("Days till collection")+ xlab ("Clade")+ylim (0, 11)+ scale_fill_manual (values=c ("#EA4025","#479D8B"), guide="none")+scale_color_manual (values=c ("#EA4025","#479D8B"), guide="none")

# ANOVA
index_symp_onset_anova <- aov ( time_from_index_sym_to_collection~ VOC + study ,data= sym_onset_index)
summary (index_symp_onset_anova)

## time from index symptom onset to contact symptom onset
sym_onset_contact <- filter (sym_onset,index ==0 )
sym_onset_contact$VOC <- factor(sym_onset_contact$VOC , levels=c("Non-VOC", "Alpha","Gamma", "Delta", "Lambda", "Omicron"))

sym_onset_contact.summary <- sym_onset_contact %>% group_by (VOC) %>% summarize (time_index_sym_contact_sypm = median (time_index_sym_contact_sypm, na.rm = TRUE))

sym_onset_contact_plot <- ggplot (sym_onset_contact, aes (VOC, time_index_sym_contact_sypm))+ geom_dotplot (stackdir="center",binaxis= "y")+theme_cowplot (12) + ylab ("Serial interval (days)")+ xlab ("Clade")+ geom_crossbar(data=sym_onset_contact.summary, aes(ymin = time_index_sym_contact_sypm, ymax = time_index_sym_contact_sypm),size=.3, width = .5)+ ylim (0,11)

# ANOVA

symp_onset_contact_anova <- aov ( time_index_sym_contact_sypm~ VOC + study ,data= sym_onset_contact)
summary (symp_onset_contact_anova )

## CT values for index cases
sequenced_indexes <- select (sym_onset_index, study_id,  Household, index_symptom_onset, VOC )
ct <- read.csv ("Processed_data/ct_values.csv")
ct <- filter (ct, study_id %in% sequenced_indexes$study_id)
ct_index <-inner_join (ct,sequenced_indexes, by = "study_id")
ct_index <- mutate (ct_index, sequenced = ifelse (ACCN %in% sym_onset$sample, "TRUE", "FALSE"))
ct_index$Collection_Date <- as.Date (ct_index$Collection_Date, "%m/%d/%Y")
ct_index$index_symptom_onset<- as.Date (ct_index$index_symptom_onset, "%m/%d/%y")
ct_index <- mutate (ct_index, day_since_sym_onset=Collection_Date-index_symptom_onset)



## Plot Figure 1
ct_index$VOC <- factor(ct_index$VOC , levels=c("Non-VOC", "Alpha","Gamma", "Delta", "Lambda", "Omicron"))

ct_plot <-  ggplot (ct_index, aes (day_since_sym_onset, N_gene_ct, group=study_id )) + geom_line (alpha =.7)+geom_point(shape = 21,aes (fill = sequenced ))+  facet_wrap(vars(VOC))+scale_fill_manual(values = c('TRUE'="black", 'FALSE'="white"))+theme_bw(12) +theme(panel.grid.major= element_blank(), panel.grid.minor= element_blank(), strip.background =element_rect(fill="grey85"), legend.position="none")+ xlab ("Days since symptom onset")+ ylab ("Cycle threshold")+xlim (0,20)+scale_y_reverse()

top <- plot_grid(sym_onset_contact_plot, index_plot2, labels= c('A', 'B'), ncol = 2, label_size = 12, rel_widths=c(1.2,1))


plot_grid(top, ct_plot, labels = c('', 'C'), label_size = 12, ncol = 1, rel_heights = c(1, 1.5))
ggsave ("Results/Plots/Figure1_infection_timeline.pdf")

## Supplementary Figure - Infection timeline

sym_onset <- arrange (sym_onset,  index)
sym_onset$index <- as.factor (sym_onset$index )
cols <- c("0"= "black", "1"= "blue")

sym_onset_Non <- filter (sym_onset, VOC== "Non-VOC")
sym_onset_Non_plot <- ggplot (sym_onset_Non)+ geom_point (aes (time_index_sym_contact_sypm, new_id, color= index), shape = 2)+  geom_point (aes( time_from_index_sym_to_collection, new_id, fill= index, color=index), shape =25)+ theme_bw () + scale_color_manual (values =cols)+ scale_fill_manual (values =cols)+ xlab ("Days since index symptom onset") + ylab ("Individual")+ facet_wrap(vars(New_household), scales="free_y", ncol=4)+xlim (0,11.5)

sym_onset_Alpha <- filter (sym_onset, VOC== "Alpha")
sym_onset_alpha_plot <- ggplot (sym_onset_Alpha)+ geom_point (aes (time_index_sym_contact_sypm, new_id, color= index), shape = 2)+  geom_point (aes( time_from_index_sym_to_collection, new_id, fill= index, color=index), shape =25)+ theme_bw () + scale_color_manual (values =cols)+ scale_fill_manual (values =cols)+ xlab ("Days since index symptom onset") + ylab ("Individual")+ facet_wrap(vars(New_household), scales="free_y", ncol=4)+xlim (0,11.5)

sym_onset_Delta <- filter (sym_onset, VOC== "Delta")
sym_onset_delta_plot <- ggplot (sym_onset_Delta)+ geom_point (aes (time_index_sym_contact_sypm, new_id, color= index), shape = 2)+  geom_point (aes( time_from_index_sym_to_collection, new_id, fill= index, color=index), shape =25)+ theme_bw () + scale_color_manual (values =cols)+ scale_fill_manual (values =cols)+ xlab ("Days since index symptom onset") + ylab ("Individual")+ facet_wrap(vars(New_household), scales="free_y", ncol=4)+xlim (0,11.5)

sym_onset_Gamma <- filter (sym_onset, VOC== "Gamma")
sym_onset_gamma_plot <- ggplot (sym_onset_Gamma)+ geom_point (aes (time_index_sym_contact_sypm, new_id, color= index), shape = 2)+  geom_point (aes( time_from_index_sym_to_collection, new_id, fill= index, color=index), shape =25)+ theme_bw () + scale_color_manual (values =cols)+ scale_fill_manual (values =cols)+ xlab ("Days since index symptom onset") + ylab ("Individual")+ facet_wrap(vars(New_household), scales="free_y", ncol=4)+xlim (0,11.5)

sym_onset_Omicron <- filter (sym_onset, VOC== "Omicron")
sym_onset_omicron_plot <- ggplot (sym_onset_Omicron)+ geom_point (aes (time_index_sym_contact_sypm, new_id, color= index), shape = 2)+  geom_point (aes( time_from_index_sym_to_collection, new_id, fill= index, color=index), shape =25)+ theme_bw () + scale_color_manual (values =cols)+ scale_fill_manual (values =cols)+ xlab ("Days since index symptom onset") + ylab ("Individual")+ facet_wrap(vars(New_household), scales="free_y", ncol=4)+xlim (0,11.5)

set_panel_size <- function(p=NULL, g=ggplotGrob(p), file=NULL, 
                           margin = unit(1,"mm"),
                           width=unit(3, "cm"), 
                           height=unit(2, "cm")){

  panels <- grep("panel", g$layout$name)
  panel_index_w<- unique(g$layout$l[panels])
  panel_index_h<- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  
  if(getRversion() < "3.3.0"){

   # the following conversion is necessary
   # because there is no `[<-`.unit method
   # so promoting to unit.list allows standard list indexing
   g$widths <- grid:::unit.list(g$widths)
   g$heights <- grid:::unit.list(g$heights)

   g$widths[panel_index_w] <-  rep(list(width),  nw)
   g$heights[panel_index_h] <- rep(list(height), nh)

} else {

   g$widths[panel_index_w] <-  rep(width,  nw)
   g$heights[panel_index_h] <- rep(height, nh)

}

  if(!is.null(file))
    ggsave(file, g, 
           width = grid::convertWidth(sum(g$widths) + margin, 
                                unitTo = "in", valueOnly = TRUE),
           height = grid::convertHeight(sum(g$heights) + margin,  
                                  unitTo = "in", valueOnly = TRUE))

  g
}

g1 <- set_panel_size(sym_onset_Non_plot)
g2 <- set_panel_size(sym_onset_alpha_plot)
g3<- set_panel_size (sym_onset_delta_plot)
g4<- set_panel_size (sym_onset_gamma_plot)
g5<- set_panel_size (sym_onset_omicron_plot)

g1plot <- gridExtra::grid.arrange(g1)
g2plot <- gridExtra::grid.arrange(g2)
g3plot <- gridExtra::grid.arrange(g3)
g4plot <- gridExtra::grid.arrange(g4)
g5plot <- gridExtra::grid.arrange(g5)

ggsave ("Results/Plots/Supl_infection_timeline_nonvoc.pdf", g1plot)
ggsave ("Results/Plots/Supl_infection_timeline_alpha.pdf", g2plot)
ggsave ("Results/Plots/Supl_infection_timeline_delta.pdf", g3plot)
ggsave ("Results/Plots/Supl_infection_timeline_gamma.pdf", g4plot)
ggsave ("Results/Plots/Supl_infection_timeline_omicron.pdf", g5plot)
