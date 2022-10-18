non_voc <- read.table ("Non-VOC/Non-VOC_likelihood", header=T)
non_voc <- non_voc %>% group_by(bottleneck_size)%>% mutate (adjusted_LL= (max(Num_sites)/Num_sites)*Log_Likelihood)
LL_sum <- summarise (non_voc, sum (adjusted_LL))
LL_sum <- rename (LL_sum, adjusted_LL = 2)

bottleneck_values_vector <- c()
for ( i in 1:200)
{  if(i%%1 == 0) {bottleneck_values_vector <- c(bottleneck_values_vector,i)}
  
  
}


Max_LL <- max(LL_sum$adjusted_LL) # Maximum value of log likelihood
Max_LL_bottleneck_index <- which(LL_sum$adjusted_LL == max(LL_sum$adjusted_LL) ) # bottleneck size at which max likelihood occurs
Max_LL_bottleneck <- bottleneck_values_vector[Max_LL_bottleneck_index] 
likelihood_ratio <- qchisq(.95, df=1) # necessary ratio of likelihoods set by confidence level
ci_tibble <- filter(LL_sum, 2*(Max_LL - adjusted_LL) <= likelihood_ratio ) 
lower_CI_bottleneck <- min(ci_tibble$bottleneck_size) #-1 # lower bound of confidence interval
upper_CI_bottleneck <- max(ci_tibble$bottleneck_size) #+1# upper bound of confidence interval






alpha <- read.table ("Alpha/Alpha_likelihood", header=T)
alpha <- alpha %>% group_by(bottleneck_size)%>% mutate (adjusted_LL= (max(Num_sites)/Num_sites)*Log_Likelihood)
LL_sum <- summarise (alpha, sum (adjusted_LL))
LL_sum <- rename (LL_sum, adjusted_LL = 2)

bottleneck_values_vector <- c()
for ( i in 1:200)
{  if(i%%1 == 0) {bottleneck_values_vector <- c(bottleneck_values_vector,i)}
  
  
}


Max_LL <- max(LL_sum$adjusted_LL) # Maximum value of log likelihood
Max_LL_bottleneck_index <- which(LL_sum$adjusted_LL == max(LL_sum$adjusted_LL) ) # bottleneck size at which max likelihood occurs
Max_LL_bottleneck <- bottleneck_values_vector[Max_LL_bottleneck_index] 
likelihood_ratio <- qchisq(.95, df=1) # necessary ratio of likelihoods set by confidence level
ci_tibble <- filter(LL_sum, 2*(Max_LL - adjusted_LL) <= likelihood_ratio ) 
lower_CI_bottleneck <- min(ci_tibble$bottleneck_size) #-1 # lower bound of confidence interval
upper_CI_bottleneck <- max(ci_tibble$bottleneck_size) #+1# upper bound of confidence interval


write.table 

delta <- read.table ("Delta/Delta_likelihood", header=T)
delta <- delta %>% group_by(bottleneck_size)%>% mutate (adjusted_LL= (max(Num_sites)/Num_sites)*Log_Likelihood)
LL_sum <- summarise (delta, sum (adjusted_LL))
LL_sum <- rename (LL_sum, adjusted_LL = 2)

bottleneck_values_vector <- c()
for ( i in 1:200)
{  if(i%%1 == 0) {bottleneck_values_vector <- c(bottleneck_values_vector,i)}
  
  
}


Max_LL <- max(LL_sum$adjusted_LL) # Maximum value of log likelihood
Max_LL_bottleneck_index <- which(LL_sum$adjusted_LL == max(LL_sum$adjusted_LL) ) # bottleneck size at which max likelihood occurs
Max_LL_bottleneck <- bottleneck_values_vector[Max_LL_bottleneck_index] 
likelihood_ratio <- qchisq(.95, df=1) # necessary ratio of likelihoods set by confidence level
ci_tibble <- filter(LL_sum, 2*(Max_LL - adjusted_LL) <= likelihood_ratio ) 
lower_CI_bottleneck <- min(ci_tibble$bottleneck_size) #-1 # lower bound of confidence interval
upper_CI_bottleneck <- max(ci_tibble$bottleneck_size) #+1# upper bound of confidence interval







omicron<- read.table ("Omicron/Omicron_likelihood", header=T)
omicron <- omicron %>% group_by(bottleneck_size)%>% mutate (adjusted_LL= (max(Num_sites)/Num_sites)*Log_Likelihood)
LL_sum <- summarise (omicron, sum (adjusted_LL))
LL_sum <- rename (LL_sum, adjusted_LL = 2)

bottleneck_values_vector <- c()
for ( i in 1:200)
{  if(i%%1 == 0) {bottleneck_values_vector <- c(bottleneck_values_vector,i)}
  
  
}


Max_LL <- max(LL_sum$adjusted_LL) # Maximum value of log likelihood
Max_LL_bottleneck_index <- which(LL_sum$adjusted_LL == max(LL_sum$adjusted_LL) ) # bottleneck size at which max likelihood occurs
Max_LL_bottleneck <- bottleneck_values_vector[Max_LL_bottleneck_index] 
likelihood_ratio <- qchisq(.95, df=1) # necessary ratio of likelihoods set by confidence level
ci_tibble <- filter(LL_sum, 2*(Max_LL - adjusted_LL) <= likelihood_ratio ) 
lower_CI_bottleneck <- min(ci_tibble$bottleneck_size) #-1 # lower bound of confidence interval
upper_CI_bottleneck <- max(ci_tibble$bottleneck_size) #+1# upper bound of confidence interval


## merged Likelihoods



non_voc <- read.table ("Non-VOC/Non-VOC_merged_likelihood", header=T)
non_voc <- non_voc %>% group_by(bottleneck_size)%>% mutate (adjusted_LL= (max(Num_sites)/Num_sites)*Log_Likelihood)
LL_sum <- summarise (non_voc, sum (adjusted_LL))
LL_sum <- rename (LL_sum, adjusted_LL = 2)

bottleneck_values_vector <- c()
for ( i in 1:200)
{  if(i%%1 == 0) {bottleneck_values_vector <- c(bottleneck_values_vector,i)}
  
  
}


Max_LL <- max(LL_sum$adjusted_LL) # Maximum value of log likelihood
Max_LL_bottleneck_index <- which(LL_sum$adjusted_LL == max(LL_sum$adjusted_LL) ) # bottleneck size at which max likelihood occurs
Max_LL_bottleneck <- bottleneck_values_vector[Max_LL_bottleneck_index] 
likelihood_ratio <- qchisq(.95, df=1) # necessary ratio of likelihoods set by confidence level
ci_tibble <- filter(LL_sum, 2*(Max_LL - adjusted_LL) <= likelihood_ratio ) 
lower_CI_bottleneck <- min(ci_tibble$bottleneck_size) #-1 # lower bound of confidence interval
upper_CI_bottleneck <- max(ci_tibble$bottleneck_size) #+1# upper bound of confidence interval






alpha <- read.table ("Alpha/Alpha_merged_likelihood", header=T)
alpha <- alpha %>% group_by(bottleneck_size)%>% mutate (adjusted_LL= (max(Num_sites)/Num_sites)*Log_Likelihood)
LL_sum <- summarise (alpha, sum (adjusted_LL))
LL_sum <- rename (LL_sum, adjusted_LL = 2)

bottleneck_values_vector <- c()
for ( i in 1:200)
{  if(i%%1 == 0) {bottleneck_values_vector <- c(bottleneck_values_vector,i)}
  
  
}


Max_LL <- max(LL_sum$adjusted_LL) # Maximum value of log likelihood
Max_LL_bottleneck_index <- which(LL_sum$adjusted_LL == max(LL_sum$adjusted_LL) ) # bottleneck size at which max likelihood occurs
Max_LL_bottleneck <- bottleneck_values_vector[Max_LL_bottleneck_index] 
likelihood_ratio <- qchisq(.95, df=1) # necessary ratio of likelihoods set by confidence level
ci_tibble <- filter(LL_sum, 2*(Max_LL - adjusted_LL) <= likelihood_ratio ) 
lower_CI_bottleneck <- min(ci_tibble$bottleneck_size) #-1 # lower bound of confidence interval
upper_CI_bottleneck <- max(ci_tibble$bottleneck_size) #+1# upper bound of confidence interval


write.table 

delta <- read.table ("Delta/Delta_merged_likelihood", header=T)
delta <- delta %>% group_by(bottleneck_size)%>% mutate (adjusted_LL= (max(Num_sites)/Num_sites)*Log_Likelihood)
LL_sum <- summarise (delta, sum (adjusted_LL))
LL_sum <- rename (LL_sum, adjusted_LL = 2)

bottleneck_values_vector <- c()
for ( i in 1:200)
{  if(i%%1 == 0) {bottleneck_values_vector <- c(bottleneck_values_vector,i)}
  
  
}


Max_LL <- max(LL_sum$adjusted_LL) # Maximum value of log likelihood
Max_LL_bottleneck_index <- which(LL_sum$adjusted_LL == max(LL_sum$adjusted_LL) ) # bottleneck size at which max likelihood occurs
Max_LL_bottleneck <- bottleneck_values_vector[Max_LL_bottleneck_index] 
likelihood_ratio <- qchisq(.95, df=1) # necessary ratio of likelihoods set by confidence level
ci_tibble <- filter(LL_sum, 2*(Max_LL - adjusted_LL) <= likelihood_ratio ) 
lower_CI_bottleneck <- min(ci_tibble$bottleneck_size) #-1 # lower bound of confidence interval
upper_CI_bottleneck <- max(ci_tibble$bottleneck_size) #+1# upper bound of confidence interval







omicron<- read.table ("Omicron/Omicron_merged_likelihood", header=T)
omicron <- omicron %>% group_by(bottleneck_size)%>% mutate (adjusted_LL= (max(Num_sites)/Num_sites)*Log_Likelihood)
LL_sum <- summarise (omicron, sum (adjusted_LL))
LL_sum <- rename (LL_sum, adjusted_LL = 2)

bottleneck_values_vector <- c()
for ( i in 1:200)
{  if(i%%1 == 0) {bottleneck_values_vector <- c(bottleneck_values_vector,i)}
  
  
}


Max_LL <- max(LL_sum$adjusted_LL) # Maximum value of log likelihood
Max_LL_bottleneck_index <- which(LL_sum$adjusted_LL == max(LL_sum$adjusted_LL) ) # bottleneck size at which max likelihood occurs
Max_LL_bottleneck <- bottleneck_values_vector[Max_LL_bottleneck_index] 
likelihood_ratio <- qchisq(.95, df=1) # necessary ratio of likelihoods set by confidence level
ci_tibble <- filter(LL_sum, 2*(Max_LL - adjusted_LL) <= likelihood_ratio ) 
lower_CI_bottleneck <- min(ci_tibble$bottleneck_size) #-1 # lower bound of confidence interval
upper_CI_bottleneck <- max(ci_tibble$bottleneck_size) #+1# upper bound of confidence interval
