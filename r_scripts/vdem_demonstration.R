# This is a demonstration of 2-way and one-way FX analysis using data provided by the varities of democracy
# Robert Kubinec 7/12/2016

require(haven)
require(dplyr)
require(magrittr)
require(sandwich)
require(tibble)
require(ggplot2)
require(stringr)

#Run once, use a faster R version
#full_dataset <- read_dta("C:/Users/bobku/Box Sync/Between Effects/V-DEM/Country_Year_V-Dem_other_STATA_v6.1/V-Dem-DS-CY+Others-v6.1.dta")
#full_dataset <- saveRDS(object = full_dataset,"data/vdem_full.rds")
# vdem_pos <- data.table::fread('data/polyarchy.csv')

full_dataset <- readRDS('data/vdem_full.rds')

full_dataset$year_factor <- as.factor(full_dataset$year)


# Load posterior estimates of the V-DEM indices

vdem_pos <- as.tbl(as.data.frame(readRDS('data/vdem_pos.rds')))
vdem_pos$V1 <- NULL

# combined posterior estimates with posterior uncertainty
merged_data <- left_join(full_dataset,vdem_pos,by=c('country_text_id','year'))
to_analyze <- merged_data[,c('e_migdppcln','year_factor','country_name')]

models <- c('Within','Between','Two-way')
varnames <- paste0('V',2:900)

over_posterior <- function(x,y) {
  to_analyze$v2x_polyarchy <- merged_data[[x]]
  model1 <- lm(formula = y,data = to_analyze)
  # Use the sandwich estimator to adjust variances
  # Need to drop coefs that come out as NA
  # This happens with the country of the Republic of Vietnam in the within-between model (model 5)
  # The sandwich estimator will automatically drop NA coefficients from the VCOV matrix, causing an error
  # With MASS
  coefs <- coef(model1)[!is.na(coef(model1))]
  sds <- vcovHC(model1,type='HC0')
  
  samples <- MASS::mvrnorm(mu=coefs,Sigma=sds)
  return(samples)
  }

# Within effect

model1 <- sapply(varnames,over_posterior,v2x_polyarchy ~ e_migdppcln + country_name)
output_model1 <- data_frame(betas=row.names(model1),coef=apply(model1,1,mean),
                        sd=apply(model1,1,sd))
# Between effect

model2 <- sapply(varnames,over_posterior,v2x_polyarchy ~ e_migdppcln + year_factor)
output_model2 <- data_frame(betas=row.names(model2),coef=apply(model2,1,mean),
                            sd=apply(model2,1,sd))
# 2-way FX

model3 <- sapply(varnames,over_posterior,v2x_polyarchy ~ e_migdppcln + year_factor + country_name)
output_model3 <- data_frame(betas=row.names(model3),coef=apply(model3,1,mean),
                            sd=apply(model3,1,sd))


# All three are statistically significant, the within effect is the largest

# Combine estimates


# Try a few more exotic models

# Between effect that varies over time
# Note that this model implicitly accounts for autocorrelation

model4 <- sapply(varnames,over_posterior,v2x_polyarchy ~ e_migdppcln + year_factor + e_migdppcln*year_factor)
output_model4 <- data_frame(betas=row.names(model4),coef=apply(model4,1,mean),
                            sd=apply(model4,1,sd))


# Within effect that varies between

model5 <- sapply(varnames,over_posterior,v2x_polyarchy ~ e_migdppcln + country_name + e_migdppcln*country_name)
output_model5 <- data_frame(betas=row.names(model5),coef=apply(model5,1,mean),
                            sd=apply(model5,1,sd))

#Models with controls
to_analyze <- merged_data[,c('e_migdppcln','year_factor','country_name','e_WORLD_DEM_DIFFUSE',
                             'e_Vanhanen_studentpercent_ipo','e_Vanhanen_familyfarm_ipo',
                             'e_cap_share_unequal','e_Total_Fuel_Income_PC',
                             'e_Vanhanen_familyfarm_ipo','e_miferrat','e_Civil_War','e_Vanhanen_urban_ipo')]

model6 <- sapply(varnames,over_posterior,v2x_polyarchy ~ e_migdppcln + e_WORLD_DEM_DIFFUSE + 
                e_cap_share_unequal + e_Total_Fuel_Income_PC  + country_name)

output_model6 <- data_frame(betas=row.names(model6),coef=apply(model6,1,mean),
                            sd=apply(model6,1,sd))

model7 <- sapply(varnames,over_posterior,v2x_polyarchy ~ e_migdppcln +
                    e_cap_share_unequal + e_Total_Fuel_Income_PC  + 
                   e_miferrat + e_Civil_War + year_factor)
output_model7 <- data_frame(betas=row.names(model7),coef=apply(model7,1,mean),
                            sd=apply(model7,1,sd))


# Calculate interaction effects and plot

model5 <- readRDS("data/model5.rds") %>% t %>% as.data.frame %>% tbl_df
countries <- select(model5,matches(":country"))
int_matrix <- tbl_df(lapply(countries,function(x) x + model5[['e_migdppcln']]))
results <- data_frame(Coef=sapply(int_matrix,mean), SD = sapply(int_matrix,sd))
results$variables <- str_extract(colnames(countries),'[A-Z].+')
results <- mutate(results,upper=Coef + 1.96*SD,lower=Coef - 1.96*SD)

ggplot(results,aes(x=Coef,y=reorder(variables,Coef))) + geom_point() + ggnetwork::theme_blank() + geom_text(aes(label=variables),hjust=-1,check_overlap=TRUE) +
  geom_vline(xintercept=0) + geom_errorbarh(aes(xmin=lower,xmax=upper))

model5 <- readRDS("data/model5.rds") %>% t %>% as.data.frame %>% tbl_df
countries <- select(model5,matches(":country"))
int_matrix <- tbl_df(lapply(countries,function(x) x + model5[['e_migdppcln']]))
results <- data_frame(Coef=sapply(int_matrix,mean), SD = sapply(int_matrix,sd))
results$variables <- str_extract(colnames(countries),'[A-Z].+')
results <- mutate(results,upper=Coef + 1.96*SD,lower=Coef - 1.96*SD)

ggplot(results,aes(x=Coef,y=reorder(variables,Coef))) + geom_point() + ggnetwork::theme_blank() + geom_text(aes(label=variables),hjust=-1,check_overlap=TRUE) +
  geom_vline(xintercept=0) + geom_errorbarh(aes(xmin=lower,xmax=upper))


 