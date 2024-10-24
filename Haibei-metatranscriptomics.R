## treatment effects by linear mixed model
## the precipitation level is noted as -0.5, 0, 0.5 in LMMs.
## the warming is noted as 0, 1 in LMMs.

library(lme4)
library(car)
library(vegan)
## environmental factors
envs_2020 <- read.csv('envs_2020.csv')
var <- c('Temperature', 'Moisture', 'pH', 'MBC', 'NEE', 'ER', 'CH4', 'Rs', 
         'Rh', 'CH4/Rh')

for (j in var){
    fm1 <- lmer(envs_2020[,j] ~ Warm*Pre + (1|block), data = envs_2020)
    presult <- car::Anova(fm1,type=2) %>% as.data.frame()
    lmm <- coef(summary(fm1))[2:4, ] %>% as.data.frame() 
    result <- cbind(lmm, presult)
    result$treat <- row.names(result)
    result$var <- colnames(envs_2020)[j]
    colnames(result) <- c('Estimate', 'SE', 'tvalue', 'chisq', 'df', 'Pvalue', 'treat', 'var')
    lmm_ENVs <- rbind(lmm_ENVs, result)
  }


## active soil microbes and metabolism
CPM_microbes <- read.csv('CPM_microbes.csv')
plot_information <- read.csv('plot_information.csv')
CPM_microbes_stand <- decostand(CPM_microbes, method = 'standardize', MARGIN = 2)
CPM_microbes_stand <- merge(plot_information, CPM_microbes_stand,
                            by.x = 'plot', by.y = 'row.names', all = T)
gene <- colnames(CPM_microbes)
  for (j in gene){
    fm1 <- lmer(CPM_microbes_stand[,j] ~ Warm*Pre + (1|block), 
                data = CPM_microbes_stand)
    presult <- car::Anova(fm1,type=2) %>% as.data.frame()
    lmm <- coef(summary(fm1))[2:4, ] %>% as.data.frame() 
    result <- cbind(lmm, presult)
    result$treat <- row.names(result)
    result$var <- colnames(CPM_microbes_stand)[j]
    colnames(result) <- c('Estimate', 'SE', 'tvalue', 'chisq', 'df', 'Pvalue', 'treat', 'var')
    lmm_CPM_microbes <- rbind(lmm_CPM_microbes, result)
  }

## Correlations between microbial metabolism and soil carbon fluxes
all_data <- merge(envs_2020, CPM_microbes, by.x = 'plot', by.y = 'row.names', all = T)
fluxes <- c('CH4','Rh', 'CH4/Rh')
  for (i in gene) {
    for (j in fluxes){
    fm1 <- lmer(all_data[,j] ~ all_data[,i] + (1|block), data = all_data)
    lmm <- coef(summary(fm1))[2, ]
    presult <- car::Anova(fm1,type=2) %>% as.data.frame()
    r2 <- r.squaredGLMM(fm1)[,'R2m']; names(r2) <- c('r2')
    r <- ifelse(lmm[1]>0,(r2)^0.5,-(r2)^0.5); names(r) <- c('r')
    res <- c(lmm, r2, r); res <- as.data.frame(res) %>% t() %>% as.data.frame()
    result <- cbind(presult, res)
    result$Y_var <- j
    result$X_var <- i
    colnames(result) <- c('chisq', 'df', 'Pvalue', 'Estimate', 'SE', 'tvalue','r2', 'r', 'Y_var', 'X_var')
    lmm_cor <- rbind(lmm_cor, result)}}

## multiple regression
## SAR refers to SAR-Transporter and SAR-Catabolism
all_data_stand <- decostand(all_data[,c(gene, fluxes)], method = 'standardize', MARGIN = 2) %>% 
    cbind(all_data[1:10], .)
mulmm_flux_ratio <- lmer(CH4_Rh ~ Moisture+Temperature+CAZy+SAR+(1|block),
                         data = all_data_stand)


## PCoA and PERMANOVA
library(ieggr)
library(ape)
species <- read.csv('species.csv')
beta.weight <- beta.g(species, dist.method = c("bray"), abundance.weighted = TRUE)
PCoA_weight <- pcoa(beta.weight, correction = "none")
PCoA_weight_vector <- as.data.frame(PCoA_weight$vectors)[,1:4] %>% 
  merge(., plot_information, by.x = 'row.names', by.y = 'plot', all = F)
set.seed(123)
adonis.webray <- adonis2(beta.weight~Warm*Pre, data=PCoA_weight_vector)





