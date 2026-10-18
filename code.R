library(TwoSampleMR)
library(tidyverse)

exposure_dat <- extract_instruments(c('ebi-a-GCST90038604',
                                      'ebi-a-GCST90014023',
                                      'ebi-a-GCST006867'))

outcome_dat <- read_outcome_data(filename = 'finngen_R13_G6_EPLEPSY.gz',
                                 sep = '\t',
                                 snps = exposure_dat$SNP,
                                 snp_col = 'rsids',
                                 beta_col = 'beta',
                                 se_col = 'sebeta',
                                 eaf_col = 'af_alt',
                                 effect_allele_col = 'alt',
                                 other_allele_col = 'ref',
                                 pval_col = 'pval',
                                 chr_col = '#chrom',
                                 pos_col = 'pos')
gc()

dat <- harmonise_data(exposure_dat = exposure_dat,
                      outcome_dat = outcome_dat)
dat$outcome <- 'Epilepsy'

# Hypertension ------------------------------------------------------------

df_1 <- dat |> 
  filter(exposure == 'Hypertension || id:ebi-a-GCST90038604')
df_1$exposure <- 'Hypertension'
mr_dat_1 <- MendelianRandomization::mr_input(bx=df_1$beta.exposure, bxse = df_1$se.exposure,
                                             by = df_1$beta.outcome, byse = df_1$se.outcome,
                                             snps = df_1$SNP, effect_allele = df_1$effect_allele.exposure,
                                             other_allele = df_1$other_allele.exposure,
                                             eaf = df_1$eaf.exposure
)

(r <- MendelianRandomization::mr_ivw(mr_dat_1,
                                     model = 'random'))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError), 3)
round(exp(r@ Estimate+1.96*r@ StdError), 3)

(r <- MendelianRandomization::mr_egger(mr_dat_1))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError.Est), 3)
round(exp(r@ Estimate+1.96*r@ StdError.Est), 3)

(r <- MendelianRandomization::mr_cML(mr_dat_1, n = 382851))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError), 3)
round(exp(r@ Estimate+1.96*r@ StdError), 3)

MRPRESSO::mr_presso(BetaExposure = 'beta.exposure', BetaOutcome = 'beta.outcome',
                    SdExposure = 'se.exposure', SdOutcome = 'se.outcome',
                    data = df_1)

# Type 1 diabetes ---------------------------------------------------------
df_2 <- dat |> 
  filter(exposure == 'Type 1 diabetes || id:ebi-a-GCST90014023')
df_2$exposure <- 'Type 1 diabetes'
mr_dat_2 <- MendelianRandomization::mr_input(bx=df_2$beta.exposure, bxse = df_2$se.exposure,
                                             by = df_2$beta.outcome, byse = df_2$se.outcome,
                                             snps = df_2$SNP, effect_allele = df_2$effect_allele.exposure,
                                             other_allele = df_2$other_allele.exposure,
                                             eaf = df_2$eaf.exposure)

(r <- MendelianRandomization::mr_ivw(mr_dat_2,
                                     model = 'fixed'))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError), 3)
round(exp(r@ Estimate+1.96*r@ StdError), 3)

(r <- MendelianRandomization::mr_egger(mr_dat_2))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError.Est), 3)
round(exp(r@ Estimate+1.96*r@ StdError.Est), 3)

(r <- MendelianRandomization::mr_cML(mr_dat_2, n = 382851))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError), 3)
round(exp(r@ Estimate+1.96*r@ StdError), 3)

MRPRESSO::mr_presso(BetaExposure = 'beta.exposure', BetaOutcome = 'beta.outcome',
                    SdExposure = 'se.exposure', SdOutcome = 'se.outcome',
                    data = df_2)

# Type 2 diabetes ---------------------------------------------------------
df_3 <- dat |> 
  filter(exposure == 'Type 2 diabetes || id:ebi-a-GCST006867')
df_3$exposure <- 'Type 2 diabetes'
mr_dat_3 <- MendelianRandomization::mr_input(bx=df_3$beta.exposure, bxse = df_3$se.exposure,
                                             by = df_3$beta.outcome, byse = df_3$se.outcome,
                                             snps = df_3$SNP, effect_allele = df_3$effect_allele.exposure,
                                             other_allele = df_3$other_allele.exposure,
                                             eaf = df_3$eaf.exposure
)

(r <- MendelianRandomization::mr_ivw(mr_dat_3,
                                     model = 'fixed'))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError), 3)
round(exp(r@ Estimate+1.96*r@ StdError), 3)

(r <- MendelianRandomization::mr_egger(mr_dat_3))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError.Est), 3)
round(exp(r@ Estimate+1.96*r@ StdError.Est), 3)

(r <- MendelianRandomization::mr_cML(mr_dat_3, n = 382851))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError), 3)
round(exp(r@ Estimate+1.96*r@ StdError), 3)

MRPRESSO::mr_presso(BetaExposure = 'beta.exposure', BetaOutcome = 'beta.outcome',
                    SdExposure = 'se.exposure', SdOutcome = 'se.outcome',
                    data = df_3)

snp <- data.table::fread('clumped_out.clumped')$SNP
exposure_dat <- format_data(dat = data.table::fread('finngen_R13_G6_EPLEPSY') |> 
                              filter(rsids %in% snp) |> 
                              as.data.frame(),
                            type = "exposure",
                            snp_col = 'rsids',
                            beta_col = 'beta',
                            se_col = 'sebeta',
                            eaf_col = 'af_alt',
                            effect_allele_col = 'alt',
                            other_allele_col = 'ref',
                            pval_col = 'pval',
                            chr_col = '#chrom',
                            pos_col = 'pos')
outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP,
                                    outcomes = c('ebi-a-GCST90038604',
                                                 'ebi-a-GCST90014023',
                                                 'ebi-a-GCST006867'))
dat <- harmonise_data(exposure_dat,
                      outcome_dat)
dat$exposure <- 'Epilepsy'
# Hypertension ------------------------------------------------------------

df_1 <- dat |> 
  filter(outcome == 'Hypertension || id:ebi-a-GCST90038604')
df_1$outcome <- 'Hypertension'
mr_dat_1 <- MendelianRandomization::mr_input(bx=df_1$beta.exposure, bxse = df_1$se.exposure,
                                             by = df_1$beta.outcome, byse = df_1$se.outcome,
                                             snps = df_1$SNP, effect_allele = df_1$effect_allele.exposure,
                                             other_allele = df_1$other_allele.exposure,
                                             eaf = df_1$eaf.exposure)

(r <- MendelianRandomization::mr_ivw(mr_dat_1,
                                     model = 'random'))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError), 3)
round(exp(r@ Estimate+1.96*r@ StdError), 3)

(r <- MendelianRandomization::mr_egger(mr_dat_1))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError.Est), 3)
round(exp(r@ Estimate+1.96*r@ StdError.Est), 3)

(r <- MendelianRandomization::mr_cML(mr_dat_1, n = 382851))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError), 3)
round(exp(r@ Estimate+1.96*r@ StdError), 3)

MRPRESSO::mr_presso(BetaExposure = 'beta.exposure', BetaOutcome = 'beta.outcome',
                    SdExposure = 'se.exposure', SdOutcome = 'se.outcome',
                    data = df_1)

# Type 1 diabetes ---------------------------------------------------------
df_2 <- dat |> 
  filter(outcome == 'Type 1 diabetes || id:ebi-a-GCST90014023')
df_2$outcome <- 'Type 1 diabetes'

mr_dat_2 <- MendelianRandomization::mr_input(bx=df_2$beta.exposure, bxse = df_2$se.exposure,
                                             by = df_2$beta.outcome, byse = df_2$se.outcome,
                                             snps = df_2$SNP, effect_allele = df_2$effect_allele.exposure,
                                             other_allele = df_2$other_allele.exposure,
                                             eaf = df_2$eaf.exposure)

(r <- MendelianRandomization::mr_ivw(mr_dat_2,
                                     model = 'random'))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError), 3)
round(exp(r@ Estimate+1.96*r@ StdError), 3)

(r <- MendelianRandomization::mr_egger(mr_dat_2))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError.Est), 3)
round(exp(r@ Estimate+1.96*r@ StdError.Est), 3)

(r <- MendelianRandomization::mr_cML(mr_dat_2, n = 382851))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError), 3)
round(exp(r@ Estimate+1.96*r@ StdError), 3)

MRPRESSO::mr_presso(BetaExposure = 'beta.exposure', BetaOutcome = 'beta.outcome',
                    SdExposure = 'se.exposure', SdOutcome = 'se.outcome',
                    data = df_2)

# Type 2 diabetes ---------------------------------------------------------
df_3 <- dat |> 
  filter(outcome == 'Type 2 diabetes || id:ebi-a-GCST006867')
df_3$outcome <- 'Type 2 diabetes'
mr_dat_3 <- MendelianRandomization::mr_input(bx=df_3$beta.exposure, bxse = df_3$se.exposure,
                                             by = df_3$beta.outcome, byse = df_3$se.outcome,
                                             snps = df_3$SNP, effect_allele = df_3$effect_allele.exposure,
                                             other_allele = df_3$other_allele.exposure,
                                             eaf = df_3$eaf.exposure
)

(r <- MendelianRandomization::mr_ivw(mr_dat_3,
                                     model = 'random'))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError), 3)
round(exp(r@ Estimate+1.96*r@ StdError), 3)

(r <- MendelianRandomization::mr_egger(mr_dat_3))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError.Est), 3)
round(exp(r@ Estimate+1.96*r@ StdError.Est), 3)

(r <- MendelianRandomization::mr_cML(mr_dat_3, n = 382851))
round(exp(r@ Estimate), 3)
round(exp(r@ Estimate-1.96*r@ StdError), 3)
round(exp(r@ Estimate+1.96*r@ StdError), 3)

MRPRESSO::mr_presso(BetaExposure = 'beta.exposure', BetaOutcome = 'beta.outcome',
                    SdExposure = 'se.exposure', SdOutcome = 'se.outcome',
                    data = df_3)

source('mr_forest_plot.R')
mr_forest_plot(mr_singlesnp(df_1))
mr_forest_plot(mr_singlesnp(df_2))
mr_forest_plot(mr_singlesnp(df_3))

