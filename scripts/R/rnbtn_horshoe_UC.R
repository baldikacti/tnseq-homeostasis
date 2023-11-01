

library(doParallel)
library(tidyverse)
library(reshape)
library(parallel)
library( foreach)
library( forcats)
library(ParallelLogger)
library(brms)


## Read TN
tn <- read.table("../dlon-tnseq/data/processed/UC_mid_batchcorrected_dclpb_included.tsv",header=T)%>%
  filter(time!=0)

genewisefits <- function(tn, i, locuslist) {
  g <- locuslist[i]
  # Selecting  gene/locus tag from list
  ord_condition<- c("none","heat","oxidative-peroxide","canavanine")
  ord_strain <- c("wild-type","DLON","dnak-dnaJ","DCLPA","DCLPB")
  ord_slevel <- c("none","LOW","MEDIUM","HIGH")

  df_g <- tn  %>% dplyr::filter(locus_tag == g)%>%
    mutate(strain=fct_relevel(strain,ord_strain))%>%
    mutate(condition=fct_relevel(condition,ord_condition))%>%
    mutate(slevel=fct_relevel(slevel,ord_slevel))

  sdf <- summary(brm(tncnt ~ strain/condition/slevel,family=negbinomial, data=df_g,
                     cores=4,iter=2000,warmup = 1000,chains=4,
                     prior =prior(horseshoe(df = 2, par_ratio=0.5), class ="b")
  ))
  fdf <- data.frame(sdf$fixed,locus_tag=g,Effect=rownames(sdf$fixed),row.names = NULL)
  return(fdf)
}

# For Each gene run the model and store results in parallel
locusresults <-  list()
locuslist <- unique(tn$locus_tag)

## parallel runs

no_cores <- 3
cl <- parallel::makeCluster(no_cores, type = "FORK")
doParallel::registerDoParallel(cl)
cat("Running model in parallel.Running with",
    no_cores, "cores and", ctype, "\n")

`%dopar%` <- foreach::`%dopar%`
# do parallel function to apply rnbtn_model_pergene in parallel
locusresults <- foreach::foreach(i = 1 : length(locuslist))  %dopar%
  genewisefits(df, i, locuslist)

#stop cluster
ParallelLogger::stopCluster(cl)

## Post-Process and storing results in a data frame
#Here L1 is dummy header given by melt function
suppressWarnings(
  mod_data <- reshape::melt(locusresults) %>% dplyr::select(-L1))

#TABLE
write.table(mod_data,"UC_results.tab",quote = FALSE,row.names = FALSE)


