

library(dplyr)
library(brms)


## Read TN
tn <- read.table("TC_mid_batchcorrected_dclpb_included.tsv",header=T)%>%
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


mod_data <- genewisefits(df, i, locuslist)


#WRITE
write.table(mod_data,paste0(i,"_results.tab"),quote = FALSE,row.names = FALSE)


