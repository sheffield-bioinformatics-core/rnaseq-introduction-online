## check if packages are found and install if not

if(!require(BiocManager)) install.packages("BiocManager")
if(!require(tidyverse)) install.packages("tidyverse")  
if(!require(org.Mm.eg.db)) BiocManager::install("org.Mm.eg.db")


## load our packages
library(org.Mm.eg.db)
library(tidyverse)

## find genes belonging to GO:0007049  - cell cycle
my_genes <- AnnotationDbi::select(org.Mm.eg.db, 
              keys = "GO:0007049",
              keytype = "GO",
              columns = c("SYMBOL","ENTREZID")) %>% pull(SYMBOL)

length(my_genes)
## 634

## load our entire set of results
results <- read_csv("background.csv")
n_sig <- sum(results$FDR < 0.05)

# 4595

n_obs <- sum(results$FDR < 0.05 & results$SYMBOL %in% my_genes)
## 233

tab <- table(results$SYMBOL %in% my_genes, results$FDR < 0.05)
tab
#
#         FALSE  TRUE
# FALSE 22196  4362
# TRUE    388   233


chisq.test(tab)

test_genes <- NULL

## set seed so the results are reproducible
set.seed(1233)

for(i in 1:1000){
  
  test_genes[[i]] <- results %>% 
    slice_sample(n = n_sig) %>% ## pick n_sig rows at random
    mutate(InPathway = SYMBOL %in% my_genes) %>% ## add extra column for how many of selected genes are in our gene set
    count(InPathway) %>% ## count up the number in our gene set
    filter(InPathway) %>% 
    pull(n) ## extract n - the number of genes in our gene set for this random set
}

data.frame(N = unlist(test_genes)) |>
  ggplot(aes(x  = N)) + geom_histogram(binwidth = 1) + geom_vline(xintercept = n_obs,col="red",lty=2) 
ggsave("geneset_randomHist.png")
