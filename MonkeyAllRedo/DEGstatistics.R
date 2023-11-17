library(dplyr)
library(tibble)

# new Partek analyses
output.dir <- "~/MonkeyAllRego/outputs/"
deg.monkeys.cd3 <- read.csv(paste0(output.dir, "deg.monkeys.all.newCD3.cd3.csv")) %>% rename(Gene.name = X)
deg.monkeys.cd4 <- read.csv(paste0(output.dir, "deg.monkeys.all.newCD3.cd4.csv")) %>% rename(Gene.name = X)
deg.monkeys.cd8 <- read.csv(paste0(output.dir, "deg.monkeys.all.newCD3.cd8.csv")) %>% rename(Gene.name = X)

# compare to Partek analyses
input.dir <- "~/partekAnalysis/"
partek.cd4.log2fc.15 <- read.table(paste0(input.dir, "gene_list CD4 T Log2FC0.15.txt"), header = TRUE, sep = "\t") %>% select(-Gene.ID) %>% filter(Gene.name != ".")
partek.cd3.adjp.05 <- read.table(paste0(input.dir, "gene_list_AdjP0.05_Cd3.txt"), header = TRUE, sep = "\t") %>% select(-Gene.ID) %>% filter(Gene.name != ".") %>% mutate(log2FC = log2(Ratio..After.vs.Before.))
partek.cd3.log2fc.15 <- partek.cd3.adjp.05 %>% filter(log2FC >= 0.15 | log2FC <= -0.15)

shared.cd4 <- inner_join(partek.cd4.log2fc.15, deg.monkeys.cd4, by = "Gene.name") %>% mutate(log2FC_partek = log2(Ratio..After.vs.Before.)) %>% rename(log2FC_seurat = avg_log2FC) %>% select(Gene.name, log2FC_partek, log2FC_seurat)
shared.cd3.all <- inner_join(partek.cd3.adjp.05, deg.monkeys.cd3, by = "Gene.name") %>% mutate(log2FC_partek = log2(Ratio..After.vs.Before.)) %>% rename(log2FC_seurat = avg_log2FC) %>% select(Gene.name, log2FC_partek, log2FC_seurat)
shared.cd3.all.signif <- shared.cd3.all %>% filter(log2FC_partek >= 0.15 | log2FC_partek <= -0.15) %>% filter(log2FC_seurat >= 0.5 | log2FC_seurat <= -0.5)
shared.cd3.all.signagree <- shared.cd3.all %>% filter(log2FC_partek * log2FC_seurat > 0)

write.csv(shared.cd3.all.signif, paste0(output.dir, "sharedwParteksignif.cd3.csv"))
write.csv(shared.cd3.all.signagree, paste0(output.dir, "sharedwParteksignagree.cd3.csv"))
write.csv(shared.cd4, paste0(output.dir, "sharedwParteksignif.cd4.csv"))

#CD8s
cd8.vs.CD4 <- inner_join(deg.monkeys.cd8, deg.monkeys.cd4, by = "Gene.name") %>% rename(log2FC_CD8 = avg_log2FC.x, log2FC_CD4 = avg_log2FC.y) %>% select(Gene.name, log2FC_CD8, log2FC_CD4)
cd8.downreg <- deg.monkeys.cd8 %>% filter(avg_log2FC <= -0.5)
cd4.downreg <- deg.monkeys.cd4 %>% filter(avg_log2FC <= -0.5)
cd3.downreg <- deg.monkeys.cd3 %>% filter(avg_log2FC <= -0.5)
temp <- cd8.vs.CD4 %>% filter(Gene.name %in% listfrompartek)

log2fc8 <- deg.monkeys.cd8 %>% filter(avg_log2FC <= -0.5 | avg_log2FC >= 0.5)
log2fc4 <- deg.monkeys.cd4 %>% filter(avg_log2FC <= -0.5 | avg_log2FC >= 0.5)
log2fc3 <- deg.monkeys.cd3 %>% filter(avg_log2FC <= -0.5 | avg_log2FC >= 0.5)
shared.cd4.signif <- shared.cd4 %>% filter(log2FC_seurat <= -0.5 | log2FC_seurat >= 0.5)

#AP1 members
listfrompartek <- c("MAF", "MAFF", "ATF3", "BATF", "FOS", "FOSL2", "FOSB", "JUN", "JUNB")
cd4.AP1members <- deg.monkeys.cd4 %>% filter(Gene.name %in% listfrompartek) %>% select(Gene.name, avg_log2FC)
cd8.AP1members <- deg.monkeys.cd8 %>% filter(Gene.name %in% listfrompartek) %>% select(Gene.name, avg_log2FC)
partek.AP1members <- shared.cd4 %>% filter(Gene.name %in% listfrompartek)
