# tree map for DEGs
library(data.table)
library(treemap)
library(UpSetR)

extract_sign <- function(f){
    tab <- fread(f)
    tab_sign <- subset(tab, padj<0.01)
    sample_set <- tab_sign$Id
    n <- length(sample_set)
}

basename <- function(name, type){
  name <- strsplit(name, split=type)[[1]][1]
  name_split <- strsplit(name, split="vs")
  name_mod <- paste(name_split[[1]][1], name_split[[1]][2], sep="_vs_")
}


interlace <- function(l1, l2, l3, common_len){
  full_list <- c()
  for (i in seq(1, common_len)){
    full_list <- c(l1[i], l2[i], l3[i], full_list)
  }
  return(full_list)
}

comb_subgroups <- function(names, vals){
  sub <- c()
  for (i in seq(1,length(names))){
    n <- paste(names[i], vals[i], sep="\n")
    sub <- c(sub, n)
  }
  return(sub)
}

# treemap
setwd("~/Desktop/report/DGE/tables/")

files_up <- list.files(pattern="up")
files_down <- list.files(pattern="down")
n_up <- sapply(files_up, extract_sign)
n_down <- sapply(files_down, extract_sign)
n_nondiff <- 13313 - (n_up + n_down)

conds <- sapply(names(n_up), basename, type=".up.")
groups <- c()
for (i in conds){groups <- c(rep(i,3), groups)}
subgroup_names <- rep(c("nondiff", "up", "down"), 6)
vals <- unname(interlace(n_nondiff, n_up, n_down, 6))
subgroups <- comb_subgroups(names = subgroup_names, vals=vals)

data <- data.frame(groups, subgroups, vals)

plot(treemap(data, index=c("groups", "subgroups"), vSize="vals", type="index",
        title="Number of unaffected, up- and down- regulated genes",
        fontcolor.labels = "black", border.lwds=0, bg.labels=255, aspRatio = 1.2,
        ymod.labels=c(1.2,0), fontface.labels=c("italic", "plain")))


# upsetr
setwd("~/Desktop/report/DGE/tables/")
files <- list.files(pattern=".complete.")
deg <- sapply(files, extract_sign)
n <- sapply(names(deg), basename, type=".complete.")
names(deg) <- n
upset(fromList(deg), nsets = 6, order.by = "freq")



