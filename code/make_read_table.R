np <- function(ct, col){
    all <- nrow(ct)
    n <- all - table(ct[,col])[[1]]
    p <- round(n/all * 100, 2)
    return(paste(n, p, sep=" "))
}


assigned <- function(col){
    return(sum(col))
}

library(matrixStats)

read_stats <- read.csv("~/Schreibtisch/read_stats", header=FALSE, sep=' ')
names(read_stats) <- c("Sample ID", "raw reads", "uniquely mapped", "um%", "multi mapped", "mm%", "assigned reads")
read_stats$`Sample ID` <- as.character(read_stats$`Sample ID`)

arp <- (read_stats$`assigned reads`/read_stats$`raw reads`)*100
read_stats$`ar%` <- arp

means <- c("Mean", round(colMeans(read_stats[,2:length(read_stats)]),2))
read_stats <- rbind(read_stats, means)

mat <- as.matrix(read_stats[1:nrow(read_stats)-1,2:length(read_stats)])
m <- mapply(mat, FUN=as.numeric)
m2 <- matrix(data=m, ncol=7, nrow=12)
sds <- c("S.d.", round(colSds(m2), 2) )

read_stats <- rbind(read_stats, sds)
read_stats

write.csv(read_stats, "~/Schreibtisch/read_stats_full")
