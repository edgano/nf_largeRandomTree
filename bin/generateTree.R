
install.packages("seqinr")
library("seqinr")


fileName <- "PF00005_large_scale_seq.fa"
fastaFile <- read.fasta(file = fileName)

#how many fasta sequences
length(fastaFile)

#get the seqs IDs
seqNames<-getName(fastaFile)

##
##
##

install.packages("ape", dependencies = TRUE)
library(ape)

#define numbers of leave
n <- length(fastaFile)
#produce the tree
rt<-rmtree(1, n, rooted = TRUE, tip.label = seqNames, br = runif)

#print the tree
# str(rt)

#convert & save to newick
write.tree(rt, file = "PF00005_rnd.dnd", append = FALSE,digits = 2, tree.names = FALSE)
