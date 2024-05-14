## change log10PVAL to PVAL for manhattan plot

args = commandArgs(trailingOnly=TRUE)

fn=args[1]
output=args[2]

mt <- read.table(fn, header=TRUE,fill=TRUE)

mt$LOG10P <- as.numeric(as.character(mt$LOG10P))
mt$PVAL <- 10**(-mt$LOG10P)

write.table(mt, output, row.names = FALSE, sep = '\t', quote=FALSE) 
