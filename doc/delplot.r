library("ggplot2");
df.l4.mix2 <- read.table(file="len4.mix2.dels.tab");
df.l4.sync64 <- read.table(file="len4.sync64.mix2.dels.tab");
df.l4.wat128 <- read.table(file="len4.wat128.ldpc.dels.tab");

df.l4.mix2$MutProb <- 1;
df.l4.sync64$MutProb <- 1;

df.l4.mix2$group <- 1;
df.l4.sync64$group <- 2;
df.l4.wat128$group <- 3;

df <- rbind(df.l4.mix2,df.l4.sync64,df.l4.wat128);
df <- df[df$DelProb > 0,];

df$Code=factor(df$group,labels=c("MIXRADAR(2) + DNASTORE(4)","SYNC(64) + MIXRADAR(2) + DNASTORE(4)","LDPC(16k,8k) + WMARK(128) + DNASTORE(4)"));
delplot <- ggplot(df,aes(x=DelProb,y=MeanEditsPerBit,color=Code)) + geom_errorbar(aes(ymin=MeanEditsPerBit-StDevEditsPerBit,ymax=MeanEditsPerBit+StDevEditsPerBit)) + geom_point(aes(shape=Code),size=2) + ylab("Edits per bit") + xlab("Deletion opening probability per base") + ggtitle("Effect on decoding accuracy of short (<5-base) deletions") + theme(legend.position=c(.73,.1)) + scale_x_log10();
ggsave("delplot.ps",plot=delplot);

