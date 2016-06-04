library("ggplot2");
df.l4 <- read.table(file="len4.dels.tab");
df.l4.mix2 <- read.table(file="len4.mix2.dels.tab");
df.l4.ham <- read.table(file="len4.ham.dels.tab");
df.l4.wat64 <- read.table(file="len4.wat64.1.ldpc.dels.tab");

df.l4$group <- 1;
df.l4.mix2$group <- 2;
df.l4.ham$group <- 3;
df.l4.wat64$group <- 4;

df <- rbind(df.l4,df.l4.mix2,df.l4.ham,df.l4.wat64);
df <- df[df$DelProb > 0,];

df$Code=factor(df$group,labels=c("DNASTORE(4)","MIXRADAR(2) + DNASTORE(4)","HAMMING(7,4) + DNASTORE(4)","LDPC(2k,1k) + WMARK(64) + DNASTORE(4)"));
delplot <- ggplot(df,aes(x=DelProb,y=MedianEditsPerBit,color=Code)) + geom_errorbar(aes(ymin=LowerQuartileEditsPerBit,ymax=UpperQuartileEditsPerBit)) + geom_point(aes(shape=Code),size=2) + ylab("Edits per bit") + xlab("Deletion opening probability per base") + ggtitle("Effect on decoding accuracy of short (<5-base) deletions") + theme(legend.position=c(.3,.85)) + scale_x_log10();
ggsave("delplot.ps",plot=delplot);

