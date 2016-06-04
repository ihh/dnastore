library("ggplot2");
df.wat64 <- read.table(file="len4.wat64.ldpc.dels.tab");
df.wat64.1 <- read.table(file="len4.wat64.1.ldpc.dels.tab");
df.wat64.16 <- read.table(file="len4.wat64.16.ldpc.dels.tab");
df.wat128.1 <- read.table(file="len4.wat128.1.ldpc.dels.tab");

df.wat64$group <- 1;
df.wat64.16$group <- 2;
df.wat64.1$group <- 3;
df.wat128.1$group <- 4;

df <- rbind(df.wat64,df.wat64.16,df.wat64.1,df.wat128.1);
df <- df[df$DelProb > 0,];

df$Code=factor(df$group,labels=c("64-bit markers, no watermark bits (radix signature only)","64-bit markers, 1 watermark bit/16 signal bits","64-bit markers, 1 watermark bit/signal bit","128-bit markers, 1 watermark bit/signal bit"));
g <- guide_legend(title="Synchronization code");

watplot <- ggplot(df,aes(x=DelProb,y=MedianEditsPerBit,color=Code)) + geom_errorbar(aes(ymin=LowerQuartileEditsPerBit,ymax=UpperQuartileEditsPerBit)) + geom_point(aes(shape=Code),size=2) + ylab("Edits per bit") + xlab("Deletion opening probability per base") + ggtitle("Effect on performance of watermark scheme") + theme(legend.position=c(.3,.85)) + scale_x_log10() + guides(colour=g,shape=g);

ggsave("watplot.ps",plot=watplot);

