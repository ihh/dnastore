library("ggplot2");
df.l4 <- read.table(file="len4.subs.tab");
df.l4.ham <- read.table(file="len4.ham.subs.tab");
df.l4.mix2 <- read.table(file="len4.mix2.subs.tab");
df.l4.ldpc <- read.table(file="len4.wat128.ldpc.subs.tab");

df.l4$group <- 1;
df.l4.mix2$group <- 2;
df.l4.ldpc$group <- 3;
df.l4.ham$group <- 4;

df <- rbind(df.l4,df.l4.mix2,df.l4.ldpc,df.l4.ham);

df$Code=factor(df$group,labels=c("DNASTORE(4)","MIXRADAR(2) + DNASTORE(4)","LDPC(16k,8k) + WMARK(128) + DNASTORE(4)","HAMMING(7,4) + DNASTORE(4)"));
subplot <- ggplot(df,aes(x=SubProb,y=MedianEditsPerBit,color=Code)) + geom_errorbar(aes(ymin=LowerQuartileEditsPerBit,ymax=UpperQuartileEditsPerBit)) + geom_point(aes(shape=Code),size=2) + ylab("Edits per bit") + xlab("Substitution probability per base") + ggtitle("Effect on decoding accuracy of point substitutions") + theme(legend.position=c(.25,.85)) + scale_x_log10();

ggsave("subplot.ps",plot=subplot);

