library("ggplot2");
df.l4 <- read.table(file="len4.mix2.dups.tab");
df.l8 <- read.table(file="len8.mix2.dups.tab");

df.l4$group <- 1;
df.l8$group <- 2;

df <- rbind(df.l4,df.l8);

df$Context=factor(df$group,labels=c("2 bases","4 bases"));
g <- guide_legend(title="Error model context");

dupplot <- ggplot(df,aes(x=DupProb,y=MedianEditsPerBit,color=Context)) + geom_errorbar(aes(ymin=LowerQuartileEditsPerBit,ymax=UpperQuartileEditsPerBit)) + geom_point(aes(shape=Context),size=2) + ylab("Edits per bit") + xlab("Duplication probability per base") + ggtitle("Effect on decoding accuracy of short (<5-base) duplications") + theme(legend.position=c(.2,.9)) + guides(colour=g,shape=g) + scale_x_log10();

ggsave("dupplot.ps",plot=dupplot);

