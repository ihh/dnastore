library("ggplot2");
df.l4 <- read.table(file="len4.mix2.dups.tab");
df.l8 <- read.table(file="len8.mix2.dups.tab");
df.l4$group <- 1;
df.l8$group <- 2;
df <- rbind(df.l4,df.l8);
df$Context=factor(df$group,labels=c("2 bases","4 bases"));
g <- guide_legend(title="Error model context");
dupplot <- ggplot(df,aes(x=DupProb,y=MeanEditsPerBit,color=Context)) + geom_errorbar(aes(ymin=MeanEditsPerBit-2*StDevEditsPerBit,ymax=MeanEditsPerBit+2*StDevEditsPerBit)) + geom_point(aes(shape=Context),size=2) + ylab("Edits per bit") + xlab("Duplication probability per base") + ggtitle("Effect on decoding accuracy of short (<5-base) duplications") + theme(legend.position=c(.2,.9)) + stat_smooth(method = "loess", formula = y ~ log(x), size = 1) + guides(colour=g,shape=g);
ggsave("dupplot.ps",plot=dupplot);

