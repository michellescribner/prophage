### Analysis of coverage for populations and clones with pf5r mutations


library(tidyverse)
library(ggpubr)
library(extrafont)

library(extrafontdb)

library(gridExtra)
library(Sushi)
library(ggrepel)
library(ggsci)


setwd("/Users/mrs/Documents/prophage/coverage/output")

### Plot Coverage

# Planktonic Population 5
anc <- read.delim("anccov_10.txt", header = FALSE)
anc$strain <- "Ancestral Strain"
p5 <- read.delim("P5cov_10.txt", header = FALSE)
p5$strain <- "Planktonic Population 5"
all <- rbind(p5, anc)
colnames(all) <- c("chrom", "Position", "end", "ReadDepth", "Strain")

p <- ggplot(all, aes(x = Position, y=ReadDepth, color=Strain)) + 
  geom_line(aes(color=Strain)) + theme_classic() + 
  scale_y_continuous(limits=c(0,9500)) + 
  scale_x_continuous(breaks=c(0, 6537648), labels = c("0", "6,537,648 bp")) +
  scale_color_manual(values=c("dark gray", "#4DBBD5FF")) +
  theme(axis.title.x = element_blank(), 
        axis.text = element_text(size = 10, color = "black"),
        axis.line = element_line(size=0.25), 
        legend.position = c(0.85,0.9),
        legend.title = element_blank(),
        plot.margin = margin(10, 20, 10, 10)) + labs(title="PA14 Genome")

all_p <- all[which(all$Position > 4335000), ]
all_p <- all_p[which(all_p$Position < 4360000), ]

p_p <- ggplot() + 
  geom_line(data = all_p, aes(x = Position, y=ReadDepth, color=Strain)) + 
  scale_y_continuous(limits=c(0,9500)) + 
  scale_x_continuous(breaks=c(4335000, 4360000), labels = c("4,335,000 bp", "4,360,000 bp")) +
  scale_color_manual(values=c("dark gray", "#4DBBD5FF")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 10, color = "black"), 
        axis.line = element_line(size=0.25), 
        plot.margin = margin(10, 20, 30, 10),
        axis.title.x = element_blank()) + 
  labs(title = "Pf5") +
  geom_segment(aes(x = 4345126, y = 8500, xend =4345126, yend=9500), color = "purple", size =1) +
  geom_segment(aes(x = 4355790, y = 8500, xend =4355790, yend=9500), color = "purple", size=1) +
  geom_segment(aes(x = 4353517, y = 8500, xend =4353517, yend=9500), color = "red", size=1) +
  geom_segment(aes(x = 4353532, y = 8500, xend =4353532, yend=9500), color = "red", size =1) 

p5 <- ggarrange(p, p_p, labels = c("A", "B"),
                ncol = 1, nrow = 2)
ggsave("PlanktonicPopulation5.png", p5, width = 7, height = 3.5)

# Biofilm Population 1 
anc <- read.delim("anccov_10.txt", header = FALSE)
anc$strain <- "Ancestral Strain"
b1 <- read.delim("B1cov_10.txt", header = FALSE)
b1$strain <- "Biofilm Population 1"

all <- rbind(b1, anc)
colnames(all) <- c("chrom", "Position", "end", "ReadDepth", "Strain")

b <- ggplot(all, aes(x = Position, y=ReadDepth, color=Strain)) + 
  geom_line(aes(color=Strain)) + theme_classic() + 
  scale_y_continuous(limits=c(0,4500)) + 
  scale_x_continuous(breaks=c(0, 6537648), labels = c("0", "6,537,648 bp")) +
  scale_color_manual(values=c("dark gray", "#4DBBD5FF")) +
  theme(axis.title.x = element_blank(), 
        axis.text = element_text(size = 10, color = "black"),
        axis.line = element_line(size=0.25), 
        legend.position = c(0.85,0.9),
        legend.title = element_blank(),
        plot.margin = margin(10, 20, 10, 10)) + labs(title="PA14 Genome")

all_p <- all[which(all$Position > 4335000), ]
all_p <- all_p[which(all_p$Position < 4360000), ]

b_p <- ggplot() + 
  geom_line(data =  all_p, aes(x = Position, y=ReadDepth, color=Strain)) + 
  scale_y_continuous(limits=c(0,4500)) + 
  scale_x_continuous(breaks=c(4335000, 4360000), labels = c("4,335,000 bp", "4,360,000 bp")) +
  scale_color_manual(values=c("dark gray", "#4DBBD5FF")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 10, color = "black"), 
        axis.line = element_line(size=0.25), 
        plot.margin = margin(10, 20, 10, 10),
        axis.title.x = element_blank()) +
  labs(title = "Pf5") +
  geom_segment(aes(x = 4345126, y = 3100, xend =4345126, yend=3400), color = "purple", size=1) +
  geom_segment(aes(x = 4355790, y = 3100, xend =4355790, yend=3400), color = "purple", size=1) +
  geom_segment(aes(x = 4338302, y = 3100, xend =4338302, yend=3400), color = "blue", size=1) +
  geom_segment(aes(x = 4353244, y = 3100, xend =4353244, yend=3400), color = "blue", size=1) +
  geom_segment(aes(x = 4338990, y = 3100, xend =4338990, yend=3400), color = "green", size=1) +
  geom_segment(aes(x = 4346095, y = 3100, xend =4346095, yend=3400), color = "green", size=1) 

b1 <- ggarrange(b, b_p, labels = c("A", "B"),
                    ncol = 1, nrow = 2)
ggsave("BiofilmPopulation1.png", b1, width = 7, height = 3.5)


###########################################################################

### Clones 

anc <- read.delim("anccov_10.txt", header = FALSE)
anc$strain <- "Ancestral Strain"
lp <- read.delim("127cov_10.txt", header = FALSE)
lp$strain <- "lasR pf5r"
mp <- read.delim("128cov_10.txt", header = FALSE)
mp$strain <- "morA pf5r"
l <- read.delim("120cov_10.txt", header = FALSE)
l$strain <- "lasR"
m <- read.delim("126cov_10.txt", header = FALSE)
m$strain <- "morA"
gp <- read.delim("G1cov_10.txt", header = FALSE)
gp$strain <- "pf5r.1"
bp <- read.delim("pf5r.2cov_10.txt", header = FALSE)
bp$strain <- "pf5r.2"
w6 <- read.delim("1306washedcov_10.txt", header = FALSE)
w6$strain <- "lasR pf5r washed"
w7 <- read.delim("1307washedcov_10.txt", header = FALSE)
w7$strain <- "morA pf5r washed"

#   lasR
all <- rbind(lp, l)
colnames(all) <- c("chrom", "Position", "end", "ReadDepth", "Strain")

p <- ggplot(all, aes(x = Position, y=ReadDepth, color=Strain)) + 
  geom_line(aes(color=Strain)) + theme_classic() + 
  scale_y_continuous(limits=c(0,10000)) + 
  scale_x_continuous(breaks=c(0, 6537640), labels = c("0", "6,537,648 bp")) +
  scale_color_manual(values=c("dark gray", "#4DBBD5FF")) +
  theme(axis.title.x = element_blank(), 
        axis.text = element_text(size = 10, color = "black", family="Arial"),
        axis.line = element_line(size=0.25), 
        legend.position = c(0.85,0.9), 
        legend.title = element_blank(),
        plot.margin = margin(10, 20, 10, 10),
        legend.text = element_text(face = "italic")) + labs(title="PA14 Genome")

all_p <- all[which(all$Position > 4340000), ]
all_p <- all_p[which(all_p$Position < 4360000), ]

p_p <- ggplot() + 
  geom_line(data = all_p, aes(x = Position, y=ReadDepth, color=Strain)) + 
  scale_y_continuous(limits=c(0,10000)) + 
  scale_x_continuous(breaks=c(4340000, 4360000), labels = c("4,340,000 bp", "4,360,000 bp")) +
  scale_color_manual(values=c("dark gray", "#4DBBD5FF")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 10, color = "black", family = "Arial"), 
        axis.line = element_line(size=0.25), 
        plot.margin = margin(10, 25, 20, 10),
        axis.title.x = element_blank()) + labs(title="Pf5") 

all_t <- all[which(all$Position > 3170000), ]
all_t <- all_t[which(all_t$Position < 3190000), ]

t <- ggplot(all_t, aes(x = Position, y=ReadDepth, color=Strain)) + geom_line(aes(color=Strain)) +
  scale_x_continuous(breaks=c(3170000, 3190000), labels = c("3,170,000 bp", "3,190,000 bp")) +
  scale_y_continuous(limits=c(0,10000)) + 
  scale_color_manual(values=c("dark gray","#4DBBD5FF")) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"), 
        axis.title = element_text(size = 11, color="black"),
        axis.title.x = element_blank(),
        axis.line = element_line(size=0.25),
        legend.position = "none",
        plot.margin = margin(10, 25, 20, 10)) + labs(title="New Circularization")
figure <- ggarrange(p, p_p, t, labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)
ggsave("lasR.png", figure, width = 5, height = 7)

# morA
all <- rbind(mp, m)
colnames(all) <- c("chrom", "Position", "end", "ReadDepth", "Strain")

p <- ggplot(all, aes(x = Position, y=ReadDepth, color=Strain)) + 
  geom_line(aes(color=Strain)) + theme_classic() + 
  scale_y_continuous(limits=c(0,10000)) + 
  scale_x_continuous(breaks=c(0, 6537648), labels = c("0", "6,537,648 bp")) +
  scale_color_manual(values=c("dark green", "purple")) +
  theme(axis.title.x = element_blank(), 
        axis.text = element_text(size = 10, color = "black", family="Arial"),
        axis.line = element_line(size=0.25), 
        legend.position = c(0.85,0.9), 
        legend.title = element_blank(),
        plot.margin = margin(10, 20, 10, 10),
        legend.text = element_text(face = "italic")) + labs(title="PA14 Genome")

all_p <- all[which(all$Position > 4340000), ]
all_p <- all_p[which(all_p$Position < 4360000), ]

p_p <- ggplot() + 
  geom_line(data = all_p, aes(x = Position, y=ReadDepth, color=Strain)) + 
  scale_y_continuous(limits=c(0,10000)) + 
  scale_x_continuous(breaks=c(4340000, 4360000), labels = c("4,340,000 bp", "4,360,000 bp")) +
  scale_color_manual(values=c("dark green", "purple")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 10, color = "black", family = "Arial"), 
        axis.line = element_line(size=0.25), 
        plot.margin = margin(10, 25, 20, 10),
        axis.title.x = element_blank()) + labs(title="Pf5") 

all_t <- all[which(all$Position > 3170000), ]
all_t <- all_t[which(all_t$Position < 3190000), ]

t <- ggplot(all_t, aes(x = Position, y=ReadDepth, color=Strain)) + geom_line(aes(color=Strain)) +
  scale_x_continuous(breaks=c(3170000, 3190000), labels = c("3,170,000 bp", "3,190,000 bp")) +
  scale_y_continuous(limits=c(0,10000)) + 
  scale_color_manual(values=c("dark green","purple")) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"), 
        axis.title = element_text(size = 11, color="black"),
        axis.title.x = element_blank(),
        axis.line = element_line(size=0.25),
        legend.position = "none",
        plot.margin = margin(10, 25, 20, 10)) + labs(title="New Circularization")
figure <- ggarrange(p, p_p, t, labels = c("D", "E", "F"),
                    ncol = 1, nrow = 3)
ggsave("morA.png", figure, width = 5, height = 7)


#   pf5r.1
all <- gp
colnames(all) <- c("chrom", "Position", "end", "ReadDepth", "Strain")

p <- ggplot(all, aes(x = Position, y=ReadDepth, color=Strain)) + 
  geom_line(aes(color=Strain)) + theme_classic() + 
  scale_y_continuous(limits=c(0,10000)) + 
  scale_x_continuous(breaks=c(0, 6537648), labels = c("0", "6,537,648 bp")) +
  scale_color_manual(values=c("#4DBBD5FF", "#4DBBD5FF")) +
  theme(axis.title.x = element_blank(), 
        axis.text = element_text(size = 10, color = "black", family="Arial"),
        axis.line = element_line(size=0.25), 
        legend.position = c(0.85,0.9), 
        legend.title = element_blank(),
        plot.margin = margin(10, 20, 10, 10),
        legend.text = element_text(face = "italic")) + labs(title="PA14 Genome")

all_p <- all[which(all$Position > 4340000), ]
all_p <- all_p[which(all_p$Position < 4360000), ]

p_p <- ggplot() + 
  geom_line(data = all_p, aes(x = Position, y=ReadDepth, color=Strain)) + 
  scale_y_continuous(limits=c(0,10000)) + 
  scale_x_continuous(breaks=c(4340000, 4360000), labels = c("4,340,000 bp", "4,360,000 bp")) +
  scale_color_manual(values=c("#4DBBD5FF", "#4DBBD5FF")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 10, color = "black", family = "Arial"), 
        axis.line = element_line(size=0.25), 
        plot.margin = margin(10, 25, 20, 10),
        axis.title.x = element_blank()) + labs(title="Pf5") 

all_t <- all[which(all$Position > 3170000), ]
all_t <- all_t[which(all_t$Position < 3190000), ]

t <- ggplot(all_t, aes(x = Position, y=ReadDepth, color=Strain)) + geom_line(aes(color=Strain)) +
  scale_x_continuous(breaks=c(3170000, 3190000), labels = c("3,170,000 bp", "3,190,000 bp")) +
  scale_y_continuous(limits=c(0,10000)) + 
  scale_color_manual(values=c("dark gray","#4DBBD5FF")) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"), 
        axis.title = element_text(size = 11, color="black"),
        axis.title.x = element_blank(),
        axis.line = element_line(size=0.25),
        legend.position = "none",
        plot.margin = margin(10, 25, 20, 10)) + labs(title="Putative IS")
figure <- ggarrange(p, p_p, t, labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)
ggsave("g1.png", figure, width = 4, height = 7)

# pf5r.2
all <- bp
colnames(all) <- c("chrom", "Position", "end", "ReadDepth", "Strain")

p <- ggplot(all, aes(x = Position, y=ReadDepth, color=Strain)) + 
  geom_line(aes(color=Strain)) + theme_classic() + 
  scale_y_continuous(limits=c(0,10000)) + 
  scale_x_continuous(breaks=c(0, 6537648), labels = c("0", "6,537,648 bp")) +
  scale_color_manual(values=c("dark blue", "#4DBBD5FF")) +
  theme(axis.title.x = element_blank(), 
        axis.text = element_text(size = 10, color = "black", family="Arial"),
        axis.line = element_line(size=0.25), 
        legend.position = c(0.85,0.9), 
        legend.title = element_blank(),
        plot.margin = margin(10, 20, 10, 10),
        legend.text = element_text(face = "italic")) + labs(title="PA14 Genome")

all_p <- all[which(all$Position > 4340000), ]
all_p <- all_p[which(all_p$Position < 4360000), ]

p_p <- ggplot() + 
  geom_line(data = all_p, aes(x = Position, y=ReadDepth, color=Strain)) + 
  scale_y_continuous(limits=c(0,10000)) + 
  scale_x_continuous(breaks=c(4340000, 4360000), labels = c("4,340,000 bp", "4,360,000 bp")) +
  scale_color_manual(values=c("dark blue", "#4DBBD5FF")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 10, color = "black", family = "Arial"), 
        axis.line = element_line(size=0.25), 
        plot.margin = margin(10, 25, 20, 10),
        axis.title.x = element_blank()) + labs(title="Pf5") 

all_t <- all[which(all$Position > 3170000), ]
all_t <- all_t[which(all_t$Position < 3190000), ]

t <- ggplot(all_t, aes(x = Position, y=ReadDepth, color=Strain)) + geom_line(aes(color=Strain)) +
  scale_x_continuous(breaks=c(3170000, 3190000), labels = c("3,170,000 bp", "3,190,000 bp")) +
  scale_y_continuous(limits=c(0,10000)) + 
  scale_color_manual(values=c("dark gray","#4DBBD5FF")) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"), 
        axis.title = element_text(size = 11, color="black"),
        axis.title.x = element_blank(),
        axis.line = element_line(size=0.25),
        legend.position = "none",
        plot.margin = margin(10, 25, 20, 10)) + labs(title="Putative IS")
figure <- ggarrange(p, p_p, t, labels = c("C", "D", "C"),
                    ncol = 1, nrow = 3)
ggsave("b1.png", figure, width = 4, height = 7)




