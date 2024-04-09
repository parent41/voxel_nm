
library(ggplot2)
library(dplyr)
library(viridis)

df = read.csv("./PubMed_Timeline_Results_by_Year.csv")
df$Year = rownames(df)
df$Count = df[,1]

df = df[seq(2, nrow(df)), c(2,3)]
df$Year = as.numeric(df$Year)
df$Count = as.numeric(df$Count)

ggplot(df, aes(x=Year, y=Count, fill=Count)) +
    geom_bar(stat="identity") + 
    scale_x_continuous(name="Year", limits=c(1999, 2025), breaks=c(2000, 2005, 2010, 2015, 2020, 2024)) + 
    scale_fill_viridis_c(option = "mako") +
    labs(title="Pubmed results for \"MRI\" and \"Microstructure\"") +
    theme_classic() + 
    theme(legend.position = "none", text=element_text(size=20), plot.title = element_text(hjust = 0.5))
ggsave("pubmed_micro.png", width=10, height=5)

