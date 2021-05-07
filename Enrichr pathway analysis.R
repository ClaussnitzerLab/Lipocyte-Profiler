#####Pathway Enrichr analysis####
#####set directory####
setwd("/Users/slaber/Desktop/")

#####load gene list for pathway enrichemnet####
##eg. connections between LP feature and gene
LP<-subset(pie_3, pie_3$variable == "Cells_Granularity_10_Mito")

library(enrichR)
#####define pathway enrichment gene list you want to use###
enrich<-enrichR::enrichr(as.character(as.character(LP$gene)), c("KEGG_2019_Human",
                                                                "BioPlanet_2019",
                                                                "WikiPathways_2019_Human",
                                                                "GO_Biological_Process_2018"))


#####check top rows of enriched pathways###
head(enrich[["KEGG_2019_Human"]])
head(enrich[["BioPlanet_2019"]])
head(enrich[["WikiPathways_2019_Human"]])
head(enrich[["GO_Biological_Process_2018"]])

######select on list of pathway###
enrichr_data<-enrich[["WikiPathways_2019_Human"]]
#####save analysis###
write.table(enrichr_data, file="Mito_Granularity enrichr data.csv", sep=",", col.names = NA, qmethod = "double")


######plot enrichr analysis####
#####define significance level####
library(dplyr)
plot_input <- enrichr_data %>% dplyr::mutate(label = ifelse(enrichr_data$Adjusted.P.value <= 0.05, "<5%FDR", ">5%FDR"))

####plot####
library(ggplot2)

plot<-ggplot(plot,aes(x=Odds.Ratio, y=-log10(Adjusted.P.value), size=-log10(Adjusted.P.value), fill=label)) +
            geom_point(alpha=0.6, color="black", pch=21) +
            geom_text(aes(label = ifelse(Adjusted.P.value <= 0.05, as.character(Term), "") ,size = -log10(Adjusted.P.value)),
                      hjust=0, vjust=0,color = "black", check_overlap = F)+
            scale_size(range = c(.01, 4), name="-log(Adjusted.P.value)")+
            scale_fill_manual(name = "Adjusted.P.value", values= c("<5%FDR"= "#4DAC26",
                                                                    ">5%FDR"= "#B9C1C7" ))+
            theme_bw()


pdf("enrichr plot Cells_Granularity_Mito.pdf")
plot_input
dev.off()
