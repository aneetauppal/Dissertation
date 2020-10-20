#Code used for PCA Analysis
#Aneeta Uppal and Additional credits to Colby Ford

library(dplyr)
library(tidyr)
library(readr)
library(dplyr)
library(mclust)
library(devtools)
library(ggbiplot)

#Uploading different data sets to run through PCA. data needs to be re-pivoted from original formats
#ALL EOUDB data
data<- read_csv("/Users/aneetauppal/Graduate_PhD/DissertationResearch/Knowledgebase/eoudb/EOUDB_forclustering.csv")
#ALL EOUDB data with family names
data<- read_csv("/Users/aneetauppal/Graduate_PhD/DissertationResearch/Knowledgebase/eoudb/EOUDBwithFamNames.csv")
#Anti-inflammtory oils of Sapindales family
data<- read_csv("/Users/aneetauppal/inflammoils_sapindales.csv")


## Pivot data to where each chemical is its own column with its respective relative percentage
data_pvt <- data %>%
  select(-Compound_ID) %>% 
  group_by(Family, Oil_Name, Oil_ID) %>%
  spread(key = Chemical_Name, value = Relative_Percentage) %>% 
  arrange(Oil_Name, Family, Oil_ID) %>% 
  replace(is.na(.), 0)

#Write transformed data to new file
#frankincense samples
write_csv(data_pvt, "/Users/aneetauppal/frankdata_trans.csv")
#"anti-inflammatory" oils of the sapindales family
write_csv(data_pvt, "/Users/aneetauppal/inflammoils_sapindales_pivot.csv")

#Upload data for using PCA analysis
data_pvt <- read_csv("/Users/aneetauppal/inflammoils_sapindales_pivot.csv")

#take all columns of just the chemical components and their relative percentages without the oil names/identifiers
eo_pca <- prcomp(data_pvt[7:169], center=TRUE, scale. = TRUE)

#Plot PCA with variables (chemicals)
ggbiplot(eo_pca, choices = 1:2, var.scale = 1,
         labels=NULL,
         ellipse = NULL,
         var.axes = TRUE, varname.size = 4.3) + ggtitle("PCA of chemical components of different EOs in Sapindales Order")+
  theme(plot.title = element_text(size=15, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")
  ) + xlim(-5,5) + ylim(-5,5)


#Plot PCA without variables, but with dots representing oil names, and colors & ellipses to signify
#groupings of the EOs based on their family names 
#choices represent PC axis that you want to plot
#labels = TRUE will plot PCA with the name of the oil on the graph instead of denoted dots 


ggbiplot(eo_pca, choices = 1:2, var.scale = 1,
         labels=NULL,
         ellipse = TRUE,
         var.axes = FALSE, groups = data_pvt$family, varname.size = 4.3) + ggtitle("PCA of chemical components of different oils in Sapindales Order")+
  theme(plot.title = element_text(size=15,hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")
  ) + xlim(-6,6) + ylim(-6,6)

