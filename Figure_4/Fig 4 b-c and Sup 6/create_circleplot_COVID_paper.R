# Libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(viridis)
circlepack_studies <- read_csv("C:/Users/Marie/Desktop/hCoCena_dev/thomas_model/circlepack_studies_noMEN_noFAT.csv")
total <- data.frame(from = c("total", "total", "total", "total", "total"), 
                    to = c("viral", "bacterial", "other", "viral and bacterial", "control"))
viral <- data.frame(from = c(rep("viral", nrow(circlepack_studies[circlepack_studies$Category == "viral",]))),
                    to = circlepack_studies[circlepack_studies$Category == "viral",]%>%dplyr::pull(., "StudyID"))
bacterial <- data.frame(from = c(rep("bacterial", nrow(circlepack_studies[circlepack_studies$Category == "bacterial",]))),
                        to = circlepack_studies[circlepack_studies$Category == "bacterial",]%>%dplyr::pull(., "StudyID"))
other <- data.frame(from = c(rep("other", nrow(circlepack_studies[circlepack_studies$Category == "other",]))),
                    to = circlepack_studies[circlepack_studies$Category == "other",]%>%dplyr::pull(., "StudyID"))
viral_and_bacterial<- data.frame(from = c(rep("viral and bacterial", nrow(circlepack_studies[circlepack_studies$Category == "viral and bacterial",]))),
                                 to = circlepack_studies[circlepack_studies$Category == "viral and bacterial",]%>%dplyr::pull(., "StudyID"))
control <- data.frame(from = c(rep("control", nrow(circlepack_studies[circlepack_studies$Category == "control",]))),
                      to = circlepack_studies[circlepack_studies$Category == "control",]%>%dplyr::pull(., "StudyID"))

edges <- rbind(total, viral, bacterial, other, viral_and_bacterial, control)
totalv <- data.frame(StudyID = "total", 
                     size = sum(circlepack_studies$size),
                     DiseaseCondition = "Disease space",
                     Category = "total",
                     name = "total")
viralv <- data.frame(StudyID = "viral", 
                     size = sum(circlepack_studies[circlepack_studies$Category == "viral",] %>% 
                                  dplyr::pull(., size)),
                     DiseaseCondition = "viral space",
                     Category = "viral",
                     name = "viral")
bacterialv <- data.frame(StudyID = "bacterial", 
                         size = sum(circlepack_studies[circlepack_studies$Category == "bacterial",] %>% 
                                      dplyr::pull(., size)),
                         DiseaseCondition = "bacterial space",
                         Category = "bacterial",
                         name = "bacterial")
otherv <- data.frame(StudyID = "other", 
                     size = sum(circlepack_studies[circlepack_studies$Category == "other",] %>% 
                                  dplyr::pull(., size)),
                     DiseaseCondition = "other space",
                     Category = "other",
                     name = "other")
vabv <- data.frame(StudyID = "viral and bacterial", 
                   size = sum(circlepack_studies[circlepack_studies$Category == "viral and bacterial",] %>% 
                                dplyr::pull(., size)),
                   DiseaseCondition = "viral and bacterial space",
                   Category = "viral and bacterial",
                   name = "viral and bacterial")
controlv <- data.frame(StudyID = "control", 
                       size = sum(circlepack_studies[circlepack_studies$Category == "control",] %>% 
                                    dplyr::pull(., size)),
                       DiseaseCondition = "control",
                       Category = "control",
                       name = "control")

circlepack_studies$name <- circlepack_studies$StudyID
vertices <- rbind(circlepack_studies, totalv, viralv, bacterialv, otherv, vabv, controlv)
# Rebuild the graph object
mygraph <- graph_from_data_frame(edges, vertices = vertices)

g <-ggraph(mygraph, layout = "circlepack", weight = size) +
  geom_node_circle(aes(fill = DiseaseCondition)) +
  scale_fill_manual(values=c("Disease space" = NA, "viral space" = NA, 
                             "bacterial space" = NA, "other space" = NA, "viral and bacterial space"= NA,
                             "control" = NA,
                             "Autoimmune disease" = "skyblue4", "Control" = "#00A251", 
                             "Autoimmune disease" = "skyblue4", "Chikungunya" = "darkgoldenrod2", "COVID-19"= "#522685",
                             "Ebola vaccine" = "lightpink2", "NLRC4-MAS and NOMID" = "turquoise3", "Rheumatoid arthritis" = "royalblue1", 
                             "Sepsis and SIRS" = "skyblue1", "Tuberculosis" = "mediumorchid3", "Tuberculosis"= "mediumorchid3",
                             "Tuberculosis" = "mediumorchid3", "Zika" = "darkgoldenrod", "Control" = "#00A251", 
                             "Tuberculosis + HIV" = "sandybrown", "Tuberculosis" = "mediumorchid3", "Tuberculosis"= "mediumorchid3"
  )) +
  
  geom_text(aes(x = x, y = y, label = name)) +
  coord_fixed()+
  theme_void()
g
# Cairo::CairoPDF(file ="C:/Users/Marie/Desktop/hCoCena_dev/thomas_model/Supp_fig4_sup_circles.pdf", 
#                 width = 10, height = 10)
# plot(g)
# dev.off()
