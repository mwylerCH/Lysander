# PCA dist matrix

library(data.table)
library(tidyverse)
library(Rtsne)
library(ggbiplot)
library(umap) 

## PCA table ----------------------

# prepare classification
classification <- fread('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/allClassified.txt',
                        data.table = F, header = F, fill = T, sep = "\t") %>% 
  na.omit

# add american
#usa <- fread('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/USAClassified.txt',
#          data.table = F, header = F, fill = T) %>% 
#  na.omit
#classification <- rbind(classification, usa)

# Issue 43
classification <- classification %>% 
  mutate(V3 = ifelse(V6 == 43 & V8 == 43 & V11 == 1, 'BB', V3))

# prepare future table for PCA
PCAtable <- classification[, c(1,3)]

# make NA instead of ''
PCAtable$V3[PCAtable$V3 == ''] <- NA

# read in segment distance OLD
# for (NR in 1:8){
#   nome <- paste0('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/Segments_', NR, '.DistMatrix')
#   pezzo <- fread(nome,
#                  data.table = F)
#   # clean names
#   pezzo$V1 <- gsub('^[^|]+[[:digit:]]\\|', '', pezzo$V1)
#   pezzo$V1 <- gsub('\\|.+$', '', pezzo$V1)
#   pezzo$V1 <- gsub('/', '_', pezzo$V1)
#   
#   
#   # order based on classification file
#   ordine <- pezzo[match(classification$V1, pezzo$V1),]
#   
#   # remove duplicated columns
#   ordine <- ordine[, !duplicated(colnames(ordine))]
# 
#   # remove ID
#   ordine$V1 <- NULL
#   
#   PCAtable <- cbind(PCAtable, ordine)
# }

# representative sequences
nomi <- fread('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/ReferenceNames.txt',
              header = F, sep = '\t')
nomi[is.na(nomi$V2),2] <- 'NA'

# read in segment distance
for (NR in 1:8){
  nameSeg <- c('PB2', 'PB1', 'PA', 'HA','NP', 'NA','MP','NS')
  nome <- paste0('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/Segment_', nameSeg[NR], '_2.DistMatrix')
  pezzo <- fread(nome,
                 data.table = F)
  
  # filter using representative sequences
  keepNames <- c('V1', 'V3', as.vector(unlist(nomi[nomi$V2 == nameSeg[NR], 1])))
  pezzo <- pezzo[, colnames(pezzo) %in% keepNames]
  
  # clean names
  pezzo$V1 <- gsub('^[^|]+[[:digit:]]\\|', '', pezzo$V1)
  pezzo$V1 <- gsub('\\|.+$', '', pezzo$V1)
  pezzo$V1 <- gsub('/', '_', pezzo$V1)
  
  # order based on classification file
  ordine <- pezzo[match(classification$V1, pezzo$V1),]
  
  # remove duplicated columns
  ordine <- ordine[, !duplicated(colnames(ordine))]
  
  # remove ID
  ordine$V1 <- NULL
  
  PCAtable <- cbind(PCAtable, ordine)
  print(dim(PCAtable))
}

PCAbackup <- PCAtable
# PCAtable <- PCAbackup

# make NA instead of empty
#sapply(PCAtable$V3, function(x){ifelse(x=='', NA, x)}) 

# random sample columns
#randomCols <- floor(runif(4000, min=3, max=ncol(PCAtable)))
#PCAtable <- PCAtable[, c(1:2, randomCols)]

# only keep References
#referenceNames <- colnames(PCAtable)[grepl('REF', colnames(PCAtable))]
#PCAtable <- PCAtable[, c('V1','V3', referenceNames)]

# remove empty lines (most are references, maybe some mismatch)
#PCAtable <- PCAtable[PCAtable$V != 1, ]

# keep only major genotypes
genoCount <- table(PCAtable$V3)
importantGeno <- names(genoCount[genoCount > 50])
#PCAtable$V3 <- ifelse(PCAtable$V3 %in% importantGeno, PCAtable$V3, '')

# one is NA (?)
PCAtable <- PCAtable[!is.na(PCAtable$V1), ]


# remove duplicated rows
noId <- PCAtable[, 3:ncol(PCAtable)]
duplicatedRows <- duplicated(noId)
PCAtable <- PCAtable[!duplicatedRows,]

# find identical columns
PCAtable <- PCAtable[!duplicated(as.list(PCAtable))]
#fwrite(PCAtable[, 3:ncol(PCAtable)], file = '//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/forPython.csv', sep = ',', quote = F)

# check reference BB (also wrong)
#PCAtable %>% 
#  filter(grepl('22P015977', V1))

# problem with A_turkey_Germany-BB_AI00868_2021
PCAtable <- PCAtable[!is.na(PCAtable$`A_/_H5N1|A/American_Crow/BC/AIVPHL-1384/2023|PB2|1|EPI_ISL_18665536`), ]

# split up for lysander
tab1 <- PCAtable[1:3000,]
tab2 <- PCAtable[3001:6000,]
tab3 <- PCAtable[6001:nrow(PCAtable),]
fwrite(tab1, file = '//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/forLysander_Tab1.csv')
fwrite(tab2, file = '//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/forLysander_Tab2.csv')
fwrite(tab3, file = '//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/forLysander_Tab3.csv')

## American Classification ----------------------------

madeUS <- fread('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/US_classification.txt',
                data.table = F, header = T)
madeUS$`File Name` <- gsub('.fa', '', madeUS$`File Name`)
# make NA all not claffisied
madeUS[grepl('Not', madeUS$Genotype), 2] <- NA

# order as PCAtable
USgeno <- madeUS[match(PCAtable$V1, madeUS$`File Name`),'Genotype']



## Geo vector ---------------------------------------

## Geography vector
library(r2country)
continentReference <- continentOf %>% unlist %>% as.data.frame()
continentReference$country <- rownames(continentReference)
colnames(continentReference)[1] <- 'continent'
continentReference$country <- gsub(' ','_', continentReference$country)

# more then 200 samples from russian regions, but also many cities
russianProvince <- grep("russia,",names(timeIn), value = T)
russianProvince <- gsub('russia, ', '', russianProvince)
russianProvince <- c(russianProvince, 'rostov', 'ossetia', 'krasnodar', 'leningrad', 'saratov', 'kostroma',
                     'voronezh', 'tyumen', 'tatarstan', 'kursk', 'kurgan', 'penza', 'tyva')
# us states
USstates <- grep("usa,",names(timeIn), value = T)
USstates <- gsub('usa, ', '', USstates)
USstates <- gsub(', .*', '', USstates)
USstates <- gsub(' ', '_', USstates)

# get ID of samples (lower capital)
geoVec <- tolower(PCAtable$V1)
geoVec2 <- tolower(PCAtable$V1)

# UK vs single states
geoVec <- gsub('england', 'united_kingdom', geoVec)
geoVec <- gsub('scotland', 'united_kingdom', geoVec)
geoVec <- gsub('wales', 'united_kingdom', geoVec)
geoVec <- gsub('nl-', 'netherlands', geoVec)
geoVec <- ifelse(grepl(paste0(russianProvince, collapse = '|'), geoVec), 'russia', geoVec)
geoVec <- ifelse(grepl(paste0(USstates, collapse = '|'), geoVec), 'united_states', geoVec)

# convert to continent
for (ROW in 1:nrow(continentReference)){
  continente <- continentReference$continent[ROW]
  paese <- continentReference$country[ROW]
  # pick IDs of interest
  geoVec[grepl(paese, geoVec)] <- continente
  
}
table(geoVec)

# remove bad samples with no geo position
geoVec[!geoVec %in% continentReference$continent] <- NA


## Year vector -----------------------------------------

yearVec <- PCAtable$V1

# extract year
yearVec <- gsub(".*_([[:digit:]]{4}$)", "\\1", yearVec)
yearVec <- gsub(".*_([[:digit:]]{4})-[[:digit:]]{2}-[[:digit:]]{2}$", "\\1", yearVec)
yearVec <- gsub(".*-([[:digit:]]{4}$)", "\\1", yearVec)

# remove remaining
yearVec[yearVec < 1900] <- NA
yearVec[grepl('[[:alpha:]]', yearVec)] <- NA


## try RTsne -----------------------------------------

# prepare matrix
matrice <- as.matrix(PCAtable[, 3:ncol(PCAtable)])
X <- normalize_input(matrice)

set.seed(12345)
tsneRes <- Rtsne(X,
                 # perplexity should been between 5 and 50
                 perplexity = 50,
                 num_threads = 20)

# alternative UMAP
#iris.umap = umap(X, n_components = 2, random_state = 15) 
#layout <- iris.umap[["layout"]] 
#layout <- data.frame(layout) 
#final <- cbind(layout, PCAtable$V3) 
#colnames(final) <- c('x','y','col')

#final %>% 
#  mutate(col = ifelse(col %in% importantGeno, col, NA)) %>% 
#  ggplot() + 
#  geom_point(aes(x=x, y=y, color=col))+
#  theme_bw()

## plot -----------------------------------------

# plot most common
importantGeno <- names(genoCount[genoCount > 50])

tsne_plot <- data.frame(x = tsneRes$Y[,1], y = tsneRes$Y[,2], col = PCAtable$V3, 
                        continent=geoVec, year = yearVec, USA = USgeno)
tsne_plot %>% 
  #filter(col != '12') %>% dim()
  mutate(col = ifelse(col %in% importantGeno, col, NA)) %>% 
  #filter(col == 'brutti') %>% 
  #filter(col %in% importantGeno) %>% 
  #filter(col != '') %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color=col))+
  theme_bw()+
  labs(title = 'Major European')
#ggsave('MWy/tsne_AIV_subsampled.pdf', width = 10, height = 10)

# all 
tsne_plot %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color=col))+
  theme_bw()+
  labs(title = 'European')+
  theme(legend.position = '')
# ggsave('MWy/tsne_AIV_subsampled_ALLgeno.pdf', width = 10, height = 10)

# plot by continent
tsne_plot %>% 
  #na.omit() %>% 
  #filter(col %in% importantGeno) %>%
  #filter(col != '') %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color=continent), alpha = 0.8)+
  theme_bw()+
  scale_color_manual(breaks = c("Africa", "Asia", "Europe", "Europe & Asia", "North America", "Ocoeania"),
                     values=c("purple", "blue", "red", "green", "black", "yellow"))+
  labs(title = 'Continent')
# ggsave('MWy/tsne_AIV_subsampled_continent.pdf', width = 10, height = 10)

# year
tsne_plot %>% 
  mutate(year = ifelse(year >= 2020, year, NA)) %>% 
  #na.omit() %>% 
  #filter(col %in% importantGeno) %>%
  #filter(col != '') %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color=as.numeric(year)), alpha = 0.8)+
  theme_bw()+
  labs(title = 'Year')


# american Genotyping
tsne_plot %>% 
  #filter(col != '12') %>% dim()
  #filter(col == 'brutti') %>% 
  #filter(col %in% importantGeno) %>% 
  #filter(col != '') %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color=USA))+
  theme_bw()+
  labs(title = 'US Genotyping')+
  theme(legend.position = '')
# ggsave('MWy/tsne_AIV_subsampled_USA.pdf', width = 10, height = 10)


## cluster plotting ------------------------------------------
library(dbscan)

tsne2Dim <- tsne_plot[, 1:2]
kNNdistplot(tsne2Dim, k = 1)
abline(h=0.3, col = "red", lty=2)

dbscan_res <- dbscan(tsne2Dim, eps = 3, minPts = 0.03)
dbscan_res
dbscan_res$cluster

# merge table
tsne2Dim$cluster <- as.factor(dbscan_res$cluster)
# add also main table
tsne_plot$newCluster <- as.factor(dbscan_res$cluster)


tsne2Dim %>% 
  #filter(col != '12') %>% dim()
  #filter(col == 'brutti') %>% 
  #filter(col %in% importantGeno) %>% 
  #filter(col != '') %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color= cluster))+
  theme_bw()+
  labs(title = 'Statistical Genotyping')+
  theme(legend.position = '')
# ggsave('MWy/tsne_AIV_subsampled_HomeMade.pdf', width = 10, height = 10)


## Statistics ----------------------------------

# count samples = 9341
nrow(tsne_plot)
# count sequences = 74728
nrow(tsne_plot)*8

# check how many exclusive clades (only the new method) = 24
onlyNew <- unique(tsne_plot[is.na(tsne_plot$col) & is.na(tsne_plot$USA), 'newCluster'])
for (CL in onlyNew){
  # table subset
  subTable <- tsne_plot[tsne_plot$newCluster == CL, ]
  NRsubTable <- nrow(subTable)
  # count unclassified
  NRunclassified <- nrow(unique(tsne_plot[is.na(tsne_plot$col) & is.na(tsne_plot$USA) & tsne_plot$newCluster == CL, ]))
  # remove if less then 5% not covered
  amountUnclassified <- NRunclassified/NRsubTable*100
  if (amountUnclassified < 5){
    onlyNew <- onlyNew[onlyNew != CL]
  }
}
length(onlyNew )

# check continents
simpleGeo <- tsne_plot
simpleGeo$continent <- ifelse(!simpleGeo$continent %in% c("Europe", "North America"), 'varia', simpleGeo$continent)
# North America=3425, Europe=4608, varia=1308
table(simpleGeo$continent)

# how much NA are classified by USDA? 74.83%
simpleGeo[simpleGeo$continent == 'North America',] %>% filter(!is.na(USA)) %>% nrow()/3425*100
# how much NA are classified by EURL? 16.32%
simpleGeo[simpleGeo$continent == 'North America',] %>% filter(!is.na(col)) %>% nrow()/3425*100

# how much Euro are classified by USDA? 4.34%
simpleGeo[simpleGeo$continent == 'Europe',] %>% filter(!is.na(USA)) %>% nrow()/4608*100
# how much Euro are classified by EURL? 78.23%
simpleGeo[simpleGeo$continent == 'Europe',] %>% filter(!is.na(col)) %>% nrow()/4608*100

# how much varia are classified by USDA? 3.89%
simpleGeo[simpleGeo$continent == 'varia',] %>% filter(!is.na(USA)) %>% nrow()/1308*100
# how much Euro are classified by EURL? 36.85%
simpleGeo[simpleGeo$continent == 'varia',] %>% filter(!is.na(col)) %>% nrow()/1308*100


### MooFlu
mammals <- tsne_plot
mammals$header <- PCAtable$V1
mammals$MooFlu <- NA
mammals[grepl('goat|dairy_cow|feline', mammals$header), 'MooFlu'] <- 'MooFloo'



mammals %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color= MooFlu))+
  theme_bw()+
  labs(title = 'MooFlu')+
  theme(legend.position = '')
# ggsave('MWy/tsne_AIV_subsampled_HomeMade.pdf', width = 10, height = 10)



## US Cow Pandemic -----------------------------

MammalVec <- PCAtable$V1
USinteresting <- MammalVec[grepl('goat|dairy_cow|feline', MammalVec)]
MammalVec[!grepl('goat|dairy_cow|feline', MammalVec)] <- 'others'
MammalVec[grepl('goat|dairy_cow|feline', MammalVec)] <- 'USmammals'

# all same references
classificationUSmammal <- fread('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/pulitiClassified.txt',
                                data.table = F, header = F, fill = T, sep = '\t') %>%
  filter(V1 %in% USinteresting)

# same references like who? at least 2022
codeUSA <- distinct(classificationUSmammal[4:ncol(classificationUSmammal)])
yearCheck <- fread('~/temp_AIVpca/pulitiClassified.txt',
                   data.table = F, header = F, fill = T, sep = '\t') %>%
  filter(V4 == 12, V5 == 33, V6 == 1, V7 == 20,
         V8 == 21, V9 == 1, V10 == 20, V11 == 28) %>%
  mutate(anno = gsub('.*_([[:digit:]]{4})$', '\\1', V1))

# plot
cbind(tsne_plot, MammalVec) %>%
  ggplot() +
  geom_point(aes(x=x, y=y, color=MammalVec))+
  theme_bw()
#ggsave('~/temp_AIVpca/tsne_AIV_subsampled_Mammals.pdf', width = 10, height = 10)

cbind(tsne_plot, MammalVec) %>%
  filter(MammalVec == 'USmammals') %>%
  ggplot() +
  geom_point(aes(x=x, y=y, color=MammalVec))+
  theme_bw()
#ggsave('~/temp_AIVpca/tsne_AIV_subsampled_onlyMammals.pdf', width = 10, height = 10)


# check which are making groups
mammalPositions <- cbind(tsne_plot, MammalVec, PCAtable$V1) %>%
  filter(MammalVec == 'USmammals')


## calculate PCA ----------------------------
pcaRes <- prcomp(PCAtable[, 3:ncol(PCAtable)])

ggbiplot(pcaRes,
         groups = PCAtable$V3, var.axes = F)

# column centering
#centeredDF <- sweep(PCAtable[, 3:ncol(PCAtable)], 2, colMeans(PCAtable[, 3:ncol(PCAtable)]), "-")
# svd
#sv <- svd(centeredDF, 2,0)

# plot results
#plot(sv$u[, 1], sv$u[, 2], col = PCAtable$V3, main = "SVD", xlab = "U1", ylab = "U2")
#plot(sv$u[, 1], sv$u[, 2], main = "SVD", xlab = "U1", ylab = "U2")



## LDA --------------------
library(MASS)

linear <- lda(V3 ~., PCAtable[PCAtable$V3 != 'BB', 2:ncol(PCAtable)])
#plot(linear)
ggbiplot(linear, groups = PCAtable$V3, var.axes = F)
#ggsave('MWy/LDA_AIV_subsampled.pdf')
ggbiplot(linear, groups = PCAtable[PCAtable$V3 != 'BB', 'V3'], var.axes = F)


## Rdimtools ----------------------------------

library(Rdimtools)

# prepare data
lab <- as.factor(as.character(PCAtable$V3))
X <- as.matrix(PCAtable[, 3:ncol(PCAtable)])

# check estimated dimensions
est.Ustat(X)


# dimension reductions
resPCA <- do.pca(X, ndim = 2)
cbind(resPCA$Y, as.character(lab)) %>% as.data.frame() %>% 
  ggplot(aes(x=as.numeric(V1), y= as.numeric(V2), color = V3)) + 
  geom_point()




risultato <- do.lscore(X, ndim = 2)
cbind(risultato$Y, as.character(lab)) %>% as.data.frame() %>% 
  filter(V3 %in% importantGeno) %>% 
  filter(V3 != '') %>% 
  ggplot(aes(x=as.numeric(V1), y= as.numeric(V2), color = V3)) + 
  geom_point()


## RF classifier -------------------------------------------------------


# split data
set.seed(67)
train.Data <- sample(nrow(PCAtable), 2/3 * nrow(PCAtable))
Data_train <- PCAtable[train.Data, -1] 
Data_test <- PCAtable[-train.Data, -1]
colnames(Data_train) <- gsub('\\|([[:digit:]])\\|','seg\\1', colnames(Data_train))
colnames(Data_test) <- gsub('\\|([[:digit:]])\\|','seg\\1', colnames(Data_test))


importance <- 'permutation'
mtry <- round(ncol(Data_train)/3)
num.trees <- 500
replace <- F
min.node.size <- 5

set.seed(67)
RFmodel_out_w <- ranger(as.factor(V3)~ ., data=Data_train,
                        importance = importance,
                        mtry = mtry,
                        num.trees = num.trees,
                        replace = replace,
                        min.node.size = min.node.size,
                        num.threads = 7)

RFmodel_out_w$r.squared
RFmodel_out_w$prediction.error
head(data.frame(V1=sort(RFmodel_out_w$variable.importance, decreasing = TRUE)),20)

# run prediction
predictedTest <- predict(RFmodel_out_w, data = Data_test[,-1])

# check how many correct
numeroRow <- predictedTest$predictions %>% length()

comparisonRF <- as.data.frame(matrix(nrow = numeroRow, ncol = 2))
colnames(comparisonRF) <- c('real', 'predicted')
comparisonRF$real <- Data_test$V3
comparisonRF$predicted <- predictedTest$predictions

# how many wrong? 1.5%
nrow(comparisonRF[comparisonRF$real != comparisonRF$predicted, ])/nrow(comparisonRF)*100

# which one? only rare
comparisonRF[comparisonRF$real != comparisonRF$predicted, ]



## US Cow Pandemic -----------------------------

MammalVec <- PCAtable$V1
USinteresting <- MammalVec[grepl('goat|dairy_cow|feline', MammalVec)]
MammalVec[!grepl('goat|dairy_cow|feline', MammalVec)] <- 'others'
MammalVec[grepl('goat|dairy_cow|feline', MammalVec)] <- 'USmammals'

# all same references
classificationUSmammal <- fread('~/temp_AIVpca/pulitiClassified.txt',
                                data.table = F, header = F, fill = T, sep = '\t') %>%
  filter(V1 %in% USinteresting)

# same references like who? at least 2022
codeUSA <- distinct(classificationUSmammal[4:ncol(classificationUSmammal)])
yearCheck <- fread('~/temp_AIVpca/pulitiClassified.txt',
                   data.table = F, header = F, fill = T, sep = '\t') %>%
  filter(V4 == 12, V5 == 33, V6 == 1, V7 == 20,
         V8 == 21, V9 == 1, V10 == 20, V11 == 28) %>%
  mutate(anno = gsub('.*_([[:digit:]]{4})$', '\\1', V1))

# plot
cbind(tsne_plot, MammalVec) %>%
  ggplot() +
  geom_point(aes(x=x, y=y, color=MammalVec))+
  theme_bw()
#ggsave('~/temp_AIVpca/tsne_AIV_subsampled_Mammals.pdf', width = 10, height = 10)

cbind(tsne_plot, MammalVec) %>%
  filter(MammalVec == 'USmammals') %>%
  ggplot() +
  geom_point(aes(x=x, y=y, color=MammalVec))+
  theme_bw()
#ggsave('~/temp_AIVpca/tsne_AIV_subsampled_onlyMammals.pdf', width = 10, height = 10)


# check which are making groups
mammalPositions <- cbind(tsne_plot, MammalVec, PCAtable$V1) %>%
  filter(MammalVec == 'USmammals')
