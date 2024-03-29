###### Biogeographical patterns of plant phylogenetic diversity at the European scale ######

# Data loading
load("/Users/sofiaprandelli/tesi/sPlotOpen.RData")

#check species --> grep / stringsplit / tidyverse
DT2.oa$Species <- gsub(DT2.oa$Species,pattern = ' ',replacement = '_')
species_level <- DT2.oa[grepl("_", DT2.oa$Species),]
NOspecies_level <- DT2.oa[!grepl("_", DT2.oa$Species),] #osservazioni che non possono essere considerate
View(table(NOspecies_level$Species))

#d$d.Species <- gsub(d$d.Species,pattern = ' ',replacement = '_')

#Unione delle 2 tabelle necessarie per le analisi e selezione delle variabili
d <- merge(species_level, header.oa, by='PlotObservationID')
d <- data.frame(d$PlotObservationID, d$Species, d$Original_abundance, d$Abundance_scale, d$Continent, d$Country, d$Latitude, d$Longitude, d$Relative_cover)

#Filtriamo per il continente Europa e teniamo solo le specie con le abbondanze x_IC (n° individui per plot) e pa
d <- d[d$d.Continent=='Europe',]
d <- d[d$d.Abundance_scale != 'x_SC',]
names(d)
unique(d$d.Original_abundance)

# set wd
library(bigreadr)
library(tidyverse)
library(CoordinateCleaner)
library(gridExtra)
library(countrycode)
library(data.table)
library(sf)
library(rnaturalearth)
library(raster)
library(RStoolbox)
library(rasterVis)
library(viridis)
library(RColorBrewer)
library(phyloregion)
library(ggplot2)
library(viridis)

#### CREATION OF GRID 0.5 RESOLUTION ####
# cambio la risoluzione della maglia da 100 km a 50 km 
world <- ne_countries(scale = "medium", returnclass = "sf")
r<-raster::extent(c(xmin=-30, xmax=54, ymin=25, ymax=74)) #extent tagliato sull'Europa
p <- as(r, "SpatialPolygons")
m <- sp::SpatialPolygonsDataFrame(p, data.frame(sp = "x"))
m2 <- fishnet(mask = m, res = 0.5) #risoluzione a 0.5 gradi = 50 km
crs(m2) <- crs(world) #set the Coordinate reference system
m2$id=1:nrow(m2)
m2=m2[,-1] #tolgo l'oggetto "grids" da m2

#### PRESENCE-ABSENCE MATRIX SPP1PA_2 ####
coordinates(d)= ~d.Longitude+d.Latitude #trasforma oggetto in spatial dataframe df 
crs(d)=crs(m1)      #set the Coordinate reference system --> coordinate in WGS84
ov2 <- over(d[,2], m2) #intersect griglia-occorrenze, double-check che colonna selezionata sia la specie
y1_2 <- cbind(as.data.frame(d[,c(1, 2, 3, 4, 5, 6)]), ov2) #add species names, genus, family, countryCode 
#add d.PlantObservationID, d.Species, d.Original_Abundance, d.Abundance_Scale, d.Continent, d.Country ??
y1_2 <- y1_2[complete.cases(y1_2), ] #adding grids value
spp1_2 <- data.frame(as.matrix(long2sparse(y1_2, grids = "id", species = "d.Species")))[,-1] #spp1 matrice di abbondanza delle specie
spp1pa_2 <- vegan::decostand(spp1_2,method="pa") #trasformo in matrice Presence-Absence
rm(spp1_2)
rm(y1_2)
rm(ov2)

#### STANDARDIZATION of plant names using TPL ####
nomiTab=c("Taxon","Genus","Hybrid.marker","Species","Abbrev", "Infraspecific.rank",
          "Infraspecific","Authority", "ID", "Plant.Name.Index", "TPL.version",
          "Taxonomic.status","Family","New.Genus","New.Hybrid.marker","New.Species",
          "New.Infraspecific.rank", "New.Infraspecific",  "New.Authority","New.ID",
          "New.Taxonomic.status", "Typo", "WFormat", "Higher.level","Date")

namesSpPa2 <- gsub(colnames(spp1pa_2),pattern = '_',replacement = ' ')
nrowSpPa2 <- length(namesSpPa2)
standNamesSpPa2 <- as.data.frame(matrix(NA, nrow = nrowSpPa2, ncol=25,
                                        dimnames=list(namesSpPa2, nomiTab)))
batchSize=700

#install.packages("Taxonstand") #requires Taxonstand (for the TPL function)
library(Taxonstand)

for(i in 1:(ceiling(nrowSpPa2/batchSize))){
  cat(paste0("\nBatch ", i," out of ", (ceiling(nrowSpPa2/batchSize)), "\n"))
  rowsSelect <- ((i-1) * batchSize + 1):min((i*batchSize), nrowSpPa2)
  taxoAux2 <- as.data.frame(TPL(namesSpPa2[rowsSelect],infra =TRUE,corr=TRUE,
                                drop.lower.level=TRUE,repeats = 10))
  standNamesSpPa2[rowsSelect,] <- taxoAux2
}

saveRDS(standNamesSpPa2,'/Users/sofiaprandelli/tesi/standnamesSpPa2.rds')

taxoPlants2<-data.frame("Species"=paste(standNamesSpPa2$New.Genus,standNamesSpPa2$New.Species),
                        "Genus"=standNamesSpPa2$Genus,"Family"=standNamesSpPa2$Family)

taxoPlants2$Species <- gsub(taxoPlants2$Species,pattern = ' ',replacement = '_')
taxoPlants2$Species <- gsub(taxoPlants2$Species,pattern = '-',replacement = '_')

#### Creation of the PHYLOGENETIC TREE ####
devtools::install_github("jinyizju/V.PhyloMaker", force=TRUE)
library(V.PhyloMaker)
phylogenyPlants2 <- phylo.maker(taxoPlants2, scenarios=c("S1","S3")) #taxoPlants tassonomia di riferimento  
myTree2 <- phylogenyPlants2$scenario.3

##### PD CALCULATION #####
install.packages("picante")
library(picante)
dim(spp1pa_2)[2]==length(myTree2$tip.label) #check che spp nella matrice di comunità siano le stesse che nell'albero
PDobs_2 <- pd(spp1pa_2, myTree2, include.root=TRUE)
write.csv(PDobs_2, '/Users/sofiaprandelli/tesi/PDobs_2.csv')

#### sesPD ####
#randomizing my community or spps in the tree to obtain a Standardized effect size (SES)
#ottenendo il valore "effettivo" di pd a prescindere dalle spp che ci sono in quella cella
require(picante)
sesPD_2 = picante::ses.pd(samp=spp1pa_2, tree=myTree2, null.model ="independentswap",runs = 999, include.root=TRUE)

#salvare sesPD come rds o come tabella cvs, basta che possa essere importato in R come data.frame
write.csv(sesPD_2, '/Users/sofiaprandelli/tesi/sesPD_2.csv')

#### WORLDCLIM DOWNLOAD ####
setwd("/Users/sofiaprandelli/tesi")
library(raster)
Worldclim2 <- raster::getData('worldclim', var='bio', res=2.5) #risoluzione più fine, a 2.5 arc-minuti
plot(Worldclim2)
plot(m2, add=TRUE)
Worldclim2=crop(Worldclim2, m2)
plot(Worldclim2)

#importo la human population density, che però ha una risoluz di 1km, mentre io ho 2.5 arc-minuti(=4km)
pop_dens = raster("/Users/sofiaprandelli/tesi/POP_DENS_2020_1km_Aggregated.tif")

#aggrego il dato della popolazione a 4 km:
#uso RESAMPLE (raster) -> aggrego il dataset e lo aggiungo allo stack di Worldclim2 e poi lo aggiungo alle variabili del modello
r=aggregate(pop_dens, fact=4, fun=mean, expand=TRUE)
r=resample(r, Worldclim2, method ="bilinear")

#### all the climatic variables + pop density + sesPD, PD ####
tmp2 = raster::extract(r, m2, fun=mean, df=TRUE) #estrai la media di tutti i valori delle 19 variabili per le celle della griglia m2
tmp3 = raster::extract(Worldclim2, m2, fun=mean, df=TRUE)

#tmp$bio1 = tmp$bio1 / 10 #i valori di bio1 su Worldclim sono moltiplicati per 10
colnames(tmp2)[1] <- "id"
colnames(tmp2)[2] <- "pop_dens"
colnames(tmp3)[1] <- "id"

#merge tmp2 e oggetto spaziale m2
mydf5 <- sp::merge(tmp2, tmp3, by='id', all=F)
mydf5=cbind(m2, mydf5)
mydf5 <- mydf5@data
# aggiungo a PDobs e a sesPD il campo 'id'
require(tidyverse)
PDobs_2 <- add_column(PDobs_2, rownames(spp1pa_2), .before='PD')
colnames(PDobs_2)[1] <- "id"
sesPD_2 <- add_column(sesPD_2, rownames(spp1pa_2), .before='ntaxa')
colnames(sesPD_2)[1] <- "id"

mydf5=sp::merge(mydf5, PDobs_2, by='id')
mydf5=sp::merge(mydf5, sesPD_2, by='id')
#mydf5 <- mydf5[,-2]
#mydf5 <- mydf5[,-23:-24]
#mydf5 <- filter(mydf5, bio1 !="NA")
#colnames(mydf5) [c(21,22,28)] <- c("PD", "SR", "sesPD")
write.csv(mydf5, '/Users/sofiaprandelli/tesi/mydf5+pop_dens.csv')
mydf5_complete <- na.omit(mydf5)
table(mydf5_complete$sesPD>="0.00") #TRUE 193, FALSE 1364

#### ggplot2::facet_wrap() SCATTERPLOT sesPD vs ALL VARIABLES ####
library(readxl)
library(ggplot2)
names(mydf5_complete)

mydf5_long=mydf5_complete %>%
  select("id", "pop_dens", "bio1","bio2",
         "bio3", "bio4", "bio5", "bio6", "bio7",
         "bio8", "bio9", "bio10", "bio11", "bio12",
         "bio13", "bio14", "bio15", "bio16", "bio17", "bio18",
         "bio19","sesPD", "SR") %>%
  rename(myID="id") %>%
  pivot_longer(-c(myID, sesPD, SR)) %>%
  mutate(name=factor(name, levels = c("pop_dens", "bio1","bio2",
                                      "bio3", "bio4", "bio5", "bio6", "bio7",
                                      "bio8", "bio9", "bio10", "bio11", "bio12",
                                      "bio13", "bio14", "bio15", "bio16", "bio17", "bio18",
                                      "bio19")))
sf3 <- mydf5_long %>%
  ggplot(aes(x=value, y=sesPD))+
#  geom_line()+
  geom_point(size = 0.8)+
  labs(title = "Standardized effect size of Phylogenetic Diversity variations",
       x="Environmental variables", y="sesPD")+
  facet_wrap(~name, scales = "free")+
  theme_light() -> pd2;pd2
tiff('/Users/sofiaprandelli/tesi/PLOTS/sesPD_BIOVAR',res=300,units='cm',width = 30,height = 25)
pd2
dev.off()

sf3 + theme(plot.title = element_text(size = 20), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) 

#### VARIABLE IMPORTANCE BOX PLOT ####
df<-data.frame(as.matrix(myRF3$variable.importance)) 
df$variable<-rownames(df)
colnames(df)[1]<-'importance'
ggplot(df, aes(x=reorder(variable,importance), y=importance,fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  ggtitle("Information Value Summary")+
  scale_fill_gradient(low="red", high="blue")+
  theme_light()+theme(legend.position = 'none')

#### Random Forest model - Extent World ####
library(raster)
Worldclim <- raster::getData('worldclim', var='bio', res=5)

pop_dens = raster("/Users/sofiaprandelli/tesi/POP_DENS_2020_1km_Aggregated.tif")
rr=aggregate(pop_dens, fact=9, fun=mean, expand=TRUE)
rr=resample(rr, Worldclim, method ="bilinear")

myRF3=ranger(sesPD~bio1+bio2+bio3+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19+pop_dens, importance = 'permutation', data=mydf5_complete)
new.stack = stack(Worldclim[[c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")]], rr)
names(new.stack)[20] <- "pop_dens"

require(utils)
myPredR = predict(new.stack, myRF3, type = "response", predict.all=FALSE, na.rm = T, progress = "text", fun = function(model, ...) predict(model, ...)$predictions)

devtools::install_github("thomasp85/scico")
install.packages("scico")
library(scico)
library(ggplot2)
library(rasterVis)

pred_world <- gplot(myPredR, maxpixels=500000) +
  geom_raster(aes(fill = value), interpolate = TRUE, color = "black") + 
  labs(x="Longitude", y="Latitude", fill="") +
  ggtitle("SesPD prediction - World extent") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14)) + 
  theme(legend.text=element_text(size=12)) +
  theme_light() +
  scale_fill_viridis_c(limits=c(-0.5, 0.5), na.value="transparent") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 12, barheight = 0.5, size=10),
         size = guide_legend(title.position="top", title.hjust = 0.5)) +
  coord_equal()

plot(pred_world)

#### RANDOM FOREST MODEL #### Extent Europe, resolution 2.5=about 4.5 km at the equator
library(ranger)

myRF3=ranger(sesPD~bio1+bio2+bio3+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19+pop_dens, importance = 'permutation', data=mydf5_complete)
print(myRF3)
importance(myRF3)

#### Prediction - bisector and intercept
data4plot=cbind.data.frame("pred"=myRF3$predictions, "obs"=mydf5_complete$sesPD)
summary(lm(data4plot$obs~data4plot$pred))
predlm <- (lm(data4plot$obs~data4plot$pred))

ggplot(data=data4plot, aes(x=pred, y=obs) ) +
  geom_point() +
  geom_abline(slope = 1) +
  geom_abline(intercept=predlm$coefficients[1], slope=predlm$coefficients[2], col="blue") +
  geom_hline(yintercept = 0,  linetype = 2, colour = "darkgray") +
  geom_vline(xintercept = 0,  linetype = 2, colour = "darkgray") +
  ylab("sesPD") +
  xlab("Predictions") +
  xlim(-0.7 , 0.4) +
  ylim(-0.7 , 0.4) +
  annotate(geom="text", x=0.23, y=0.30, label="y=x",color="black")+
  annotate(geom="text", x=0.23, y=0.06, label=paste("R2=",round(tmpp$r.squared, 2)), color="blue")+
  ggtitle("Prediction of Standardized Effect Size of Phylgenetic Diversity 
          based on Bioclimatic and Anthropogenic Variables") +
  #  annotate("R squared (OOB): 0.1701877") +     paste("R2 (OBB):", round(myRF3$R2, 3))
#  scale_fill_continuous(type = "viridis") +
  #  geom_abline(intercept = 0, slope = 1) +
  theme_bw()
sesPD_prediction_3 + geom_abline(slope = 1)

#### Prediction - hex bins geom
ggplot(data=data4plot, aes(x=pred, y=obs) ) +
  geom_hex(bins = 70) +
  ylab("sesPD") +
  xlab("Predictions") +
  xlim(-0.7 , 0.4) +
  ylim(-0.7 , 0.4) +
  ggtitle("Prediction of Standardized Effect Size of Phylgenetic Diversity 
     based on Bioclimatic and Anthropogenic Variables") +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

#### PREDICTION ####
new.stack3 = stack(Worldclim2[[c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")]], r)
names(new.stack3)[20] <- "pop_dens"

require(utils)
#myPredR2 = predict(new.stack2, myRF2, type = "response", predict.all=FALSE, na.rm = T, progress = "text", fun = function(model, ...) predict(model, ...)$predictions)
#plot(myPredR2, main = "Prediction of sesPD based on
#     Random Forest model")

myPredR3 = predict(new.stack3, myRF3, type = "response", predict.all=FALSE, na.rm = T, progress = "text", fun = function(model, ...) predict(model, ...)$predictions)
myPredR3

#### FINAL RASTER ####
# 1 (prova)
gplot(myPredR3, maxpixels=500000) +
  geom_raster(aes(fill = value), interpolate = TRUE, color = "black") + 
  labs(x="Longitude",y="Latitude", fill="")+
  theme_light()+
  theme(legend.position = "bottom")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 20, barheight = 0.8),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  coord_equal()

#### 2 (prova)
gplot(myPredR3, maxpixels=500000) +
  geom_raster(aes(fill = value), interpolate = TRUE, color = "black") + 
  labs(x="Longitude", y="Latitude", fill="")+
  theme_light()+
  scale_fill_gradientn(colors=rainbow(100)) +
  theme(legend.position = "bottom")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 20, barheight = 0.8),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  coord_equal()

#### 3 (prova)
gplot(myPredR3, maxpixels=500000) +
  geom_raster(aes(fill = value), interpolate = TRUE, color = "black") + 
  labs(x="Longitude", y="Latitude", fill="")+
  theme_light()+
  scale_fill_viridis() +
  theme(legend.position = "bottom")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 20, barheight = 0.8),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  coord_equal()

#### 4 (SCALE FILL SCICO or VIRIDIS)
devtools::install_github("thomasp85/scico")
install.packages("scico")
library(scico)
library(ggplot2)
library(rasterVis)
library(viridis)
#scico_palette_show()

scico_pred_europe <- gplot(myPredR3, maxpixels=500000) +
  geom_raster(aes(fill = value), interpolate = TRUE, color = "black") + 
  labs(x="Longitude", y="Latitude", fill="", title="sesPD prediction - European extent")+
  ggtitle("SesPD prediction - European extent")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14)) + 
  theme(legend.text=element_text(size=12)) +
  theme_light()+
#  scale_fill_scico(palette = 'roma', limits=c(-0.5, 0.5), na.value="transparent") +
  scale_fill_viridis_c(limits=c(-0.5, 0.5), na.value="transparent") +
  theme_bw() +
  theme(legend.position = "bottom")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 12, barheight = 0.5, size=10),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  coord_equal()

grid.arrange(pred_world, scico_pred_europe, ncol=1, nrow=2)
                   
