# setup -------------------------------------------------------------------

rm(list=ls())

library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(funspace)
library(TPD)
library(ggpubr)
library(broom)
library(nlme)
library(ape)

print(citation("funspace"), bibtex=T)
print(citation("TPD"), bibtex=T)
print(citation("nlme"), bibtex=T)
print(citation("ape"), bibtex=T)

# data import -------------------------------------------------------------

read_me<-list.files(path="Data/plots", full.names = TRUE)
read_me2<-strsplit(read_me, "/" ) |>
  sapply("[",3)
read_me2<-gsub(".csv","",read_me2)
dat=list()
for(i in seq_along(read_me)){
  df<-read.csv(read_me[i], header = T, sep=",", col.names = c("N", "Species", "uh", "x", "y", "soil depth", "height",
                                                              "LS", "SLA", "LA", "notes"), colClasses =
                 c("NULL","factor","NULL","numeric", "numeric",
                   "numeric", "numeric", "numeric", "numeric",
                   "numeric", "factor"),skip = 6)
  df <- subset(df, notes!="seedling")
  df$Species <- str_trim(df$Species)
  df$Species <- as.factor(df$Species)
  df <- df[,-9]
  df <- df |> select(x,y,Species,soil.depth:LA)
  df <- df |> drop_na(x,y)
  df <- df |> mutate(LA = na_if(LA, 0))
  dat[[i]]<-df
  names(dat)[i]<-read_me2[i]
}
rm(i, df)

specieslookup <- read.csv("Data/specieslist.csv",header=T,sep=",")
specieslookup[2:5] <- lapply(specieslookup[2:5], as.logical)

dat <- lapply(dat, function(x) left_join(x, specieslookup))
rm(specieslookup, read_me, read_me2)

# multi-trait with PCA ----------------------------------------------------

datPCA <- lapply(dat, na.omit)
datPCA <- map(datPCA, ~.x[,c(-4,-6)])
dattt <- do.call(rbind, datPCA)
dattt <- dattt |> rename(Height=height)

#asymmetric distribution, log transform
hist(log10(dattt$Height))
hist(log10(dattt$SLA))
hist(log10(dattt$LA))

dattt <- dattt |> 
  mutate(across(4:6, log10))

#how many dimensions?

funspaceDim(dattt[4:6]) #2

pca_results <- princomp(dattt[,c(4:6)], cor = T)  #same results as prcomp, down to numerical accuracy
#it is allowed to invert both loadings and scores
pca_results$loadings[,1] <- -pca_results$loadings[,1]
pca_results$scores[,1] <- -pca_results$scores[,1]
summary(pca_results) #92% variance

dattt <- cbind(dattt, pca_results$scores)
dattt <- dattt |> mutate(plot=str_extract(rownames(dattt),".*?(?=\\.)")) |> 
  mutate(plot=as.numeric(as.factor(plot))) |> 
  mutate(Species=as.character(Species)) |> 
  mutate(Comm=ifelse(plot <= 5, 0, 1)) |> 
  mutate(Comm=factor(Comm, labels = c("Closed", "Open")))

# funspace ---------------------------------------------------------------------

#use princomp
#global TPD

plot(funspace(pca_results), quant.plot=T, arrows=T)
dev.off()

## by community----

funcommunity <- funspace(pca_results, group.vec = dattt$Comm)

#for submission

postscript('funspacelog.eps', onefile = T, bg = "white", paper = "special",
           width = 7, height = 3.5, pointsize = 8, horizontal = F)
par(mfrow=c(1,2), mar = c(par("mar")[1], par("mar")[2], par("mar")[4], par("mar")[4]), bg= "white")
plot(funcommunity, type = "groups", globalContour=T, globalContour.lwd = 2,quant.plot = T,
     quant = c(0.99,0.95,0.5), arrows = T, arrows.length = 0.8, arrows.label.pos = 1.2)
dev.off()

#for manuscript

png('funspacelog.png', bg = "white", width = 7, height = 3.5, units = "in",
    pointsize = 8, res = 1000)
par(mfrow=c(1,2), mar = c(par("mar")[1], par("mar")[2], par("mar")[4], par("mar")[4]), bg= "white")
plot(funcommunity, type = "groups", globalContour=T, globalContour.lwd = 2,quant.plot = T,
     quant = c(0.99,0.95,0.5), arrows = T, arrows.length = 0.8, arrows.label.pos = 1.2)
dev.off()

#for slides

png('funspacelog_nobg.png', width = 14.61, height = 13.75, units = "cm", res=300)
par(mfrow=c(2,1), mar=c(3, par("mar")[2], 1, par("mar")[4]), bg= NA) #set to NA for slides
plot(funcommunity, type = "groups", globalContour=T, globalContour.lwd = 2,quant.plot = T,
     quant = c(0.99,0.95,0.5), arrows = T, arrows.length = 0.7, arrows.label.pos = 1.2)
dev.off()

##growth form-----

comb_form <- data.frame(expand.grid(forb=F:T, Comm=1:2), condition= c("Grass, Closed", "Forb, Closed",
                                                                      "Grass, Open", "Forb, Open"))
comb_form$forb <- as.logical(comb_form$forb)
comb_form$Comm <- factor(comb_form$Comm,labels = c("Closed", "Open"))

dattt1 <- left_join(dattt,comb_form)
dattt1 <- dattt1 |> filter(woody==F)
pca_results1 <- princomp(dattt1[,c(4:6)], cor = T)

funform <- funspace(pca_results1,group.vec = dattt1$condition)
par(mfrow=c(2,2))
plot(funform, type = "groups", globalContour=T, globalContour.lwd = 2,quant.plot = T,
     quant = c(0.99,0.95,0.5), arrows = T)
rm(comb_form)

summary(funform)

##clonality-----

comb_clon <- data.frame(expand.grid(clonal=F:T, Comm=1:2), condition= c("Nonclonal, Closed", "Clonal, Closed",
                                                                        "Nonclonal, Open", "Clonal, Open"))
comb_clon$clonal <- as.logical(comb_clon$clonal)
comb_clon$Comm <- factor(comb_clon$Comm,labels = c("Closed", "Open"))
dattt2 <- left_join(dattt,comb_clon)
rm(comb_clon)
pca_results <- princomp(dattt2[,c(4:6)], cor = T)

funclone <- funspace(pca_results,group.vec = dattt2$condition)
par(mfrow=c(2,2))
plot(funclone, type = "groups", globalContour=T, globalContour.lwd = 2,quant.plot = T,
     quant = c(0.99,0.95,0.5), arrows = T)


# TPD ---------------------------------------------------------------------

speciesn_PCA <- dattt |> 
  group_by(Comm, Species) |> 
  summarise(count = n()) |> 
  ungroup()

#TPD in two dimensions requires at least 3 observations
dattt_TPD <- dattt |> 
  left_join(speciesn_PCA, by = c("Comm", "Species")) |> 
  filter(count > 2) |> 
  select(-count)

#for single traits

dattt_norm <- do.call(rbind, dat) |> 
  relocate(LS, .before = height)

#log10-transform
dattt_norm <- dattt_norm |> 
  mutate(across(6:8, log10))

dattt_norm <- dattt_norm |> mutate(plot=str_extract(rownames(dattt_norm),".*?(?=\\.)")) |> 
  mutate(plot=as.numeric(as.factor(plot))) |> 
  mutate(Species=as.character(Species)) |> 
  mutate(Comm=ifelse(plot <= 5, 0, 1)) |> 
  mutate(Comm=factor(Comm, labels = c("Closed", "Open")))

speciesn <- dattt_norm |> 
  group_by(Comm, Species) |> 
  summarise(count = n()) |> 
  ungroup()


##overall----

TPDs <- TPDs(dattt_TPD$Species, dattt_TPD[,c("Comp.1","Comp.2")], samples = dattt_TPD$Comm)

sampUnit <- dattt_TPD |> select(Species, Comm) |> group_by(Comm,Species) |>
  summarise(n=n()) |> 
  pivot_wider(names_from = Species, values_from = n, values_fill = 0) |> 
  as.data.frame()
row.names(sampUnit) <- sampUnit$Comm
sampUnit <- sampUnit[,-1]

TPDc <- TPDc(TPDs, sampUnit)

dissim_comm <- dissim(TPDc)

sink(file="traits/dissim_redun.txt")
dissim_comm
redundancy(TPDc)
sink()
rm(dissim_comm,TPDc,TPDs, sampUnit)

##closed----

###PCA----

dattt_CL <- dattt_TPD |> filter(Comm=="Closed")
species <- dattt_CL$Species
traits <- dattt_CL[,c("Comp.1","Comp.2")]

TPDs_CL <- TPDs(species = species, traits = traits)
rm(traits)

sampUnit <- data.frame(Species=species, Comm="Closed") |>
  count(Species) |> 
  pivot_wider(names_from = Species, values_from = n, values_fill = 0) |> 
  mutate(across(everything(), as.integer))
row.names(sampUnit) <- "Closed"

TPDc_CL <- TPDc(TPDs_CL, sampUnit)
rm(sampUnit)

#pairwise overlap

pairwiseOL_CL <- 1-dissim(TPDs_CL)[["populations"]][["dissimilarity"]]
pairwiseOL_CL[lower.tri(pairwiseOL_CL, diag = T)] <- 0
pairwiseOL_CL <- as.data.frame(as.table(pairwiseOL_CL)) |> 
  filter(Freq > 0)

props <- speciesn_PCA |>
  filter(count > 2, Comm=="Closed") |> 
  mutate(prop = count/sum(count)) |> 
  select(Species, prop)
weightedOL_PCA_CL <- pairwiseOL_CL |> 
  left_join(props, by=join_by(Var1 == Species)) |> 
  left_join(props, by=join_by(Var2 == Species)) |> 
  mutate(product = prop.x*prop.y) |> 
  summarise(sum(product)) |> 
  pull()

rm(pairwiseOL_CL,props)

#community overlap

overlap_CL <- 1-uniqueness(TPDs_CL, TPDc_CL)
overlap_CL <- data.frame(Species=colnames(overlap_CL),"Trait space overlap"=t(overlap_CL)) |> 
  rename("Trait space"="Closed") |> arrange(desc(`Trait space`))


rm(TPDc_CL,TPDs_CL, species)

###single traits----
#single traits require at least 2 observations

singletrait_CL <- dattt_norm |> 
  left_join(speciesn, by = c("Comm", "Species")) |> 
  filter(count > 1) |> 
  select(-count) |> 
  filter(Comm=="Closed")

#### height----
species <- singletrait_CL$Species[is.finite(singletrait_CL$height)]
traits <- singletrait_CL[is.finite(singletrait_CL$height),"height", drop=F]
TPDs_CL_H <- TPDs(species = species, traits = traits)
spatstat.geom::integral.density(density(traits$height[species=="Medicago lopulina"])) #0.999
rm(traits)

sampUnit <- data.frame(Species=species, Comm="Closed") |> 
  count(Species) |> 
  pivot_wider(names_from = Species, values_from = n, values_fill = 0) |> 
  mutate(across(everything(), as.integer))
row.names(sampUnit) <- "Closed"

TPDc_CL_H <- TPDc(TPDs_CL_H, sampUnit)
rm(sampUnit)

#pairwise overlap
pairwiseOL_CL <- 1-dissim(TPDs_CL_H)[["populations"]][["dissimilarity"]]
pairwiseOL_CL[lower.tri(pairwiseOL_CL, diag = T)] <- 0
pairwiseOL_CL <- as.data.frame(as.table(pairwiseOL_CL)) |> 
  filter(Freq > 0)

props <- speciesn |>
  filter(count > 1, Comm=="Closed") |> 
  mutate(prop = count/sum(count)) |> 
  select(Species, prop)
weightedOL_H_CL <- pairwiseOL_CL |> 
  left_join(props, by=join_by(Var1 == Species)) |> 
  left_join(props, by=join_by(Var2 == Species)) |> 
  mutate(product = prop.x*prop.y) |> 
  summarise(sum(product)) |> 
  pull()

rm(pairwiseOL_CL,props)

#community overlap
overlap_CL_H <- 1-uniqueness(TPDs_CL_H, TPDc_CL_H)
overlap_CL_H <- data.frame(Species=colnames(overlap_CL_H),"Trait space overlap"=t(overlap_CL_H)) |> 
  rename(Height="Closed") |> arrange(desc(Height))

overlaps_CL <- full_join(overlap_CL,overlap_CL_H)

rm(TPDs_CL_H,TPDc_CL_H, species)

#### LA----
species <- singletrait_CL$Species[is.finite(singletrait_CL$LA)]
traits <- singletrait_CL[is.finite(singletrait_CL$LA),"LA", drop=F]
TPDs_CL_LA <- TPDs(species = species, traits = traits)
spatstat.geom::integral.density(density(traits$LA[species=="Eryngium ametistinum"])) #0.999
rm(traits)

sampUnit <- data.frame(Species=species, Comm="Closed") |> 
  count(Species) |> 
  pivot_wider(names_from = Species, values_from = n, values_fill = 0) |> 
  mutate(across(everything(), as.integer))
row.names(sampUnit) <- "Closed"

TPDc_CL_LA <- TPDc(TPDs_CL_LA, sampUnit)
rm(sampUnit)

#pairwise overlap
pairwiseOL_CL <- 1-dissim(TPDs_CL_LA)[["populations"]][["dissimilarity"]]
pairwiseOL_CL[lower.tri(pairwiseOL_CL, diag = T)] <- 0
pairwiseOL_CL <- as.data.frame(as.table(pairwiseOL_CL)) |> 
  filter(Freq > 0)

props <- speciesn |>
  filter(count > 1, Comm=="Closed") |> 
  mutate(prop = count/sum(count)) |> 
  select(Species, prop)
weightedOL_LA_CL <- pairwiseOL_CL |> 
  left_join(props, by=join_by(Var1 == Species)) |> 
  left_join(props, by=join_by(Var2 == Species)) |> 
  mutate(product = prop.x*prop.y) |> 
  summarise(sum(product)) |> 
  pull()

rm(pairwiseOL_CL,props)

#community overlap
overlap_CL_LA <- 1-uniqueness(TPDs_CL_LA, TPDc_CL_LA)
overlap_CL_LA <- data.frame(Species=colnames(overlap_CL_LA),"Trait space overlap"=t(overlap_CL_LA)) |> 
  rename(LA="Closed") |> arrange(desc(LA))

overlaps_CL <- full_join(overlaps_CL,overlap_CL_LA)

rm(TPDs_CL_LA,TPDc_CL_LA, species)

####SLA----

species <- singletrait_CL$Species[is.finite(singletrait_CL$SLA)]
traits <- singletrait_CL[is.finite(singletrait_CL$SLA),"SLA", drop=F]
TPDs_CL_SLA <- TPDs(species = species, traits = traits)
spatstat.geom::integral.density(density(traits$SLA[species=="Plantago argentea"])) #0.999
rm(traits)

sampUnit <- data.frame(Species=species, Comm="Closed") |> 
  count(Species) |> 
  pivot_wider(names_from = Species, values_from = n, values_fill = 0) |> 
  mutate(across(everything(), as.integer))
row.names(sampUnit) <- "Closed"

TPDc_CL_SLA <- TPDc(TPDs_CL_SLA, sampUnit)
rm(sampUnit)

#pairwise overlap
pairwiseOL_CL <- 1-dissim(TPDs_CL_SLA)[["populations"]][["dissimilarity"]]
pairwiseOL_CL[lower.tri(pairwiseOL_CL, diag = T)] <- 0
pairwiseOL_CL <- as.data.frame(as.table(pairwiseOL_CL)) |> 
  filter(Freq > 0)

props <- speciesn |>
  filter(count > 1, Comm=="Closed") |> 
  mutate(prop = count/sum(count)) |> 
  select(Species, prop)
weightedOL_SLA_CL <- pairwiseOL_CL |> 
  left_join(props, by=join_by(Var1 == Species)) |> 
  left_join(props, by=join_by(Var2 == Species)) |> 
  mutate(product = prop.x*prop.y) |> 
  summarise(sum(product)) |> 
  pull()

rm(pairwiseOL_CL,props)


#community overlap
overlap_CL_SLA <- 1-uniqueness(TPDs_CL_SLA, TPDc_CL_SLA)
overlap_CL_SLA <- data.frame(Species=colnames(overlap_CL_SLA),"Trait space overlap"=t(overlap_CL_SLA)) |> 
  rename(SLA="Closed") |> arrange(desc(SLA))

overlaps_CL <- full_join(overlaps_CL,overlap_CL_SLA)

rm(TPDs_CL_SLA,TPDc_CL_SLA, species)

###wrapping up----

#community overlaps
overlaps_CL <- speciesn |> filter(Comm=="Closed") |> 
  right_join(overlaps_CL) |> 
  select(-Comm) |> 
  arrange(desc(count)) |> 
  rename(Abundance=count)

weights_CL <- data.frame(
  "Trait space"=
  overlaps_CL$Abundance/sum(overlaps_CL$Abundance[!is.na(overlaps_CL$`Trait space`)]),
  "Height"=overlaps_CL$Abundance/sum(overlaps_CL$Abundance[!is.na(overlaps_CL$`Height`)]),
  "LA"=overlaps_CL$Abundance/sum(overlaps_CL$Abundance[!is.na(overlaps_CL$`LA`)]),
  "SLA"=overlaps_CL$Abundance/sum(overlaps_CL$Abundance[!is.na(overlaps_CL$`SLA`)])
)


wm_CL <- map2(overlaps_CL[3:6], weights_CL, ~ weighted.mean(.x, w = .y, na.rm=T))
wm_CL <- c(list(Species= "Weighted mean", Abundance=NA), wm_CL)

overlaps_CL <- rbind(wm_CL,overlaps_CL)

#pairwise overlaps
pairwise_CL <- c(list("Species" = "Mean pairwise overlap",
                      "Abundance" = NA,
                      "Trait space"=weightedOL_PCA_CL,
                      "Height" = weightedOL_H_CL,
                      "LA" = weightedOL_LA_CL,
                      "SLA" = weightedOL_SLA_CL))

overlaps_CL <- rbind(pairwise_CL, overlaps_CL)
overlaps_CL$Abundance <- as.numeric(overlaps_CL$Abundance)

rm(sampUnit,dattt_CL,singletrait_CL,overlap_CL,overlap_CL_H,overlap_CL_LA,overlap_CL_SLA,
   wm_CL, weights_CL, pairwise_CL)
rm(list = ls(pattern = "weightedOL"))


##open----

###PCA----
dattt_OP <- dattt_TPD |> filter(Comm=="Open")
species <- dattt_OP$Species
traits <- dattt_OP[,c("Comp.1","Comp.2")]

TPDs_OP <- TPDs(species = species, traits = traits)
rm(traits)

sampUnit <- data.frame(Species=species, Comm="Open") |>
  count(Species) |> 
  pivot_wider(names_from = Species, values_from = n, values_fill = 0) |> 
  mutate(across(everything(), as.integer))
row.names(sampUnit) <- "Open"

TPDc_OP <- TPDc(TPDs_OP, sampUnit)
rm(sampUnit)

#pairwise overlap
pairwiseOL_OP <- 1-dissim(TPDs_OP)[["populations"]][["dissimilarity"]]
pairwiseOL_OP[lower.tri(pairwiseOL_OP, diag = T)] <- 0
pairwiseOL_OP <- as.data.frame(as.table(pairwiseOL_OP)) |> 
  filter(Freq > 0)

props <- speciesn_PCA |>
  filter(count > 2, Comm=="Open") |> 
  mutate(prop = count/sum(count)) |> 
  select(Species, prop)
weightedOL_PCA_OP <- pairwiseOL_OP |> 
  left_join(props, by=join_by(Var1 == Species)) |> 
  left_join(props, by=join_by(Var2 == Species)) |> 
  mutate(product = prop.x*prop.y) |> 
  summarise(sum(product)) |> 
  pull()

rm(pairwiseOL_OP,props)

#community overlap
overlap_OP <- 1-uniqueness(TPDs_OP, TPDc_OP)
overlap_OP <- data.frame(Species=colnames(overlap_OP),"Trait space overlap"=t(overlap_OP)) |> 
  rename("Trait space"="Open") |> arrange(desc(`Trait space`))

rm(TPDc_OP,TPDs_OP, species)

###single traits----

singletrait_OP <- dattt_norm |> 
  left_join(speciesn, by = c("Comm", "Species")) |> 
  filter(count > 1) |> 
  select(-count) |> 
  filter(Comm=="Open")

####height----
species <- singletrait_OP$Species[is.finite(singletrait_OP$height)]
traits <- singletrait_OP[is.finite(singletrait_OP$height),"height", drop=F]
TPDs_OP_H <- TPDs(species = species, traits = traits)
rm(traits)

sampUnit <- data.frame(Species=species, Comm="Open") |> 
  count(Species) |> 
  pivot_wider(names_from = Species, values_from = n, values_fill = 0) |> 
  mutate(across(everything(), as.integer))
row.names(sampUnit) <- "Open"

TPDc_OP_H <- TPDc(TPDs_OP_H, sampUnit)
rm(sampUnit)

#pairwise overlap
pairwiseOL_OP <- 1-dissim(TPDs_OP_H)[["populations"]][["dissimilarity"]]
pairwiseOL_OP[lower.tri(pairwiseOL_OP, diag = T)] <- 0
pairwiseOL_OP <- as.data.frame(as.table(pairwiseOL_OP)) |> 
  filter(Freq > 0)

props <- speciesn |>
  filter(count > 1, Comm=="Open") |> 
  mutate(prop = count/sum(count)) |> 
  select(Species, prop)
weightedOL_H_OP <- pairwiseOL_OP |> 
  left_join(props, by=join_by(Var1 == Species)) |> 
  left_join(props, by=join_by(Var2 == Species)) |> 
  mutate(product = prop.x*prop.y) |> 
  summarise(sum(product)) |> 
  pull()

rm(pairwiseOL_OP,props)

#community overlap
overlap_OP_H <- 1-uniqueness(TPDs_OP_H, TPDc_OP_H)
overlap_OP_H <- data.frame(Species=colnames(overlap_OP_H),"Trait space overlap"=t(overlap_OP_H)) |> 
  rename(Height="Open") |> arrange(desc(Height))

overlaps_OP <- full_join(overlap_OP,overlap_OP_H)

rm(TPDs_OP_H,TPDc_OP_H, species)

####LA----
species <- singletrait_OP$Species[is.finite(singletrait_OP$LA)]
traits <- singletrait_OP[is.finite(singletrait_OP$LA),"LA", drop=F]
traits <- traits[species!="Dianthus sylvestris",, drop=F] #only one observation without NA
species <- species[species!="Dianthus sylvestris"]
TPDs_OP_LA <- TPDs(species = species, traits = traits)
rm(traits)

sampUnit <- data.frame(Species=species, Comm="Open") |> 
  count(Species) |> 
  pivot_wider(names_from = Species, values_from = n, values_fill = 0) |> 
  mutate(across(everything(), as.integer))
row.names(sampUnit) <- "Open"

TPDc_OP_LA <- TPDc(TPDs_OP_LA, sampUnit)
rm(sampUnit)

#pairwise overlap
pairwiseOL_OP <- 1-dissim(TPDs_OP_LA)[["populations"]][["dissimilarity"]]
pairwiseOL_OP[lower.tri(pairwiseOL_OP, diag = T)] <- 0
pairwiseOL_OP <- as.data.frame(as.table(pairwiseOL_OP)) |> 
  filter(Freq > 0)

props <- speciesn |>
  filter(count > 1, Comm=="Open", Species != "Dianthus sylvestris") |> 
  mutate(prop = count/sum(count)) |> 
  select(Species, prop)
weightedOL_LA_OP <- pairwiseOL_OP |> 
  left_join(props, by=join_by(Var1 == Species)) |> 
  left_join(props, by=join_by(Var2 == Species)) |> 
  mutate(product = prop.x*prop.y) |> 
  summarise(sum(product)) |> 
  pull()

rm(pairwiseOL_OP,props)

#community overlap
overlap_OP_LA <- 1-uniqueness(TPDs_OP_LA, TPDc_OP_LA)
overlap_OP_LA <- data.frame(Species=colnames(overlap_OP_LA),"Trait space overlap"=t(overlap_OP_LA)) |> 
  rename(LA="Open") |> arrange(desc(LA))

overlaps_OP <- full_join(overlaps_OP,overlap_OP_LA)

rm(TPDs_OP_LA,TPDc_OP_LA, species)

####SLA----

species <- singletrait_OP$Species[is.finite(singletrait_OP$SLA)]
traits <- singletrait_OP[is.finite(singletrait_OP$SLA),"SLA", drop=F]
traits <- traits[species!="Dianthus sylvestris",, drop=F] #only one observation without NA
species <- species[species!="Dianthus sylvestris"]
TPDs_OP_SLA <- TPDs(species = species, traits = traits)
rm(traits)

sampUnit <- data.frame(Species=species, Comm="Open") |> 
  count(Species) |> 
  pivot_wider(names_from = Species, values_from = n, values_fill = 0) |> 
  mutate(across(everything(), as.integer))
row.names(sampUnit) <- "Open"

TPDc_OP_SLA <- TPDc(TPDs_OP_SLA, sampUnit)
rm(sampUnit)

#pairwise overlap
pairwiseOL_OP <- 1-dissim(TPDs_OP_SLA)[["populations"]][["dissimilarity"]]
pairwiseOL_OP[lower.tri(pairwiseOL_OP, diag = T)] <- 0
pairwiseOL_OP <- as.data.frame(as.table(pairwiseOL_OP)) |> 
  filter(Freq > 0)

props <- speciesn |>
  filter(count > 1, Comm=="Open", Species!="Dianthus sylvestris") |> 
  mutate(prop = count/sum(count)) |> 
  select(Species, prop)
weightedOL_SLA_OP <- pairwiseOL_OP |> 
  left_join(props, by=join_by(Var1 == Species)) |> 
  left_join(props, by=join_by(Var2 == Species)) |> 
  mutate(product = prop.x*prop.y) |> 
  summarise(sum(product)) |> 
  pull()

rm(pairwiseOL_OP,props)

#community overlap
overlap_OP_SLA <- 1-uniqueness(TPDs_OP_SLA, TPDc_OP_SLA)
overlap_OP_SLA <- data.frame(Species=colnames(overlap_OP_SLA),"Trait space overlap"=t(overlap_OP_SLA)) |> 
  rename(SLA="Open") |> arrange(desc(SLA))

overlaps_OP <- full_join(overlaps_OP,overlap_OP_SLA)

rm(TPDs_OP_SLA,TPDc_OP_SLA, species)

###wrapping up----

#community overlaps
overlaps_OP <- speciesn |> filter(Comm=="Open") |> 
  right_join(overlaps_OP) |> 
  select(-Comm) |> 
  arrange(desc(count)) |> 
  rename(Abundance=count)

weights_OP <- data.frame(
  "Trait space"=
    overlaps_OP$Abundance/sum(overlaps_OP$Abundance[!is.na(overlaps_OP$`Trait space`)]),
  "Height"=overlaps_OP$Abundance/sum(overlaps_OP$Abundance[!is.na(overlaps_OP$`Height`)]),
  "LA"=overlaps_OP$Abundance/sum(overlaps_OP$Abundance[!is.na(overlaps_OP$`LA`)]),
  "SLA"=overlaps_OP$Abundance/sum(overlaps_OP$Abundance[!is.na(overlaps_OP$`SLA`)])
)


wm_OP <- map2(overlaps_OP[3:6], weights_OP, ~ weighted.mean(.x, w = .y, na.rm=T))
wm_OP <- c(list(Species= "Weighted mean", Abundance=NA), wm_OP)

overlaps_OP <- rbind(wm_OP,overlaps_OP)

#pairwise overlaps
pairwise_OP <- tibble("Species" = "Mean pairwise overlap",
                      "Abundance" = NA,
                      "Trait space"=weightedOL_PCA_OP,
                      "Height" = weightedOL_H_OP,
                      "LA" = weightedOL_LA_OP,
                      "SLA" = weightedOL_SLA_OP)

overlaps_OP <- rbind(pairwise_OP, overlaps_OP)
overlaps_OP$Abundance <- as.numeric(overlaps_OP$Abundance)

rm(sampUnit,dattt_OP,singletrait_OP,overlap_OP,overlap_OP_H,overlap_OP_LA,overlap_OP_SLA,
   wm_OP, weights_OP, pairwise_OP)
rm(list = ls(pattern = "weightedOL"))

##export to excel----
writexl::write_xlsx(list("Closed"=overlaps_CL,"Open"=overlaps_OP), "traits/Overlaps_TPD.xlsx", col_names = T, format_headers = T)


# trait correlation -------------------------------------------------------

ggqqplot(dattt_norm$height)
ggqqplot(dattt_norm$LA)
ggqqplot(dattt_norm$SLA)

dattt_norm <- dattt_norm |> rename("Soil depth" = soil.depth)

correlations <- as.data.frame(t(combn(c("LA","SLA","height","LS", "Soil depth"),2))) |> 
  mutate(correlation = map2(.data$V1, .data$V2, 
                            \(uno,due) cor.test(dattt_norm[[uno]], 
                                                dattt_norm[[due]], 
                                                method = "spearman") |> 
                              tidy())
  ) |> 
  unnest_wider(correlation) |> 
  select(V1, V2, estimate, p.value, method) |> 
  mutate(Variables = paste0(V1, ", ", V2), .before = 1) |> 
  select(!c(V1,V2)) |> 
  arrange(desc(abs(estimate)))

write.table(correlations, "traits/correlations.txt")

# variance partitioning ----

#height
height_lm_obs <- lme(fixed = height~1, data = dattt_norm, na.action = na.omit,
                 random = ~ 1|Comm/plot/Species)

height_varcomp_obs <- varcomp(height_lm_obs, scale = T)

#LA
LA_lm <- lme(fixed = LA~1, data = dattt_norm, na.action = na.omit,
                 random = ~ 1|Comm/plot/Species)

LA_varcomp <- varcomp(LA_lm, scale = T)

#SLA
SLA_lm <- lme(fixed = SLA~1, data = dattt_norm, na.action = na.omit,
                 random = ~ 1|Comm/plot/Species)

SLA_varcomp <- varcomp(SLA_lm, scale = T)

sink("varcomp.txt")
"height"
round(height_varcomp, digits = 3)
"LA"
round(LA_varcomp, digits = 3)
"SLA"
round(SLA_varcomp, digits = 3)
sink()

# average traits----

conf_int <- function (x){
  x <- na.omit(x)
  se <- plotrix::std.error(x, na.rm = TRUE)
  tscore <- qt(p=0.025, df=length(x)-1, lower.tail = F)
  tscore*se
}

summary_fun <- list(
  mean = ~mean(.x, na.rm = TRUE), 
  CI = ~conf_int(.x)
)

dattt_nolog <- bind_rows(dat, .id = "plot") |>
  mutate(plot=as.numeric(as.factor(plot))) |> 
  mutate(Species=as.character(Species)) |> 
  mutate(Comm=ifelse(plot <= 5, 0, 1)) |> 
  mutate(Comm=factor(Comm, labels = c("Closed", "Open")))


trait_summary <- dattt_nolog |> group_by(Comm) |> summarise(across(5:9, summary_fun))
trait_summary <- trait_summary |> 
  pivot_longer(cols = !Comm, names_to = c("trait", ".value"), names_sep = "_")

write.table(trait_summary,"traits/traitsummary.txt")

library(ggplot2)

#pointrange
ggplot(trait_summary)+
  geom_pointrange(aes(x=trait, y=mean, ymin=CIl, ymax=CIu, group=Comm, color=Comm),
                   position = position_dodge(width=0.2))+
  facet_wrap(~trait)

#boxplot
dattt_nolog |> 
  select(4:8 | 14) |> 
  pivot_longer(!Comm, names_to = "trait") |> 
ggplot()+
  geom_boxplot(aes(x=Comm, y=value))+
  facet_wrap(vars(trait), scales = "free")
