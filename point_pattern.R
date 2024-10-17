#libraries----
library(stringr)
library(dplyr)
library(tidyr)
library(spatstat)
library(purrr)
library(GET)
library(ggplot2)
library(ggfortify)

rm(list=ls())

print(citation("spatstat"), bibtex=T)
print(citation("GET"), bibtex=T)

# read data ---------------------------------------------------------------

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
rm(i, df, read_me)

specieslookup <- read.csv("Data/specieslist.csv",header=T,sep=",")
specieslookup[2:5] <- lapply(specieslookup[2:5], as.logical)

dat <- lapply(dat, function(x) left_join(x, specieslookup))
rm(specieslookup)

# spatstat list -----------------------------------------------------------

#list of all
datpp <- list()
rejects <- list()
for (i in 1:10){
  datpp[[i]] <- as.ppp(dat[[i]], owin(c(0,50), c(0,50), unitname = "cm"))
  rejects[[i]] <- attr(datpp[[i]], "rejects")
  names(datpp)[i]<-read_me2[i]
}
rm(i,rejects)
datpp <- as.solist(datpp)

#ggplot of one unit----
library(ggplot2)

dat[[3]] |>
  ggplot(aes(x, y)) +
  xlim(c(0,50))+
  ylim(c(0,50))+
  geom_point(aes(color=Species, size=height))+
  scale_y_continuous(position = "right")+
  theme_light()+
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill='transparent', color=NA))+
  labs(x= "x (cm)", y= "y (cm)")
#ggsave('closed_3_right.png', bg="transparent", width = 14, height = 14, units = "cm", dpi = "print")

# find duplicates -------------------------------
dup <- list()
for (i in 1:10){
  temp <- as.ppp(dat[[i]], owin(c(0,50), c(0,50)))
  dup[[i]] <- which(duplicated(unmark(temp))==T)
}
rm(dup,temp,i)

#quadrat test----
quadrats <- lapply(datpp, function(x) {quadrat.test(x, 4, 4, method = "MonteCarlo", nsim = 1999)})
quadrattest <- purrr::map_dbl(quadrats, "p.value")
which(quadrattest<=0.05) #7 is inhomogeneous
plot(unmark(datpp[[7]])); plot(quadrats[[7]], add=T)
rm(quadrats)

# isolated marks ----

soildepth <- solist()
for(i in 1:10){
  soildepth[[i]] <- datpp[[i]] %mark% marks(datpp[[i]])["soil.depth"]
  names(soildepth)[i] <- read_me2[i]
}
rm(i)

height <- solist()
for(i in 1:10){
  height[[i]] <- datpp[[i]] %mark% marks(datpp[[i]])[c("Species","height")]
  names(height)[i] <- read_me2[i]
}
rm(i)

SLA <- solist()
for(i in 1:10){
  SLA[[i]] <- datpp[[i]] %mark% marks(datpp[[i]])[c("Species","SLA")]
  names(SLA)[i] <- read_me2[i]
}
rm(i)
SLA <- lapply(SLA, function(x) subset(x, SLA >= 0))

LA <- solist()
for(i in 1:10){
  LA[[i]] <- datpp[[i]] %mark% marks(datpp[[i]])[c("Species","LA")]
  names(LA)[i] <- read_me2[i]
}
rm(i)
LA <- lapply(LA, function(x) subset(x, LA >= 0))


# NAs in marks----

marks <- list()
for (i in 1:10){
  marks[[i]] <- marks(datpp[[i]])
}
rm(i)

NAmarks <- bind_rows(marks, .id = "plot") |> 
  filter(if_any(everything(), is.na))

rm(marks)

# kernel smoother of marks ----

##height-------------

onlyheight <- lapply(dat, function (df) df |> 
                   select(x, y, height) |> 
                   drop_na())

for (i in 1:10){
  onlyheight[[i]] <- as.ppp(onlyheight[[i]], owin(c(0,50), c(0,50), unitname = "cm"))
  names(onlyheight)[i]<-read_me2[i]
}
rm(i)

heightsmooth <- lapply(onlyheight[setdiff(1:10, 7)], function(x) Smooth(x, hmax=50))

heightdev <- lapply(heightsmooth, summary)

heightsum <- lapply(dat[setdiff(1:10, 7)], function (x) x |> 
                  select(height) |> 
                  summarise(mean=mean(height, na.rm = T), sd=sd(height, na.rm = T))) |> 
  bind_rows()

heightscores <- vector(mode = "numeric", length = length(heightsmooth))
for(i in seq_along(heightsmooth)){
  heightscores[[i]] <- max(abs(heightdev[[i]][["min"]]-heightsum[i,"mean"]),
                       abs(heightdev[[i]][["max"]]-heightsum[i,"mean"]))/heightsum[i,"sd"]
}

rm(onlyheight, i, heightdev, heightsum)

##LA-------------

onlyLA <- lapply(dat, function (df) df |> 
  select(x, y, LA) |> 
  drop_na())

for (i in 1:10){
  onlyLA[[i]] <- as.ppp(onlyLA[[i]], owin(c(0,50), c(0,50), unitname = "cm"))
  names(onlyLA)[i]<-read_me2[i]
}
rm(i)

LAsmooth <- lapply(onlyLA[setdiff(1:10, 7)], function(x) Smooth(x, hmax=50))

LAdev <- lapply(LAsmooth, summary)

LAsum <- lapply(dat[setdiff(1:10, 7)], function (x) x |> 
                      select(LA) |> 
                      summarise(mean=mean(LA, na.rm = T), sd=sd(LA, na.rm = T))) |> 
  bind_rows()

LAscores <- vector(mode = "numeric", length = length(LAsmooth))
for(i in seq_along(LAsmooth)){
  LAscores[[i]] <- max(abs(LAdev[[i]][["min"]]-LAsum[i,"mean"]),
                           abs(LAdev[[i]][["max"]]-LAsum[i,"mean"]))/LAsum[i,"sd"]
}

rm(onlyLA, i, LAdev, LAsum)

##SLA-------------

onlySLA <- lapply(dat, function (df) df |> 
                   select(x, y, SLA) |> 
                   drop_na())

for (i in 1:10){
  onlySLA[[i]] <- as.ppp(onlySLA[[i]], owin(c(0,50), c(0,50), unitname = "cm"))
  names(onlySLA)[i]<-read_me2[i]
}
rm(i)

SLAsmooth <- lapply(onlySLA[setdiff(1:10, 7)], function(x) Smooth(x, hmax=50))

SLAdev <- lapply(SLAsmooth, summary)

SLAsum <- lapply(dat[setdiff(1:10, 7)], function (x) x |> 
                  select(SLA) |> 
                  summarise(mean=mean(SLA, na.rm = T), sd=sd(SLA, na.rm = T))) |> 
  bind_rows()

SLAscores <- vector(mode = "numeric", length = length(SLAsmooth))
for(i in seq_along(SLAsmooth)){
  SLAscores[[i]] <- max(abs(SLAdev[[i]][["min"]]-SLAsum[i,"mean"]),
                       abs(SLAdev[[i]][["max"]]-SLAsum[i,"mean"]))/SLAsum[i,"sd"]
}

rm(onlySLA, i, SLAdev, SLAsum)

## results output----
smoothplot <- function (image, score){
  plotnames <- names(image) |> str_replace("-plot","")
  plot.listof(image, main = "",
              main.panel = paste0(plotnames, ", ", "z = ", round(score, digits = 2)))
}

png('heightsmooth.png', bg = "white", width = 7, height = 7, units = "in",
    pointsize = 12, res = 1000)
smoothplot(heightsmooth, heightscores)
dev.off()

png('LAsmooth.png', bg = "white", width = 7, height = 7, units = "in",
    pointsize = 12, res = 1000)
smoothplot(LAsmooth, LAscores)
dev.off()

png('SLAsmooth.png', bg = "white", width = 7, height = 7, units = "in",
    pointsize = 12, res = 1000)
smoothplot(SLAsmooth, SLAscores)
dev.off()

# soil depth variogram ----

soils <- lapply(soildepth, function(x) markvario(x, normalise = T, r=seq(from=0, to=25, length.out=512)))
plot.listof(soils, cbind(trans,theo) ~r, legend=F, xlim=c(0,25))

soilvario_CL <- pool(as.anylist(soils[1:5]), variance=T)
soilvario_OP <- pool(as.anylist(soils[c(6,8,9,10)]), variance=T)

png('soilvario.png', bg = "white", width = 7, height = 1.6, units = "in",
    pointsize = 8, res = 1000)
plot_soilvario <- plot(anylist(soilvario_CL,soilvario_OP), cbind(pooliso, pooltheo)~r, shade = c("loiso", "hiiso"), 
     xlim = c(0,25), legend = F, main="", main.panel = c("Closed", "Open"), equal.scales = T)
dev.off()

rm(soils,soilvario_CL,soilvario_OP)


# multi-trait with PCA ----------------------------------------------------

datPCA <- lapply(dat, na.omit)
datPCA <- purrr::map(datPCA,~.x[,c(-4,-6)])
dattt <- do.call(rbind, datPCA)

## log-transform ----------------------------------------------------------
dattt <- dattt |> 
  mutate(across(4:6, log10))

pca_results <- prcomp(dattt[,c(4:6)], scale. = T, rank. = 2)
pca_results$rotation[,2] <- -pca_results$rotation[,2]
pca_results$x[,2] <- -pca_results$x[,2]

autoplot(pca_results, loadings = T, loadings.label = T)

dattt <- cbind(dattt,pca_results$x)
dattt <- dattt |> mutate(plot=str_extract(rownames(dattt),".*?(?=\\.)")) |> 
  mutate(plot=as.numeric(as.factor(plot))) |> 
  mutate(Species=as.character(Species)) |> 
  mutate(Comm=ifelse(plot <= 5, 0, 1)) |> 
  mutate(Comm=factor(Comm, labels = c("Closed", "Open")))

#reassign to plots
principal_components <- pca_results$x
principal_components <- as.data.frame(principal_components)
num_rows_per_dataframe <- sapply(datPCA, function(df) nrow(df))
PCA_list <- split(principal_components, rep(seq_along(datPCA), num_rows_per_dataframe))

for (i in seq_along(datPCA)) {
  datPCA[[i]] <- cbind(datPCA[[i]], PCA_list[[i]])
}
rm(PCA_list,principal_components,num_rows_per_dataframe,dattt, pca_results, i)

#programita export ----

programheight <- lapply(height, as.data.frame) |> 
  lapply(function(x) select(x,-Species))

programLA <- lapply(LA, as.data.frame) |> 
  lapply(function(x) select(x,-Species))

programSLA <- lapply(SLA, as.data.frame) |> 
  lapply(function(x) select(x,-Species))

## height-----
programheight <- programheight |> 
  map(~ .x |>  mutate(pattern=1, dummy= 0) |> 
        relocate(pattern, .after=y))

for (i in 1:10){
  writeLines(paste("0", "50", "0", "50", nrow(programheight[[i]]), sep="\t"), paste0("height",i,".mcf"))
  write.table(programheight[[i]], paste0("height",i,".mcf"),col.names=FALSE, row.names=FALSE, append=TRUE)
}

## SLA-----
programSLA <- programSLA |> 
  map(~ .x |>  mutate(pattern=1, dummy= 0) |> 
        relocate(pattern, .after=y))

for (i in 1:10){
  writeLines(paste("0", "50", "0", "50", nrow(programSLA[[i]]), sep="\t"), paste0("SLA",i,".mcf"))
  write.table(programSLA[[i]], paste0("SLA",i,".mcf"),col.names=FALSE, row.names=FALSE, append=TRUE)
}

## LA-----
programLA <- programLA |> 
  map(~ .x |>  mutate(pattern=1, dummy= 0) |> 
        relocate(pattern, .after=y))

for (i in 1:10){
  writeLines(paste("0", "50", "0", "50", nrow(programLA[[i]]), sep="\t"), paste0("LA",i,".mcf"))
  write.table(programLA[[i]], paste0("LA",i,".mcf"),col.names=FALSE, row.names=FALSE, append=TRUE)
}

## PCA----
programPC1 <- datPCA |>
  map(~ .x |>  
        select(x,y, PC1) |> 
        mutate(pattern=1, dummy= 0) |> 
        relocate(pattern, .after=y))

for (i in 1:10){
  writeLines(paste("0", "50", "0", "50", nrow(programPC1[[i]]), sep="\t"), paste0("PC1_",i,".mcf"))
  write.table(programPC1[[i]], paste0("PC1_",i,".mcf"),col.names=FALSE, row.names=FALSE, append=TRUE)
}

programPC2 <- datPCA |>
  map(~ .x |>  
        select(x,y, PC2) |> 
        mutate(pattern=1, dummy= 0) |> 
        relocate(pattern, .after=y))

for (i in 1:10){
  writeLines(paste("0", "50", "0", "50", nrow(programPC2[[i]]), sep="\t"), paste0("PC2_",i,".mcf"))
  write.table(programPC2[[i]], paste0("PC2_",i,".mcf"),col.names=FALSE, row.names=FALSE, append=TRUE)
}

##only heterospecifics----

program_bivar2 <- function(x){
  for (k in seq_along(x)){
    df <- x[[k]]
  df$mark <- df[[4]]
  species <- unique(df$Species)
  for(i in seq_along(species)){
    output <- df |> 
      mutate(pattern = if_else(Species == species[i], 1, 2))|> 
      relocate(pattern, .after=y) |> 
      select(-Species)
    writeLines(paste("0", "50", "0", "50", nrow(df), sep="\t"), 
               paste0(getwd(),"/heterospecifics/",names(df[4]), names(x[k]), species[i],".mcf"))
    write.table(output, paste0(getwd(),"/heterospecifics/",names(df[4]), names(x[k]), species[i],".mcf"),
                col.names=FALSE, row.names=FALSE, append=TRUE)
  }
}
}
dir.create(file.path(getwd(),"heterospecifics"), recursive = TRUE)

###height----
program_bivar2(programheight[1:5])

programheight[[8]][81,3] <- "not determined"
program_bivar2(programheight[c(6,8,9,10)])
###SLA----
program_bivar2(programSLA[1:5])
program_bivar2(programSLA[c(6,8,9,10)])

###LA----
program_bivar2(programLA[1:5])
program_bivar2(programLA[c(6,8,9,10)])

###PCA----
programPC1 <- lapply(datPCA, function(x) select(x, x, y, Species, PC1))
programPC2 <- lapply(datPCA, function(x) select(x, x, y, Species, PC2))

program_bivar2(programPC1[1:5])
program_bivar2(programPC1[c(6,8,9,10)])

program_bivar2(programPC2[1:5])
program_bivar2(programPC2[c(6,8,9,10)])

#intensity for ring width----
lamb <- sapply(datpp, intensity)
lam_cl <- mean(lamb[1:5])
lam_op <- mean(lamb[6:10])
ring_cl <- 0.2/(sqrt(lam_cl))
ring_op <- 0.2/(sqrt(lam_op))

sink("ring.txt")
cat(lamb)
cat("lam_cl",lam_cl)
cat("lam_op",lam_op)
cat("ring_cl",ring_cl)
cat("ring_op",ring_op)
sink()
rm(lamb,lam_cl,lam_op)

# unmarked pair correlation-------------------------------

## closed ------------------------------------------------------------------

gr_closed <- anylapply(datpp[1:5], envelope, fun=pcf, nsim=999, nrank=25, fix.n=T, savefuns=T)
#isotropic edge corr

gr_closed_pool <- pool(gr_closed, savefuns=T)

bw <- max(sapply(datpp[1:5], bw.stoyan))

gr_closed_pool <-  crop_curves(gr_closed_pool, allfinite = T, r_min = bw)

plot_gr_closed <- plot(global_envelope_test(gr_closed_pool, alpha=0.05, type='area'))+
  ylab(expression(italic(g(r))))+
  xlab(expression(italic(r)(cm)))+
  theme(legend.position="none")

rm(gr_closed, gr_closed_pool, bw)

## open ------------------------------------------------------------------

gr_open <- anylapply(datpp[c(6,8,9,10)], envelope, fun=pcf, nsim=999, nrank=25, fix.n=T, savefuns=T)
#isotropic edge corr

gr_open_pool <- pool(gr_open, savefuns=T)

bw <- max(sapply(datpp[c(6,8,9,10)], bw.stoyan))

gr_open_pool <-  crop_curves(gr_open_pool, allfinite = T, r_min = bw)

plot_gr_open <- plot(global_envelope_test(gr_open_pool, alpha=0.05, type='area'))+
  ylab(expression(italic(g(r))))+
  xlab(expression(italic(r)(cm)))+
  theme(legend.position="none")

rm(gr_open, gr_open_pool, bw)

## patchwork ---------------------------------------------------------------

library(patchwork)
plot_gr_closed+plot_gr_open+
  plot_layout(ncol=2, axis_titles = "collect", guides = "collect")+
  plot_annotation(tag_levels = "a", tag_suffix = ")")&
  theme_minimal()&
  xlab(expression(italic(r)~(cm)))&
  theme(legend.position='none', plot.title = element_text(size=10),
        plot.tag = element_text(size = 10))
ggsave('gr.png', width = 7, height = 1.6, units = "in",
       bg = 'white', scale = 1.2, dpi = 1000)