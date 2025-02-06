#setup --------------------------------------------------------------------

library(GET)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)

rm(list=ls())

# functions ---------------------------------------------------------------

#from programita to GET
programita_names <- c("kmm_uni"="MCF11_t1", 
                      "vario_uni"="MCF11_t4", 
                      "cordens"="CorDens11",
                      "rmark2_uni"="MCF11_t3",
                      "rmark1_biv"="MCF12_t2",
                      "rmark2_biv"="MCF12_t3",
                      "kmm_biv"="MCF12_t1",
                      "vario_biv"="MCF12_t4")

make_envelope <- function(x, fun){
  if(missing(fun)) stop("summary function is required")
  fun <- unname(programita_names[fun])
  temp <- x |> 
    select(simnr, r, {{fun}}) |> 
    pivot_wider(names_from = simnr,values_from ={{fun}}, names_prefix = "rep") |> 
    data.matrix()
  r <- temp[,1]
  obs <- temp[,2]
  sim <- temp[,-c(1:2)]
  env <- curve_set(obs=obs,sim=sim, r=r)
  env
}

make_envelope_PCA <- function(x,y, fun){
  if(missing(fun)) stop("summary function is required")
  fun <- unname(programita_names[fun])
  PC1 <- x |> 
    select(simnr, r, {{fun}}) |> 
    pivot_wider(names_from = simnr,values_from ={{fun}}, names_prefix = "rep") |> 
    data.matrix()
  PC2 <- y |> 
    select(simnr, r, {{fun}}) |> 
    pivot_wider(names_from = simnr,values_from ={{fun}}, names_prefix = "rep") |> 
    data.matrix()
  temp <- PC1+PC2
  temp <- cbind(PC1[,1], temp[,-1])
  r <- temp[,1]
  obs <- temp[,2]
  sim <- temp[,-c(1:2)]
  env <- curve_set(obs=obs,sim=sim, r=r)
  env
}

  
#global envelope plot

plot_GE <- function(x){
  plot <- plot(global_envelope_test(x, alpha=0.05, type='area'))+
    scale_x_continuous(labels=function(x)x/10)
  plot$labels$title <- gsub("Global envelope test: ", "", plot$labels$title)
  plot+
    theme(title = element_text(face = if(as.numeric(gsub("p = ", "", plot$labels$title))>0.05)
    "plain"
    else "bold"))
}


#plot patchwork

patch <- function(x, ncol=2, pref="", tag_loc = "plot"){
  wrap_plots(x, ncol = ncol, byrow = F, guides = "collect")+
  plot_layout(axis_titles = "collect")+
  plot_annotation(tag_levels = "a", tag_suffix = ")", tag_prefix = pref)&
  theme_minimal()&
  xlab(expression(italic(r)~(cm)))&
  theme(legend.position='none', plot.title = element_text(size=10),
        plot.tag = element_text(size = 10), plot.tag.location = tag_loc)
}

#Height---------------------------------------------------------------------

Hclosed <- read.table("Univariate/H_CLOSED.rep", nrows=251000, sep="", fill=T, header=T, dec=".")
Hopen <- read.table("Univariate/H_OPEN_NO2.rep", nrows=251000, sep="", fill=T, header=T, dec=".")

##mark correlation----
#closed

Hclosedmc <- make_envelope(Hclosed, "kmm_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_Hclosedmc <- plot_GE(Hclosedmc)+
  ylab(expression(italic(k[mm](r))))

#open

Hopenmc <- make_envelope(Hopen, "kmm_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_Hopenmc <- plot_GE(Hopenmc)+
  ylab(expression(italic(k[mm](r))))

##variogram----
#closed

Hclosedvario <- make_envelope(Hclosed, "vario_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_Hclosedvario <- plot_GE(Hclosedvario)+
  ylab(expression(italic(gamma[m](r))))

#open

Hopenvario <- make_envelope(Hopen, "vario_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_Hopenvario <- plot_GE(Hopenvario)+
  ylab(expression(italic(gamma[m](r))))

##denscorr----
#closed

Hcloseddens <- make_envelope(Hclosed, "cordens") |> 
  crop_curves(r_min=0, r_max = 125)

plot_Hcloseddens <- plot_GE(Hcloseddens)+
  ylab(expression(italic(C[mK](r))))

#open

Hopendens <- make_envelope(Hopen, "cordens") |> 
  crop_curves(r_min=0, r_max = 125)

plot_Hopendens <- plot_GE(Hopendens)+
  ylab(expression(italic(C[mK](r))))

rm(Hclosed,Hopen, Hclosedmc, Hclosedvario, Hcloseddens,
   Hopenmc, Hopenvario, Hopendens)


#SLA---------------------------------------------------------------------------
SLAclosed <- read.table("Univariate/SLA_CLOSED.rep", nrows=251000, sep="", fill=T, header=T, dec=".")
SLAopen <- read.table("Univariate/SLA_OPEN_NO2.rep", nrows=251000, sep="", fill=T, header=T, dec=".")

##mark correlation----
#closed
SLAclosedmc <- make_envelope(SLAclosed, "kmm_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_SLAclosedmc <- plot_GE(SLAclosedmc)+
  ylab(expression(italic(k[mm](r))))

#open
SLAopenmc <- make_envelope(SLAopen, "kmm_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_SLAopenmc <- plot_GE(SLAopenmc)+
  ylab(expression(italic(k[mm](r))))

##variogram----
#closed
SLAclosedvario <- make_envelope(SLAclosed, "vario_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_SLAclosedvario <- plot_GE(SLAclosedvario)+
  ylab(expression(italic(gamma[m](r))))

#open
SLAopenvario <- make_envelope(SLAopen, "vario_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_SLAopenvario <- plot_GE(SLAopenvario)+
  ylab(expression(italic(gamma[m](r))))

##denscorr----
#closed
SLAcloseddens <- make_envelope(SLAclosed, "cordens") |> 
  crop_curves(r_min=0, r_max = 125)

plot_SLAcloseddens <- plot_GE(SLAcloseddens)+
  ylab(expression(italic(C[mK](r))))

#open
SLAopendens <- make_envelope(SLAopen, "cordens") |> 
  crop_curves(r_min=0, r_max = 125)

plot_SLAopendens <- plot_GE(SLAopendens)+
  ylab(expression(italic(C[mK](r))))

rm(SLAclosed,SLAopen, SLAclosedmc, SLAclosedvario, SLAcloseddens,
   SLAopenmc, SLAopenvario, SLAopendens)


#LA---------------------------------------------------------------------------
LAclosed <- read.table("Univariate/LA_CLOSED.rep", nrows=251000, sep="", fill=T, header=T, dec=".")
LAopen <- read.table("Univariate/LA_OPEN_NO2.rep", nrows=251000, sep="", fill=T, header=T, dec=".")

##mark correlation----
#closed

LAclosedmc <- make_envelope(LAclosed, "kmm_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_LAclosedmc <- plot_GE(LAclosedmc)+
  ylab(expression(italic(k[mm](r))))

#open

LAopenmc <- make_envelope(LAopen, "kmm_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_LAopenmc <- plot_GE(LAopenmc)+
  ylab(expression(italic(k[mm](r))))

##variogram----
#closed

LAclosedvario <- make_envelope(LAclosed, "vario_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_LAclosedvario <- plot_GE(LAclosedvario)+
  ylab(expression(italic(gamma[m](r))))

#open

LAopenvario <- make_envelope(LAopen, "vario_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_LAopenvario <- plot_GE(LAopenvario)+
  ylab(expression(italic(gamma[m](r))))

##denscorr----
#closed

LAcloseddens <- make_envelope(LAclosed, "cordens") |> 
  crop_curves(r_min=0, r_max = 125)

plot_LAcloseddens <- plot_GE(LAcloseddens)+
  ylab(expression(italic(C[mK](r))))

#open

LAopendens <- make_envelope(LAopen, "cordens") |> 
  crop_curves(r_min=0, r_max = 125)

plot_LAopendens <- plot_GE(LAopendens)+
  ylab(expression(italic(C[mK](r))))

rm(LAclosed,LAopen, LAclosedmc, LAclosedvario, LAcloseddens,
   LAopenmc, LAopenvario, LAopendens)


#PCA ----
##closed----

PC1_closed<- read.table("Univariate/PC1_CLOSED_LOG.rep", nrows=126000, sep="", fill=T, header=T, dec=".")
PC2_closed <- read.table("Univariate/PC2_CLOSED_LOG.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

PCAclosedvario <- make_envelope_PCA(PC1_closed, PC2_closed, fun = "vario_uni")

plot_PCAclosedvario <- plot_GE(PCAclosedvario)+
  ylab(expression(italic(gamma[m](r))))

rm(PC1_closed,PC2_closed, PCAclosedvario)

##open----

PC1_open<- read.table("Univariate/PC1_OPEN_LOG.rep", nrows=126000, sep="", fill=T, header=T, dec=".")
PC2_open <- read.table("Univariate/PC2_OPEN_LOG.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

PCAopenvario <- make_envelope_PCA(PC1_open, PC2_open, fun = "vario_uni")

plot_PCAopenvario <- plot_GE(PCAopenvario)+
  ylab(expression(italic(gamma[m](r))))

rm(PC1_open, PC2_open, PCAopenvario)


#Heterospecific --------------------------------------------------------------------

##height closed---------------------------------------------------------------------

Hclosedhet <- read.table("hetero/mcf_CLOSEDheight.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation----

Hclosedmc_het <- make_envelope(Hclosedhet, fun = "kmm_biv")

plot_Hclosedmc_het <- plot_GE(Hclosedmc_het)+
  ylab(expression(italic(k[mm](r))))

###mark variogram----

Hclosedvario_het <- make_envelope(Hclosedhet, fun = "vario_biv")

plot_Hclosedvario_het <- plot_GE(Hclosedvario_het)+
  ylab(expression(italic(gamma[m](r))))

###rmark----

Hclosedrm_het <- make_envelope(Hclosedhet, fun = "rmark2_biv")

plot_Hclosedrm_het <- plot_GE(Hclosedrm_het)+
  ylab(expression(italic(k[.m](r))))

rm(Hclosedhet, Hclosedmc_het, Hclosedvario_het, Hclosedrm_het)

##height open---------------------------------------------------------------------

Hopenhet <- read.table("hetero/mcf_OPENheight.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation----

Hopenmc_het <- make_envelope(Hopenhet, fun = "kmm_biv")

plot_Hopenmc_het <- plot_GE(Hopenmc_het)+
  ylab(expression(italic(k[mm](r))))

###mark variogram----

Hopenvario_het <- make_envelope(Hopenhet, fun = "vario_biv")

plot_Hopenvario_het <- plot_GE(Hopenvario_het)+
  ylab(expression(italic(gamma[m](r))))

###rmark----

Hopenrm_het <- make_envelope(Hopenhet, fun = "rmark2_biv")

plot_Hopenrm_het <- plot_GE(Hopenrm_het)+
  ylab(expression(italic(k[.m](r))))

rm(Hopenhet, Hopenmc_het, Hopenvario_het, Hopenrm_het)


##LA closed---------------------------------------------------------------------

LAclosedhet <- read.table("hetero/mcf_CLOSEDla.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation----

LAclosedmc_het <- make_envelope(LAclosedhet, fun = "kmm_biv")

plot_LAclosedmc_het <- plot_GE(LAclosedmc_het)+
  ylab(expression(italic(k[mm](r))))

###mark variogram----

LAclosedvario_het <- make_envelope(LAclosedhet, fun = "vario_biv")

plot_LAclosedvario_het <- plot_GE(LAclosedvario_het)+
  ylab(expression(italic(gamma[m](r))))

###rmark----

LAclosedrm_het <- make_envelope(LAclosedhet, fun = "rmark2_biv")

plot_LAclosedrm_het <- plot_GE(LAclosedrm_het)+
  ylab(expression(italic(k[.m](r))))

rm(LAclosedhet, LAclosedmc_het, LAclosedvario_het, LAclosedrm_het)

##LA open---------------------------------------------------------------------

LAopenhet <- read.table("hetero/mcf_OPENla.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation----

LAopenmc_het <- make_envelope(LAopenhet, fun = "kmm_biv")

plot_LAopenmc_het <- plot_GE(LAopenmc_het)+
  ylab(expression(italic(k[mm](r))))

###mark variogram----

LAopenvario_het <- make_envelope(LAopenhet, fun = "vario_biv")

plot_LAopenvario_het <- plot_GE(LAopenvario_het)+
  ylab(expression(italic(gamma[m](r))))

###rmark----

LAopenrm_het <- make_envelope(LAopenhet, fun = "rmark2_biv")

plot_LAopenrm_het <- plot_GE(LAopenrm_het)+
  ylab(expression(italic(k[.m](r))))

rm(LAopenhet, LAopenmc_het, LAopenvario_het, LAopenrm_het)


##SLA closed---------------------------------------------------------------------

SLAclosedhet <- read.table("hetero/mcf_CLOSEDsla.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation----

SLAclosedmc_het <- make_envelope(SLAclosedhet, fun = "kmm_biv")

plot_SLAclosedmc_het <- plot_GE(SLAclosedmc_het)+
  ylab(expression(italic(k[mm](r))))

###mark variogram----

SLAclosedvario_het <- make_envelope(SLAclosedhet, fun = "vario_biv")

plot_SLAclosedvario_het <- plot_GE(SLAclosedvario_het)+
  ylab(expression(italic(gamma[m](r))))

###rmark----

SLAclosedrm_het <- make_envelope(SLAclosedhet, fun = "rmark2_biv")

plot_SLAclosedrm_het <- plot_GE(SLAclosedrm_het)+
  ylab(expression(italic(k[.m](r))))

rm(SLAclosedhet, SLAclosedmc_het, SLAclosedvario_het, SLAclosedrm_het)

##SLA open---------------------------------------------------------------------

SLAopenhet <- read.table("hetero/mcf_OPENsla.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation----

SLAopenmc_het <- make_envelope(SLAopenhet, fun = "kmm_biv")

plot_SLAopenmc_het <- plot_GE(SLAopenmc_het)+
  ylab(expression(italic(k[mm](r))))

###mark variogram----

SLAopenvario_het <- make_envelope(SLAopenhet, fun = "vario_biv")

plot_SLAopenvario_het <- plot_GE(SLAopenvario_het)+
  ylab(expression(italic(gamma[m](r))))

###rmark----

SLAopenrm_het <- make_envelope(SLAopenhet, fun = "rmark2_biv")

plot_SLAopenrm_het <- plot_GE(SLAopenrm_het)+
  ylab(expression(italic(k[.m](r))))

rm(SLAopenhet, SLAopenmc_het, SLAopenvario_het, SLAopenrm_het)

##PCA ----

###closed----

PC1_closed_het <- read.table("hetero/mcf_CLOSEDpc1_log.rep", nrows=126000, sep="", fill=T, header=T, dec=".")
PC2_closed_het <- read.table("hetero/mcf_CLOSEDpc2_log.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

PCAclosedvario_het <- make_envelope_PCA(PC1_closed_het, PC2_closed_het, fun = "vario_biv")

plot_PCAclosedvario_het <- plot_GE(PCAclosedvario_het)+
  ylab(expression(italic(gamma[m](r))))

rm(PC1_closed_het, PC2_closed_het, PCAclosedvario_het)

###open----

PC1_open_het <- read.table("hetero/mcf_OPENpc1_log.rep", nrows=126000, sep="", fill=T, header=T, dec=".")
PC2_open_het <- read.table("hetero/mcf_OPENpc2_log.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

PCAopenvario_het <- make_envelope_PCA(PC1_open_het, PC2_open_het, fun = "vario_biv")

plot_PCAopenvario_het <- plot_GE(PCAopenvario_het)+
  ylab(expression(italic(gamma[m](r))))

rm(PC1_open_het, PC2_open_het, PCAopenvario_het)


# patchworks ------

## height ------------------------------------------------------------------

H_closed <- list(plot_Hclosedmc, plot_Hclosedvario, plot_Hcloseddens,
               plot_Hclosedmc_het, plot_Hclosedvario_het, plot_Hclosedrm_het)

H_closed <- patch(H_closed, pref = "I")+
    plot_annotation(title = "Closed",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                                            margin = margin(0,0,0,0, unit = "in"))))

H_open <- list(plot_Hopenmc, plot_Hopenvario, plot_Hopendens,
                 plot_Hopenmc_het, plot_Hopenvario_het, plot_Hopenrm_het)

H_open <- patch(H_open, pref = "II")+
  plot_annotation(title = "Open",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                                          margin = margin(0,0,0,0, unit = "in"))))

wrap_elements(H_closed) + plot_spacer() + wrap_elements(H_open)+
  plot_layout(ncol = 3, widths = c(10, -1.01, 10))

ggsave('Plots/height.png', width = 7, height = 3.5, units = "in",
       bg = 'white', scale = 1.5, dpi = 1000)

ggsave('Plots/Height.eps', width = 7, height = 3.5, units = "in", bg = 'white',
       scale = 1.5)

## LA ------------------------------------------------------------------

LA_closed <- list(plot_LAclosedmc, plot_LAclosedvario, plot_LAcloseddens,
                 plot_LAclosedmc_het, plot_LAclosedvario_het, plot_LAclosedrm_het)

LA_closed <- patch(LA_closed, pref = "I")+
  plot_annotation(title = "Closed",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                                          margin = margin(0,0,0,0, unit = "in"))))

LA_open <- list(plot_LAopenmc, plot_LAopenvario, plot_LAopendens,
               plot_LAopenmc_het, plot_LAopenvario_het, plot_LAopenrm_het)

LA_open <- patch(LA_open, pref = "II")+
  plot_annotation(title = "Open",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                                          margin = margin(0,0,0,0, unit = "in"))))

wrap_elements(LA_closed) + plot_spacer() + wrap_elements(LA_open)+
  plot_layout(ncol = 3, widths = c(10, -1.01, 10))

ggsave('Plots/LA.png', width = 7, height = 3.5, units = "in",
       bg = 'white', scale = 1.5, dpi = 1000)

ggsave('Plots/LA.eps', width = 7, height = 3.5, units = "in", bg = 'white',
       scale = 1.5)

## SLA ------------------------------------------------------------------

SLA_closed <- list(plot_SLAclosedmc, plot_SLAclosedvario, plot_SLAcloseddens,
                  plot_SLAclosedmc_het, plot_SLAclosedvario_het, plot_SLAclosedrm_het)

SLA_closed <- patch(SLA_closed, pref = "I")+
  plot_annotation(title = "Closed",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                                          margin = margin(0,0,0,0, unit = "in"))))

SLA_open <- list(plot_SLAopenmc, plot_SLAopenvario, plot_SLAopendens,
                plot_SLAopenmc_het, plot_SLAopenvario_het, plot_SLAopenrm_het)

SLA_open <- patch(SLA_open, pref = "II")+
  plot_annotation(title = "Open",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                                          margin = margin(0,0,0,0, unit = "in"))))

wrap_elements(SLA_closed) + plot_spacer() + wrap_elements(SLA_open)+
  plot_layout(ncol = 3, widths = c(10, -1.01, 10))

ggsave('Plots/SLA.png', width = 7, height = 3.5, units = "in",
       bg = 'white', scale = 1.5, dpi = 1000)

ggsave('Plots/SLA.eps', width = 7, height = 3.5, units = "in", bg = 'white',
       scale = 1.5)

## PCA ------------------------------------------------------------------

PCA_closed <- list(plot_PCAclosedvario, plot_PCAclosedvario_het)

PCA_closed <- patch(PCA_closed, ncol = 1, pref = "I", tag_loc = "margin")+
  plot_annotation(title = "Closed",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                                          margin = margin(0,0,0,0, unit = "in"))))

PCA_open <- list(plot_PCAopenvario, plot_PCAopenvario_het)

PCA_open <- patch(PCA_open, ncol = 1, pref = "II", tag_loc = "margin")+
  plot_annotation(title = "Closed",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                                          margin = margin(0,0,0,0, unit = "in"))))

wrap_elements(PCA_closed) + plot_spacer() + wrap_elements(PCA_open)+
  plot_layout(ncol = 3, widths = c(10, -1.01, 10))

ggsave('Plots/PCA.png', width = 7, height = 3.5, units = "in",
       bg = 'white', scale = 1.5, dpi = 1000)

ggsave('Plots/PCA.eps', width = 7, height = 3.5, units = "in", bg = 'white',
       scale = 1.5)

# #for slides
# patch_nobg(H_list)
# ggsave('Plots/Height_nobg.png', bg='transparent', width = 29.21, height = 14, units = "cm",
#        dpi = "print", scale = 0.8)

#HS---------------------------------------------------------------------------

LSclosed <- read.table("Univariate/LS_CLOSED.rep", nrows=126000, sep="", fill=T, header=T, dec=".")
LSopen <- read.table("Univariate/LS_OPEN_NO2.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

##denscorr----
LScloseddens <- make_envelope(LSclosed, "cordens")|> 
  crop_curves(r_min=0, r_max = 124) #last bin is correlation with NN

plot_LScloseddens <- plot_GE(LScloseddens)+
  ylab(expression(italic(C[mK](r))))

LSopendens <- make_envelope(LSopen, "cordens")|> 
  crop_curves(r_min=0, r_max = 124) #last bin is correlation with NN

plot_LSopendens <- plot_GE(LSopendens)+
  ylab(expression(italic(C[mK](r))))

rm(LSclosed,LSopen)

patch(list(plot_LScloseddens,plot_LSopendens), ncol = 2)

ggsave('Plots/LS.png', width = 7, height = 1.6, units = "in",
       bg = 'white', scale = 1.2, dpi = 1000)

#Local 15cm -----

##height closed---------------------------------------------------------------------

Hclosedloc <- read.table("local/H_CLOSED_loc.rep", nrows=126000, sep="", fill=T, header=T, dec=".")


###mark correlation----
Hclosedmc_loc <- make_envelope(Hclosedloc, fun = "kmm_uni")
  
plot_Hclosedmark_loc <- plot_GE(Hclosedmc_loc)+
  ylab(expression(italic(k[mm](r))))

###mark variogram----

Hclosedvario_loc <- make_envelope(Hclosedloc, fun = "vario_uni")

plot_Hclosedvario_loc <- plot_GE(Hclosedvario_loc)+
  ylab(expression(italic(gamma[m](r))))

###denscorr----

Hcloseddens_loc <- make_envelope(Hclosedloc, fun = "cordens") |> 
  crop_curves(r_min=0, r_max = 124)
  
plot_Hcloseddens_loc <- plot_GE(Hcloseddens_loc)+
  ylab(expression(italic(C[mK](r))))

rm(Hclosedloc, Hclosedmc_loc,Hclosedvario_loc, Hcloseddens_loc)

##LA closed---------------------------------------------------------------------

LAclosedloc <- read.table("local/LA_CLOSED_loc.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation----

LAclosedmc_loc <- make_envelope(LAclosedloc, fun = "kmm_uni")

plot_LAclosedmark_loc <- plot_GE(LAclosedmc_loc)+
  ylab(expression(italic(k[mm](r))))

###mark variogram----

LAclosedvario_loc <- make_envelope(LAclosedloc, fun = "vario_uni")

plot_LAclosedvario_loc <- plot_GE(LAclosedvario_loc)+
  ylab(expression(italic(gamma[m](r))))

###denscorr----

LAcloseddens_loc <- make_envelope(LAclosedloc, fun = "cordens") |> 
  crop_curves(r_min=0, r_max = 124)

plot_LAcloseddens_loc <- plot_GE(LAcloseddens_loc)+
  ylab(expression(italic(C[mK](r))))

rm(LAclosedloc, LAclosedmc_loc, LAclosedvario_loc, LAcloseddens_loc)

##PCA-----

###open----

PC1_open_loc <- read.table("local/PC1_OPEN_LOC.rep", nrows=126000, sep="", fill=T, header=T, dec=".")
PC2_open_loc <- read.table("local/PC2_OPEN_LOC.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

PCAopenvario_loc <- make_envelope_PCA(PC1_open_loc, PC2_open_loc, fun = "vario_uni")

plot_PCAopenvario_loc <- plot_GE(PCAopenvario_loc)+
  ylab(expression(italic(gamma[m](r))))

rm(PC1_open_loc, PC2_open_loc, PCAopenvario_loc)


#Local hetero 15cm -----------------------------------------------------------------

##height closed---------------------------------------------------------------------

Hclosedlochet <- read.table("local/H_CLOSED_HET_LOC.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation

Hclosedmc_het_loc <- make_envelope(Hclosedlochet, fun = "kmm_biv")

plot_Hclosedmc_het_loc <- plot_GE(Hclosedmc_het_loc)+
  ylab(expression(italic(k[mm](r))))

###mark variogram----

Hclosedvario_het_loc <- make_envelope(Hclosedlochet, fun = "vario_biv")

plot_Hclosedvario_het_loc <- plot_GE(Hclosedvario_het_loc)+
  ylab(expression(italic(gamma[m](r))))

###rmark----

Hclosedrm_het_loc <- make_envelope(Hclosedlochet, fun = "rmark2_biv")  

plot_Hclosedrm_het_loc <- plot_GE(Hclosedrm_het_loc)+
  ylab(expression(italic(k[.m](r))))


rm(Hclosedlochet, Hclosedmc_het_loc, Hclosedvario_het_loc, Hclosedrm_het_loc)

##LA closed---------------------------------------------------------------------

LAclosedlochet <- read.table("local/LA_CLOSED_HET_LOC.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation

LAclosedmc_het_loc <- make_envelope(LAclosedlochet, fun = "kmm_biv")

plot_LAclosedmc_het_loc <- plot_GE(LAclosedmc_het_loc)+
  ylab(expression(italic(k[mm](r))))

###mark variogram----

LAclosedvario_het_loc <- make_envelope(LAclosedlochet, fun = "vario_biv")

plot_LAclosedvario_het_loc <- plot_GE(LAclosedvario_het_loc)+
  ylab(expression(italic(gamma[m](r))))

###rmark----

LAclosedrm_het_loc <- make_envelope(LAclosedlochet, fun = "rmark2_biv")  

plot_LAclosedrm_het_loc <- plot_GE(LAclosedrm_het_loc)+
  ylab(expression(italic(k[.m](r))))


rm(LAclosedlochet, LAclosedmc_het_loc, LAclosedvario_het_loc, LAclosedrm_het_loc)

##SLA closed---------------------------------------------------------------------

SLAclosedlochet <- read.table("local/SLA_CLOSED_HET_LOC.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation

SLAclosedmc_het_loc <- make_envelope(SLAclosedlochet, fun = "kmm_biv")

plot_SLAclosedmc_het_loc <- plot_GE(SLAclosedmc_het_loc)+
  ylab(expression(italic(k[mm](r))))

###mark variogram----

SLAclosedvario_het_loc <- make_envelope(SLAclosedlochet, fun = "vario_biv")

plot_SLAclosedvario_het_loc <- plot_GE(SLAclosedvario_het_loc)+
  ylab(expression(italic(gamma[m](r))))

###rmark----

SLAclosedrm_het_loc <- make_envelope(SLAclosedlochet, fun = "rmark2_biv")  

plot_SLAclosedrm_het_loc <- plot_GE(SLAclosedrm_het_loc)+
  ylab(expression(italic(k[.m](r))))


rm(SLAclosedlochet, SLAclosedmc_het_loc, SLAclosedvario_het_loc, SLAclosedrm_het_loc)

##SLA open---------------------------------------------------------------------

SLAopenlochet <- read.table("local/SLA_OPEN_HET_LOC.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation

SLAopenmc_het_loc <- make_envelope(SLAopenlochet, fun = "kmm_biv")

plot_SLAopenmc_het_loc <- plot_GE(SLAopenmc_het_loc)+
  ylab(expression(italic(k[mm](r))))

###mark variogram----

SLAopenvario_het_loc <- make_envelope(SLAopenlochet, fun = "vario_biv")

plot_SLAopenvario_het_loc <- plot_GE(SLAopenvario_het_loc)+
  ylab(expression(italic(gamma[m](r))))

###rmark----

SLAopenrm_het_loc <- make_envelope(SLAopenlochet, fun = "rmark2_biv")  

plot_SLAopenrm_het_loc <- plot_GE(SLAopenrm_het_loc)+
  ylab(expression(italic(k[.m](r))))


rm(SLAopenlochet, SLAopenmc_het_loc, SLAopenvario_het_loc, SLAopenrm_het_loc)

##PCA-----

###closed----

PC1_closed_het_loc <- read.table("local/PC1_CLOSED_HET_LOC.rep", nrows=126000, sep="", fill=T, header=T, dec=".")
PC2_closed_het_loc <- read.table("local/PC2_CLOSED_HET_LOC.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

PCAclosedvario_het_loc <- make_envelope_PCA(PC1_closed_het_loc, PC2_closed_het_loc, fun = "vario_biv")

plot_PCAclosedvario_het_loc <- plot_GE(PCAclosedvario_het_loc)+
  ylab(expression(italic(gamma[m](r))))

rm(PC1_closed_het_loc, PC2_closed_het_loc, PCAclosedvario_het_loc)

# local patchworks --------------------------------------------------------
## Height -----
H_closed_loc <- list(plot_Hclosedmark_loc, plot_Hclosedvario_loc, plot_Hcloseddens_loc,
                 plot_Hclosedmc_het_loc, plot_Hclosedvario_het_loc, plot_Hclosedrm_het_loc)

H_closed_loc <- patch(H_closed_loc)
H_closed_loc
ggsave('Plots/Height_loc.png', width = 7, height = 3.5, units = "in",
       bg = 'white', scale = 1.2, dpi = 1000)


## LA-----

LA_closed_loc <- list(plot_LAclosedmark_loc, plot_LAclosedvario_loc, plot_LAcloseddens_loc,
                     plot_LAclosedmc_het_loc, plot_LAclosedvario_het_loc, plot_LAclosedrm_het_loc)

LA_closed_loc <- patch(LA_closed_loc)
LA_closed_loc
ggsave('Plots/LA_loc.png', width = 7, height = 3.5, units = "in",
       bg = 'white', scale = 1.2, dpi = 1000)

## SLA----

SLA_loc <- list(plot_SLAclosedmc_het_loc, plot_SLAclosedvario_het_loc, plot_SLAclosedrm_het_loc,
                plot_SLAopenmc_het_loc, plot_SLAopenvario_het_loc, plot_SLAopenrm_het_loc)

SLA_loc <- patch(SLA_loc)
SLA_loc
ggsave('Plots/SLA_loc.png', width = 7, height = 3.5, units = "in",
       bg = 'white', scale = 1.2, dpi = 1000)

## PCA----

PCA_loc <- patch(list(plot_PCAclosedvario_het_loc, plot_PCAopenvario_loc), ncol = 1, tag_loc = "margin")
PCA_loc
ggsave('Plots/PCA_loc.png', width = 7, height = 3.5, units = "in",
       bg = 'white', scale = 1.2, dpi = 1000)


# species means -----------------------------------------------------------

## Height closed-----------------------------------------------------------

Hclosedmean <- read.table("Mean/H_CLOSEDmean.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

##mark correlation

Hclosedmcmean <- make_envelope(Hclosedmean, "kmm_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_Hclosedmcmean <- plot_GE(Hclosedmcmean)+
  ylab(expression(italic(k[mm](r))))


##variogram

Hclosedvariomean <- make_envelope(Hclosedmean, "vario_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_Hclosedvariomean <- plot_GE(Hclosedvariomean)+
  ylab(expression(italic(gamma[m](r))))

##denscorr

Hcloseddensmean <- make_envelope(Hclosedmean, "cordens") |> 
  crop_curves(r_min=0, r_max = 124)

plot_Hcloseddensmean <- plot_GE(Hcloseddensmean)+
  ylab(expression(italic(C[mK](r))))


### heterospecific -------------------------------------------------

Hclosedhetmean <- read.table("Mean/H_CLOSEDmeanhet.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation

Hclosedmc_hetmean <- make_envelope(Hclosedhetmean, fun = "kmm_biv")

plot_Hclosedmc_hetmean <- plot_GE(Hclosedmc_hetmean)+
  ylab(expression(italic(k[mm](r))))

###mark variogram

Hclosedvario_hetmean <- make_envelope(Hclosedhetmean, fun = "vario_biv")

plot_Hclosedvario_hetmean <- plot_GE(Hclosedvario_hetmean)+
  ylab(expression(italic(gamma[m](r))))

###rmark

Hclosedrm_hetmean <- make_envelope(Hclosedhetmean, fun = "rmark2_biv")

plot_Hclosedrm_hetmean <- plot_GE(Hclosedrm_hetmean)+
  ylab(expression(italic(k[.m](r))))


H_closedmean <- list(plot_Hclosedmcmean, plot_Hclosedvariomean, plot_Hcloseddensmean,
                 plot_Hclosedmc_hetmean, plot_Hclosedvario_hetmean, plot_Hclosedrm_hetmean)

patch(H_closedmean)+
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                                          margin = margin(0,0,0,0, unit = "in"))))



ggsave('Plots/heightmean.png', width = 7, height = 3.5, units = "in",
       bg = 'white', scale = 1.5, dpi = 1000)

## LA closed-----------------------------------------------------------

LAclosedmean <- read.table("Mean/LA_CLOSEDmean.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

##mark correlation

LAclosedmcmean <- make_envelope(LAclosedmean, "kmm_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_LAclosedmcmean <- plot_GE(LAclosedmcmean)+
  ylab(expression(italic(k[mm](r))))


##variogram

LAclosedvariomean <- make_envelope(LAclosedmean, "vario_uni") |> 
  crop_curves(r_min=0, r_max = 125)

plot_LAclosedvariomean <- plot_GE(LAclosedvariomean)+
  ylab(expression(italic(gamma[m](r))))

##denscorr

LAcloseddensmean <- make_envelope(LAclosedmean, "cordens") |> 
  crop_curves(r_min=0, r_max = 124)

plot_LAcloseddensmean <- plot_GE(LAcloseddensmean)+
  ylab(expression(italic(C[mK](r))))


### heterospecific -------------------------------------------------

LAclosedhetmean <- read.table("Mean/LA_CLOSEDmeanhet.rep", nrows=126000, sep="", fill=T, header=T, dec=".")

###mark correlation

LAclosedmc_hetmean <- make_envelope(LAclosedhetmean, fun = "kmm_biv")

plot_LAclosedmc_hetmean <- plot_GE(LAclosedmc_hetmean)+
  ylab(expression(italic(k[mm](r))))

###mark variogram

LAclosedvario_hetmean <- make_envelope(LAclosedhetmean, fun = "vario_biv")

plot_LAclosedvario_hetmean <- plot_GE(LAclosedvario_hetmean)+
  ylab(expression(italic(gamma[m](r))))

###rmark

LAclosedrm_hetmean <- make_envelope(LAclosedhetmean, fun = "rmark2_biv")

plot_LAclosedrm_hetmean <- plot_GE(LAclosedrm_hetmean)+
  ylab(expression(italic(k[.m](r))))


LA_closedmean <- list(plot_LAclosedmcmean, plot_LAclosedvariomean, plot_LAcloseddensmean,
                     plot_LAclosedmc_hetmean, plot_LAclosedvario_hetmean, plot_LAclosedrm_hetmean)

patch(LA_closedmean)+
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                                          margin = margin(0,0,0,0, unit = "in"))))



ggsave('Plots/LAmean.png', width = 7, height = 3.5, units = "in",
       bg = 'white', scale = 1.5, dpi = 1000)
