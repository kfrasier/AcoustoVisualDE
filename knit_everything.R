setwd('F:/NASData/AcoustoVisualDE/AcoustoVisualDE')
#
# myList =  c("Pm_model_GAM_density_GRIIDC.Rmd","Zc_model_GAM_density_GRIIDC.Rmd","Gg_model_GAM_density_GRIIDC.Rmd",
#           ,"Zc_model_NN_density_GRIIDC.Rmd","Kspp_model_NN_density_GRIIDC.Rmd",
#"Ssp_model_NN_density_GRIIDC.Rmd","Me_model_NN_density_GRIIDC.Rmd", "Pm_model_NN_density_GRIIDC.Rmd","Gg_model_NNet_density_GRIIDC.Rmd",
 myList =  c(           "Gmsp_model_NN_density_GRIIDC.Rmd")

#"Pm_model_runs_NN_density.Rmd","Zc_model_runs_NN_density.Rmd","Gg_model_runs_NN_density.Rmd",
#"Pm_model_runs_NN.Rmd","Gg_model_runs_NN.Rmd","Zc_model_runs_NN.Rmd",
# myList = c("Gg_model_runs_NN_density.Rmd","Pm_model_runs_NN_density.Rmd","Pm_model_runs_gam.Rmd", "Zc_model_runs_gam.Rmd", "Gg_model_runs_gam.Rmd")
#            "Pm_model_runs_NN.Rmd","Gg_model_runs_NN.Rmd","Gmsp_model_runs_NN_density.Rmd")"Gg_model_runs_gam_density.Rmd",
#          "Zc_model_runs_NN.Rmd",
# myList =  c("Gg_model_runs_gam_density.Rmd",
#             "Zc_model_runs_gam_density.Rmd",
#             "Ssp_model_runs_gam_density.Rmd", "Pm_model_runs_gam_density.Rmd",
#             "Kspp_model_runs_gam_density.Rmd","Gmsp_model_runs_gam_density.Rmd","Me_model_runs_gam_density.Rmd")
# 
# myList =  c("Zc_model_runs_NN.Rmd",
#              "Ssp_model_runs_NN.Rmd", "Pm_model_runs_gam.Rmd",
#              "Kspp_model_runs_NN.Rmd", "Me_model_runs_NN.Rmd")#"Gg_model_runs_NN_density.Rmd", "Gg_model_runs_NN.Rmd",
#               # don't run because of NAs?
# myList =  c("Gg_model_runs_gam.Rmd","Zc_model_runs_gam.Rmd",
#             "Ssp_model_runs_gam.Rmd","Pm_model_runs_gam.Rmd",
#             "Kspp_model_runs_gam.Rmd","Gmsp_model_runs_gam.Rmd")# no sightings:"Me_model_runs_gam.Rmd"

# "Gmsp_model_runs_NN_density.Rmd", 
#              "Zc_model_runs_NN_density.Rmd", "Me_model_runs_NN_density.Rmd",
#              "Ssp_model_runs_NN_density.Rmd", "Pm_model_runs_gam_density.Rmd",
#              "Kspp_model_runs_NN_density.Rmd",
# myList =  c( 
#              "Zc_model_runs_NN.Rmd", 
#              "Ssp_model_runs_NN.Rmd", "Pm_model_runs_NN.Rmd", "Pm_model_runs_NN_density.Rmd", 
#              "Kspp_model_runs_NN.Rmd")#"Me_model_runs_NN_density.Rmd","Me_model_runs_NN.Rmd","Gmsp_model_runs_NN.Rmd", "Gg_model_runs_NN.Rmd",
#    "Kspp_model_runs_gam_density.Rmd",
#   "Kspp_model_runs_NN.Rmd",  "Kspp_model_runs_NN_density.Rmd",  
#   "Me_model_runs_NN.Rmd",    "Me_model_runs_NN_density.Rmd",
#   "Pm_model_runs_gam.Rmd",   "Pm_model_runs_gam_density.Rmd",
#   "Pm_model_runs_NN.Rmd",    "Pm_model_runs_NN_density.Rmd",     
#   "Zc_model_runs_gam.Rmd",   "Zc_model_runs_gam_density.Rmd",        
#   "Zc_model_runs_NN.Rmd",    "Zc_model_runs_NNdensity.Rmd",
#  "Gg_model_runs_gam.Rmd",  
#  "Gg_model_runs_gam_density.Rmd",        
#  "Gg_model_runs_NN.Rmd",    
#                             "Gg_model_runs_NN_density.Rmd",
#  "Ssp_model_runs_gam.Rmd",  
#  "Ssp_model_runs_gam_density.Rmd",        
#  "Ssp_model_runs_NN.Rmd",  
#  "Ssp_model_runs_NN_density.Rmd",
#  "Gmsp_model_runs_gam.Rmd",
# "Gmsp_model_runs_gam_density.Rmd",        
#  "Gmsp_model_runs_NN.Rmd", 
# "Gmsp_model_runs_NN_density.Rmd")

for (f in myList) rmarkdown::render(f)
