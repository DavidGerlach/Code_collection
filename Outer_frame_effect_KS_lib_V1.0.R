# This function analyses and corrects outerframe effects for KS library

#Version 1.0: development



#======= Import libraries =======

library(ggplot2)
library(stringr)
library(RColorBrewer)
library(scales) #for simpsons color palette
library(ggsci) #for Simpsons color palette
library(beeswarm) #for acurate dotplots
library(dplyr)
library(plyr)
library(prodlim)
library(raincloudplots)
library(plot.matrix) #for easy heatmaps
library(vioplot)

#======== Set paths ==========
here_path <- "C:/Users/xuehua/Documents/Wuerzburg_temp/Students/Kristine_Schweinshaut/KS_lib/"
Out_dir <- paste0(here_path,"Output/")
In_dir <- paste0(here_path,"Input/")
source_path <- "C:/Users/xuehua/Documents/Wuerzburg_temp/Scripts/"

#======= Source functions ======

source(paste0(source_path,"Auxs/","Growth_analyzer_func_V2.0.R",collapse=NULL)) #source standard growth analyze
source(paste0(source_path,"Auxs/","t_lag_function_V6.2.R",collapse=NULL)) # use latest version
source(paste0(source_path,"Auxs/","Aux_functions_DG_V1.2.R",collapse=NULL))
source(paste0(source_path,"Auxs/","Aux_functions_2.R" ,collapse=NULL))
source(paste0(source_path,"Auxs/","SAMI_convert_DG1.R")) #call ARB's growth converter | function Ext_SAMIplates
source(paste0(source_path,"Auxs/","t_lag_function_sEc_lib_V1.0.R"))
source(paste0(source_path,"Auxs/","SAMI_organize_V1.0.R")) #master function to organize SAMI data

####==== Parameters #######

current_version <- "1.0"

#===== Load Data ====

#load KS IRIS + meta table
Data_df <- read.csv(paste0(In_dir, "KS_IRIS_data_V2.0_colonygrowth.csv"),sep = ";", row.names = 1)


##==== Check if plate effect exists ===== ###
#outer wells grow better due to more access to nutrients
plot_plate_effect <- function(input_table, plate_ID = c("Lib_plate","Plate_medium","Tech_Rep","Plate_date","Passage"),target_col="size", SD_cutoff=NULL){
  uni_df <- unique(input_table[,plate_ID])
  sapply(1:nrow(uni_df),FUN=function(X){
    current_plate <- input_table[ident_rows(input_table[,plate_ID],search_vec = uni_df[X,]),]
    outer_wells <- current_plate[,"Well"][str_detect(current_plate[,"Well"],"(A)|(O)|([:upper:]1\\b)|(\\w{1}24)")] #get out wells of 384 well plate
    plot_df <- data.frame(Well=current_plate[,"Well"],Well_pos= ifelse(current_plate[,"Well"] %in% outer_wells,"Outer","Inner"),current_plate[,target_col])
    names(plot_df)[ncol(plot_df)] <- target_col
    ## exclude data points that are out of range one SD of
    if(!is.null(SD_cutoff)){ plot_df<- plot_df[(plot_df[,target_col] > median(plot_df[,target_col]) -SD_cutoff* sd(plot_df[,target_col]) ) & 
                                      (plot_df[,target_col] < median(plot_df[,target_col]) +SD_cutoff* sd(plot_df[,target_col])),]
      }
    boxplot( plot_df[,target_col]~plot_df[,"Well_pos"],ylab = target_col, col=c("sienna","skyblue1"),
             main=apply(uni_df[X,],1,paste,sep="",collapse="_"))
    

  })#end of sapply
  
}


pdf(paste0(Out_dir,"Analyse_outer.frame_V",current_version,".pdf"),width = 15,height = 15)
par(mfrow=c(4,4),mar=c(2,2,5,1),oma=c(1,1,1,1))
plot_plate_effect(input_table = Data_df,target_col = "size",SD_cutoff=1,
                  plate_ID = c("Lib_plate","Plate_medium","Tech_Rep","Plate_date","Passage","strain"))
dev.off()

##### Correct plate ######

correct_outer_frame <- function(input_table, plate_ID = c("Lib_plate","Plate_medium","Tech_Rep","Plate_date","Passage","strain"), 
                                correct_cols =c("opacity","size") ){
  uni_df <- unique(input_table[,plate_ID])
  out_df <- input_table
  outer_wells <- input_table[,"Well"][str_detect(input_table[,"Well"],"(A)|(O)|([:upper:]1\\b)|(\\w{1}24)")] #get out wells
  Well_pos <- input_table[,"Well"] %in% outer_wells
  #loop trough plates
  for (p in 1:nrow(uni_df)){
    target_rows <- ident_rows2(input_table[,plate_ID],search_vec = uni_df[p,])
    
    for ( i in 1: length(correct_cols)){
      
      #first correct outer wells
      current_rows <- (1:nrow(out_df) %in% target_rows) & Well_pos == T 
      out_df[current_rows,correct_cols[i]] <- out_df[current_rows,correct_cols[i]]/ median(out_df[current_rows,correct_cols[i]],na.rm = T )
      # correct inner  wells
      current_rows <- (1:nrow(out_df) %in% target_rows) & Well_pos == F 
      out_df[current_rows,correct_cols[i]] <- out_df[current_rows,correct_cols[i]]/ median(out_df[current_rows,correct_cols[i]],na.rm = T )
 
    }#loop trough columns to correct
  }
  return(out_df)
  
}

rim.correct_df <- correct_outer_frame(input_table = Data_df,correct_cols =c("opacity","size"))


### check if correction worked

pdf(paste0(Out_dir,"check_Analyse_outer.frame_V",current_version,".pdf"),width = 15,height = 15)
par(mfrow=c(4,4),mar=c(2,2,5,1),oma=c(1,1,1,1))
plot_plate_effect(input_table = rim.correct_df,target_col = "size",SD_cutoff=1,
                  plate_ID =  c("Lib_plate","Plate_medium","strain","Tech_Rep","Plate_date","Passage") )
dev.off()


#### Remove KS_2A SG1531 ######

#it appears that mixed Plate KS_2(A) causes problems with the correction. Better omit this condition from further analysis
rim.correct_df2<- rim.correct_df[!(rim.correct_df$Lib_plate == "KS_2A" &  rim.correct_df$strain == "SG1531"),]

###### Export #####
write.csv(rim.correct_df2,paste0(Out_dir,"KS_IRIS_data_frame_correct_V",current_version,".csv"))











