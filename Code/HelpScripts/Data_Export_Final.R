
upload_to_drive <- function(repo = "Data_Novotny/AmpliconSeqAnalysis") {
  
  require(googledrive)
  
  #########################
  # Upload COI to the drive
  
  drive_rm(file.path(repo, "Taxtab_COI.rds"))
  drive_upload("Data_export/Taxtab_COI.rds", file.path(repo, "Taxtab_COI.rds"))
  
  drive_rm(file.path(repo, "ASV_COI.rds"))
  drive_upload("Data_export/ASV_COI.rds", file.path(repo, "ASV_COI.rds"))
  
  drive_rm(file.path(repo, "Samptab_COI.rds"))
  drive_upload("Data_export/Samptab_COI.rds", file.path(repo, "Samptab_COI.rds"))
  
  
  #########################
  # Upload 18S to the drive
  
  drive_rm(file.path(repo, "Taxtab_18S.rds"))
  drive_upload("Data_export/Taxtab_18S.rds", file.path(repo, "Taxtab_18S.rds"))
  
  drive_rm(file.path(repo, "ASV_18S.rds"))
  drive_upload("Data_export/ASV_18S.rds", file.path(repo, "ASV_18S.rds"))
  
  drive_rm(file.path(repo, "Samptab_18S.rds"))
  drive_upload("Data_export/Samptab_18S.rds", file.path(repo, "Samptab_18S.rds"))
  
  
  #########################
  # Upload 12 to the drive
  
  drive_rm(file.path(repo, "Taxtab_12S.rds"))
  drive_upload("Data_export/Taxtab_12S.rds", file.path(repo, "Taxtab_12S.rds"))
  
  drive_rm(file.path(repo, "ASV_12S.rds"))
  drive_upload("Data_export/ASV_12S.rds", file.path(repo, "ASV_12S.rds"))
  
  drive_rm(file.path(repo, "Samptab_12S.rds"))
  drive_upload("Data_export/Samptab_12S.rds", file.path(repo, "Samptab_12S.rds"))
} 
  repo = "Data_Novotny/AmpliconSeqAnalysis"
  drive_rm(file.path(repo, "Samptab_12S.rds"))
  drive_upload("Data_export/Samptab_12S.rds", file.path(repo, "Samptab_12S.rds"))
