# source files
source("scripts/general_functions.R")

config = list("methods" = c("fastcar"), "process" = c("decontaminate"))
load_libraries()

dir = "/data/Perkins/Mouse_Kidney/data/CellRangerCompressed"


for (s_id in c("BG5_BG22C", "BG56_BG22C_MeOH")) {
  #cellMatrix = read.cell.matrix(paste(dir, s_id, "filtered_gene_bc_matrices/", sep="/"))
  fullMatrix = read.full.matrix(paste(dir, s_id, "raw_gene_bc_matrices","mm10", sep="/"))
  
  for (i in seq(5)) {
    ambProfile = describe.ambient.RNA.sequence(fullCellMatrix = fullMatrix, 
                                               start = 10, 
                                               stop = 500, 
                                               by = 10, 
                                               contaminationChanceCutoff = i * 0.02) # 0.02, 0.04, 0.06, 0.08, 0.10
     
    p = plot.ambient.profile(ambProfile)
  }
}
