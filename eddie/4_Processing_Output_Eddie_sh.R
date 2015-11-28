RunID <- "a"
runOnEddie = T
#~~ Parameters to be tested


test.frame <- expand.grid(list(
  pop.size            = c(20, 50, 100, 150, 200, 500),
  loci.number         = c(20, 50, 100, 150, 200, 500),
  sel.thresh          = seq(0, 1, 0.2),
  modifier.found.freq = seq(0, 1, 0.2),
  generations         = 200,
  n.found.hap         = 100,
  iterations          = 100,
  modifier            = c("RNF212", "PRDM9")))

test.frame <- unique(test.frame)
test.frame$Model.Name <- paste0(RunID, 
                                "_n", test.frame$pop.size,
                                "_l", test.frame$loci.number,
                                "_s", test.frame$sel.thresh,
                                "_m", test.frame$modifier.found.freq,
                                "_g", test.frame$generations,
                                "_i", test.frame$iterations,
                                "_" , test.frame$modifier)

for(i in 1:nrow(test.frame)){
  
  writeLines(paste0("x <- \"", test.frame$Model.Name[i], "\"\n",
                    "mod.freq <- ", test.frame$modifier.found.freq[i], "\n"),
             paste0(test.frame$Model.Name[i], ".R"))

  system(paste("cat 4_model_template.R >>", paste0(test.frame$Model.Name[i],".R"))) 
  
  if(runOnEddie == T){
    
    writeLines(paste0("#!/bin/sh\n", 
                      "\n",
                      "#$ -cwd\n",
                      "#$ -l h_rt=04:00:00\n",
                      "#$ -V\n",
                      "#$ -l h_vmem=5200M\n",
                      "\n",
                      ". /etc/profile.d/modules.sh\n",
                      "module load R/3.1.1\n",
                      "R CMD BATCH ", test.frame$Model.Name[i], ".R --no-restore --no-save\n"),
               paste0(test.frame$Model.Name[i], ".sh"))
    
    system(paste0("qsub ", test.frame$Model.Name[i], ".sh"))
    
  }
  
  
}

