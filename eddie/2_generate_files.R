


RunID <- "a"
runOnEddie <- T
 

#~~ Parameters to be tested


test.frame <- expand.grid(list(
  pop.size            = c(10, 20, 50, 100, 150, 200, 500, 1000),
  loci.number         = c(10, 20, 50, 100, 150, 200, 500, 1000),
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
  
  writeLines(c(paste0("n.females           <- ", test.frame$pop.size[i]), 
               paste0("n.males             <- ", test.frame$pop.size[i]),   
               paste0("n.offspring         <- 2"),   
               paste0("male.sel.thresh     <- ", test.frame$sel.thresh[i]),
               paste0("female.sel.thresh   <- 1"),
               paste0("generations         <- ", test.frame$generations[i]), 
               paste0("modifier.found.freq <- ", test.frame$modifier.found.freq[i]),
               paste0("iterations          <- ", test.frame$iterations[i]),
               paste0("n.loci              <- ", test.frame$loci.number[i]),
               paste0("n.found.hap         <- 100"),
               paste0("rate.increase       <- 1.5"),
               paste0("modifier            <- \"", test.frame$modifier[i], "\""),
               paste0("model.name          <- \"", test.frame$Model.Name[i], "\"")),
             paste0(test.frame$Model.Name[i],".R"))
  
  system(paste("cat 1_model_template.R >>", paste0(test.frame$Model.Name[i],".R"))) 
    
  if(runOnEddie == T){
    
    writeLines(paste0("#!/bin/sh\n", 
                      "\n",
                      "#$ -cwd\n",
                      ifelse(test.frame$pop.size[i] >= 500, "#$ -l h_rt=10:00:00\n", "#$ -l h_rt=04:00:00\n"),
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


