#' Restart a simulation on extinction.
#' 
#' @param expr
#' @param n number of times to repeat until success


restartOnExtinct <- function(expr, n = 10) {
  suppressWarnings({
  success <- T
  for (i in 1:n) {
    res <- tryCatch(expr,
                    error = function(e) {
                      print(sprintf("Population extinct. Restarting...", conditionMessage(e)))
                      if(i == n) print(sprintf("Population extinct ten times. Consider changing population size or selection parameters. Aborting...", conditionMessage(e)))
                      
                      success <<- F
                      e
                    }
    )
    if (success) break
  }
  res
  })
}




