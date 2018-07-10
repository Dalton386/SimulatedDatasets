# interactive()
# 
# pttnum = as.numeric(readline("patient number is: "))
# pwynum = as.numeric(readline("pathway number is: "))
# sbtnum = as.numeric(readline("subtype number is: "))

library(MCMCpack)

dir.create("./simulatedData/")
dir.create("./simulatedData/Uniform/")
dir.create("./simulatedData/Dirichlet/")

#patient number and pathway number
pttnum = c(100,500)
pwynum = c(3,6,10)

for (ptt in pttnum) {
  for (pwy in pwynum) {
    # subtype number is 3 times more than pathway number
    sbt = 3 * pwy
    drchlt = rdirichlet(1, rep(0.7, sbt))[1,]
    # patient-subtype distribution with respect to Dirichlet distribution
    ddstri = matrix(data = sample(x = sbt, size = ptt, replace = TRUE, prob = drchlt), ncol = 1, nrow = ptt)
    # patient-subtype distribution with respect to uniform distribution
    udstri = matrix(data = sample(x = sbt, size = ptt, replace = TRUE, prob = rep(1/sbt, sbt)), ncol = 1, nrow = ptt)
    
    isinteresting = FALSE
    while(!isinteresting){
      # simulate subtype-pathway distribution
      sbt_pwy = matrix(0L, ncol = pwy, nrow = sbt)
      dec = 1 / sbt
      p = 0.5
      # generate columns in an order of decreasing expections
      for (i in 1:pwy){
        col = sample(x = c(0,1), size = sbt, replace = TRUE, prob = c(1-p, p))
        p = p - dec
        sbt_pwy[,i] = col
      }
      
      # terminate if some subtypes have the same pathway combination;
      # otherwise, simulate it again
      for (i in 1:(sbt-1)){
        for (j in (i+1):sbt){
          if (identical(sbt_pwy[i,], sbt_pwy[j,])){
            isinteresting = TRUE
            break
          }
        }
        if (isinteresting)
          break
      }
    }

    # record outputs in separate files
    write.matrix(sbt_pwy, file = paste0("./simulatedData/Uniform/usim",ptt,"_sub-",sbt,"_path-",pwy,".profile"))
    write.matrix(udstri, file = paste0("./simulatedData/Uniform/usim",ptt,"_sub-",sbt,".subs"))

    write.matrix(sbt_pwy, file = paste0("./simulatedData/Dirichlet/dsim",ptt,"_sub-",sbt,"_path-",pwy,".profile"))
    write.matrix(ddstri, file = paste0("./simulatedData/Dirichlet/dsim",ptt,"_sub-",sbt,".subs"))
  }
}
