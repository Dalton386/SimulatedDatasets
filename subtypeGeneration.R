# interactive()
# 
# pttnum = as.numeric(readline("patient number is: "))
# pwynum = as.numeric(readline("pathway number is: "))
# sbtnum = as.numeric(readline("subtype number is: "))

library(MCMCpack)

dir.create("./simulatedData/")
dir.create("./simulatedData/Uniform/")
dir.create("./simulatedData/Dirichlet/")

pttnum = c(100,500)
pwynum = c(3,6,10)

for (ptt in pttnum) {
  for (pwy in pwynum) {
    sbt = 3 * pwy
    drchlt = rdirichlet(1, rep(0.7, sbt))[1,]
    ddstri = matrix(data = sample(x = sbt, size = ptt, replace = TRUE, prob = drchlt), ncol = 1, nrow = ptt)
    udstri = matrix(data = sample(x = sbt, size = ptt, replace = TRUE, prob = rep(1/sbt, sbt)), ncol = 1, nrow = ptt)
    
    # no additional noise
    isinteresting = FALSE
    while(isinteresting){
      sbt_pwy = matrix(0L, ncol = pwy, nrow = sbt)
      dec = 1 / sbt
      p = 0.5
      for (i in 1:pwy){
        col = sample(x = c(0,1), size = sbt, replace = TRUE, prob = c(1-p, p))
        p = p - dec
        sbt_pwy[,i] = col
      }
      
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

    write.matrix(sbt_pwy, file = paste0("./simulatedData/Uniform/usim",ptt,"_sub-",sbt,"_path-",pwy,".profile"))
    write.matrix(udstri, file = paste0("./simulatedData/Uniform/usim",ptt,"_sub-",sbt,".subs"))

    write.matrix(sbt_pwy, file = paste0("./simulatedData/Dirichlet/dsim",ptt,"_sub-",sbt,"_path-",pwy,".profile"))
    write.matrix(ddstri, file = paste0("./simulatedData/Dirichlet/dsim",ptt,"_sub-",sbt,".subs"))
  }
}
