conditions <- data.frame()
for (gam in c(3,15,50,150,500,1000,2500,5000)){
  for (snp in c(5000, 30000, 100000)){
    for (cov in c(0.001, 0.01, 0.1, 0.223, 0.357, 0.511, 0.693, 1.204, 2.303)){
      for (se in c(0.001, 0.005, 0.05)){
        for (avgr in c(0.6, 1, 3)){
          for (rsd in c(42, 357, 1848)){
            new_row <- data.frame(gam = gam, snp=as.character(as.integer(snp)), cov=cov, se=se, avgr=avgr, rsd=rsd)
            conditions <- rbind(conditions, new_row)
          }
        }
      }
    }
  }
}

write.table(conditions, file="~/mccoyLab/rhapsodi_work/hapcut2sims/conditions.txt", quote=FALSE,row.names = FALSE, col.names=FALSE, sep="\t")