panel.cor = function(x, y, digits = 2, prefix = "", cex.cor){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, method="pearson",use="pairwise.complete.obs"),2)
  txt <- format(c(r, 0.12), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  text(0.5, 0.5, txt, cex=1.5)
}

panel.points = function(x,y){
  points(x,y,pch=20)
}