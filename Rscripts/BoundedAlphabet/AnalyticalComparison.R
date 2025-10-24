M = 1e5
s = 1.02
pj = sim_TruncatedZipfs_features(M = M, s = s)

(1 - pj)^nn
# Vario n

nn = seq(1,100,by = 1)
Bonf = 1/nn * log(M/nn)

S = sum(pj)
Or14 = rep(0,length(nn))
Or26 = rep(0,length(nn))
for( n in nn ){
  Or14[n] = 1/(n-log(n)) * log( (1/0.05) * sum((1 - pj)^n) )
  Or26[n] = compute_UB_intersection(n, 0.05, 1e-36, S)
}


plot(x = nn, y = Bonf, type = "l", lwd = 3, col = "darkred",
     ylab = "bound", xlab = "n")
points(x = nn, y = Or14, type = "l", lwd = 3, col = "darkgreen")
points(x = nn, y = Or26, type = "l", lwd = 3, col = "darkblue")
legend("topright", c("Bonf.","Eq.(14)", "Eq.(26)"),
       col = c("darkred","darkgreen","darkblue"), 
      lwd = 3)
