#'@export
calculateRSquare<- function(graph) {
  # calculate degree
  #library("igraph");
  d = degree(graph, mode = "all") # node'ların dugumlerinin derecelerini bulur
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE) # herbir derecenin gelme olasiligini kac kere o dereceden dugum var ona gore hesaplar. mesela 20 dugumun 3 tanesinin derecesi 1 ise, dugum derecesinin 1 olma olasiligi % 0.15 gibi
  degree = 1:max(d) # 1'den maksimumum dugun sayisina kadar git
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  #   print(paste("Alpha =", round(alpha, 3)))
  #   print(paste("R square =", round(R.square, 3)))
  rsquare=round(R.square, 3);
  
  #   plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", 
  #        col = 1, main = "Degree Distribution")
  #   curve(power.law.fit, col = "red", add = T, n = length(d))
  
  return (rsquare);
}

