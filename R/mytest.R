mytest = function(){
  bd = readData(system.file("extdata", "data_BB04.inp", package = "rbayesfst"))
  bf = init(bd)
  setRunParams(bf, numOut = 21)
  r = sample(bf, seed = 123)
  plot(r$beta)
}