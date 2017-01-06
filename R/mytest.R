mytest = function(){
  d = readData(system.file("extdata", "data_BB04.inp", packge = "rbayesfst"))
  b = new(BayesFst)
  b$setData(d)
  b$interaction = FALSE
  b$printData()
  b$printCounts()
  b$setRunParameters(2001, 0.2, 0.05, 0.02, TRUE)
  b$run(25437) #, 11746, 21291)
}