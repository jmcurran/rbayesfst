mytest = function(){
  d = readData("data/data_BB04.inp")
  b = new(BayesFst)
  b$setData(d)
  b$interaction = TRUE
  b$printData()
  b$setRunParameters(2001, 0.2, 0.05, 0.02, TRUE)
  b$run(1)
}