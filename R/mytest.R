mytest = function(){
  d = readData("/Users/jcur002/Dropbox/Code/git/rbayesfst/data/data_BB04.inp")
  b = new(BayesFst)
  b$setData(d)
  b$interaction = FALSE
  b$printData()
  b$printCounts()
  b$setRunParameters(101, 0.2, 0.05, 0.02, TRUE)
  b$run(25437, 11746, 21291)
  # newpi = c(0.34504698665347278, 0.63224240033200596, 0.02271061301452117)
  # fpar = c(353.01488819034944, 617.01675796349809, 29.968353846152503)
  # print(b$ldiriTest(fpar, newpi))
}