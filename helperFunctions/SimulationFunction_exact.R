SimulationFunction_exact = function(number_taxa, number_quads, repeats, data1, data2, verbose=F){
  # number_taxa = 6
  # number_quads = 5
  # repeats = 100
  # data1 = copy(test1)
  # data2 = myTab_n06
  # verbose = T
  
  # Step 1: define variables
  n = number_taxa
  k = number_quads
  m = dim(data1$data)[1]
  
  # Step 2: some checks
  allQuads = data.table::as.data.table(t(combn(n,4)))
  stopifnot(dim(data1$data)[1] == dim(allQuads)[1])
  if(verbose == T) message("       Input data matches to given taxa size")
  
  # Step 3: creae all combinations 
  allCombis = data.table::as.data.table(t(combn(m,k)))
  rep2 = dim(allCombis)[1]
  if(rep2<repeats){
    message("       Testing all combination only once ...")
    repeats = rep2
  }else{
    message("       Working on ",repeats," repeats of ", dim(allCombis)[1]," combinations")
  }
  
  # Step 4: set up simulation
  set.seed(2015)
  x = sample(x = 1:rep2,size = repeats,replace = F)
  y = seq(1,repeats,by=ceiling(repeats/100))
  
  # Step 5: test all simulated sets
  dumTab = foreach(i = 1:repeats)%do%{
    #i=1
    if(is.element(i,y) & verbose==T){
      message("              Working on combination ",i)
    } 
    myX = x[i]
    myCombi = as.numeric(allCombis[myX,])
    data3 = copy(data1)
    data3$data[myCombi,status := "input"]
    data3$data[is.na(status),status := "unresolved"]
    
    # Step 5.1: exact test with 4-way partition property
    time0 = Sys.time()
    
    dummy1 = data3$data[status=="input",]
    myTab1 = copy(data2)
    myTab1[,count :=0]
    
    for(i in 1:dim(dummy1)[1]){
      #i=1
      myInput = dummy1[i,ctuple]
      filt = grepl(myInput,myTab1$allQuads)
      myTab1[filt==T,status := "covered"]
      myTab1[filt==T,count := count + 1]
      myTab1
    }
    
    filt3 = myTab1$count == 0
    if(sum(filt3)==0){
      test_1_res = "PhyloDec"
    } else {
      test_1_res = "NOT PhyloDec"
    }  
    time1 = Sys.time()
    diff_1 =  round(difftime(time1,time0,units = "sec"),3)

    # Step 5.2: FFTs
    time0 = Sys.time()
    test_FTT_sim = FTT_algorithmRed(data = data3$data, n=n, c=4,verbose = F)
    filt1 = test_FTT_sim$status == "unresolved"
    if(sum(filt1)==0){
      test_2_res = "FTT"
      test_2_rounds = max(test_FTT_sim$round)
    } else {
      test_2_res = "not FTT"
      test_2_rounds = max(test_FTT_sim$round) + 1
    } 
    time1 = Sys.time()
    diff_2 =  round(difftime(time1,time0,units = "sec"),3)
    
    # Step 5.3: summary
    quads = data3$data[status == "input", paste(ctuple, collapse ="|")]
    res = data.table(FWPP = test_1_res,
                     FWPP_time = diff_1,
                     FTT = test_2_res,
                     FTT_time = diff_2,
                     FTT_round = test_2_rounds,
                     input = quads)
    res
  }
  dumTab = rbindlist(dumTab)
  return(dumTab)
  
}
