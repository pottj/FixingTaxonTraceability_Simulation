Evaluation 4: get tight upper bound
================
Janne Pott
2024-01-09

# Introduction

The upper bound not decisive sets is given as

$$\binom{n}{4} - (n-4)+1 $$

E.g., for $n=6$, the bound is $\binom{6}{4} - 3 = 12$, so there are sets
of size 12 which are neither decisive nor fixing taxon traceable. Here I
just want to check if this bound is always met.

(In the simulation, we did not test all possible combinations for the
maximal $k_n$, hence we might have missed those sets in the direct test,
but tested them with a lower starting $k$ with additional resolved
quadruples.)

``` r
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
#> Warning: package 'ggplot2' was built under R version 4.2.3
#> Warning: package 'cowplot' was built under R version 4.2.3
```

# Load data

``` r
load("../results/02_SimulationResults_n06.RData")
load("../results/02_SimulationResults_n07.RData")
load("../results/02_SimulationResults_n08.RData")
load("../results/02_SimulationResults_n09.RData")
load("../results/02_SimulationResults_n10.RData")

SimulationResults_n06[,n := 6]
SimulationResults_n07[,n := 7]
SimulationResults_n08[,n := 8]
SimulationResults_n09[,n := 9]
SimulationResults_n10[,n := 10]

sim = rbind(SimulationResults_n06,SimulationResults_n07,SimulationResults_n08,SimulationResults_n09,SimulationResults_n10)
Sim = copy(sim)

dumTab = foreach(i=1:dim(Sim)[1])%do%{
  #i=357
  myRow = Sim[i,]
  myk = myRow$k
  myn = myRow$n
  
  if(myn<10){
      filefn = paste0("../results/02_SimulationsData_n0",myn,"/SimulationResults_n",myn,"_k",myk,"_rep10000.RData")
  }else{
      filefn = paste0("../results/02_SimulationsData_n",myn,"/SimulationResults_n",myn,"_k",myk,"_rep10000.RData")
  }
  load(filefn)
  myTest = myTest[FWPP == "NOT PhyloDec",]
  if(dim(myTest)[1]==0){
    res = data.table(n=myn, 
                     k=myk, 
                     FWPP = T, 
                     FixingTaxa = T,
                     input = NA)
  }else{
    res = data.table(n=myn, 
                     k=myk, 
                     FWPP = F, 
                     FixingTaxa = F,
               input = myTest$input)

      }
  res
}

PhyloDecSets_notFTT = rbindlist(dumTab)
PhyloDecSets_notFTT[,table(n,FWPP)]
#>     FWPP
#> n     FALSE   TRUE
#>   6   27403      0
#>   7  139172      0
#>   8  284743      3
#>   9  505899      9
#>   10 822868     23
PhyloDecSets_notFTT = PhyloDecSets_notFTT[FWPP==F,]
```

# Test data

I want to test if I have sets starting with some input $k$ and ending
with additional resolved quadruples so that the sum of input and solved
quadruples equals the tight bound.

``` r
dumTab = foreach(i = 6:10)%do%{
  #i=6
  dumTab1 = copy(PhyloDecSets_notFTT)
  dumTab1 = dumTab1[n==i,]
  
  setorder(dumTab1,-k)
  
  tracer = FALSE
  theo_threshold = choose(i,4) - (i-3)
    
  while(tracer == F & dim(dumTab1)[1]>0){
    dumTab2 = dumTab1[1,]
    myInput = dumTab1[1,input]
    myInput = unlist(strsplit(myInput,split = "[|]"))
    myOutput = data.table(myInput = myInput)
    myFn = paste0("../results/03_4_GetUpperBound/Example.txt")
    fwrite(myOutput,file=myFn,col.names = F)
    TestData = FTT_createInput(fn = myFn, sepSym = "_",c = 4)
    TestFTT = FTT_algorithmRed(data=TestData$data,
                               verbose = F,c = 4,n = myRow$n)
        
    testStatus = TestFTT[,.N,by=status]
    sim_res = testStatus[status != "unresolved",sum(N)]
    
    if(theo_threshold == sim_res){
        myResult = dumTab1[1,]
        myResult[,green := sim_res]
        myResult[,rounds := max(TestFTT$round)]
        print(myResult)
        myResult
        tracer = TRUE
    }else{
      dumTab1 = dumTab1[-1,]
    }
  }
  myResult
}
#> Input contains 12 sets with 6 different taxa. 
#> The largest set has 4 taxa.
#>    n  k  FWPP FixingTaxa
#> 1: 6 12 FALSE      FALSE
#>                                                                                              input
#> 1: 1_2_3_4|1_2_3_5|1_2_3_6|1_2_4_5|1_2_4_6|1_2_5_6|1_3_4_5|1_3_4_6|1_3_5_6|2_3_4_5|2_3_4_6|2_3_5_6
#>    green rounds
#> 1:    12      0
#> Input contains 31 sets with 7 different taxa. 
#> The largest set has 4 taxa.
#>    n  k  FWPP FixingTaxa
#> 1: 7 31 FALSE      FALSE
#>                                                                                                                                                                                                                                                      input
#> 1: 1_2_4_5|1_2_4_6|1_2_4_7|1_2_5_6|1_2_5_7|1_2_6_7|1_3_4_5|1_3_4_6|1_3_4_7|1_3_5_6|1_3_5_7|1_3_6_7|1_4_5_6|1_4_5_7|1_4_6_7|1_5_6_7|2_3_4_5|2_3_4_6|2_3_4_7|2_3_5_6|2_3_5_7|2_3_6_7|2_4_5_6|2_4_5_7|2_4_6_7|2_5_6_7|3_4_5_6|3_4_5_7|3_4_6_7|3_5_6_7|4_5_6_7
#>    green rounds
#> 1:    31      0
#> Input contains 62 sets with 8 different taxa. 
#> The largest set has 4 taxa.
#>    n  k  FWPP FixingTaxa
#> 1: 8 62 FALSE      FALSE
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              input
#> 1: 1_2_3_5|1_2_3_6|1_2_3_7|1_2_3_8|1_2_4_5|1_2_4_6|1_2_4_7|1_2_4_8|1_2_5_6|1_2_5_7|1_2_5_8|1_2_6_7|1_2_6_8|1_2_7_8|1_3_5_6|1_3_5_7|1_3_5_8|1_3_6_7|1_3_6_8|1_3_7_8|1_4_5_6|1_4_5_7|1_4_5_8|1_4_6_7|1_4_6_8|1_4_7_8|1_5_6_7|1_5_6_8|1_5_7_8|1_6_7_8|2_3_4_5|2_3_4_6|2_3_4_7|2_3_4_8|2_3_5_6|2_3_5_7|2_3_6_7|2_3_6_8|2_3_7_8|2_4_5_6|2_4_5_7|2_4_5_8|2_4_6_7|2_4_6_8|2_4_7_8|2_5_6_7|2_5_6_8|2_5_7_8|2_6_7_8|3_4_5_6|3_4_5_8|3_4_6_7|3_4_6_8|3_4_7_8|3_5_6_7|3_5_6_8|3_5_7_8|3_6_7_8|4_5_6_7|4_5_6_8|4_5_7_8|4_6_7_8
#>    green rounds
#> 1:    65      1
#> Input contains 113 sets with 9 different taxa. 
#> The largest set has 4 taxa.
#>    n   k  FWPP FixingTaxa
#> 1: 9 113 FALSE      FALSE
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      input
#> 1: 1_2_3_4|1_2_3_5|1_2_3_7|1_2_3_8|1_2_3_9|1_2_4_5|1_2_4_6|1_2_4_7|1_2_4_8|1_2_4_9|1_2_5_6|1_2_5_7|1_2_5_9|1_2_6_7|1_2_6_8|1_2_6_9|1_2_7_8|1_2_7_9|1_2_8_9|1_3_4_5|1_3_4_6|1_3_4_7|1_3_4_8|1_3_4_9|1_3_5_6|1_3_5_7|1_3_5_8|1_3_5_9|1_3_6_7|1_3_6_8|1_3_6_9|1_3_7_8|1_3_7_9|1_3_8_9|1_4_5_6|1_4_5_7|1_4_5_8|1_4_6_7|1_4_6_8|1_4_6_9|1_4_7_8|1_4_7_9|1_4_8_9|1_5_6_7|1_5_6_8|1_5_6_9|1_5_7_8|1_5_7_9|1_5_8_9|1_6_7_8|1_6_7_9|1_7_8_9|2_3_4_5|2_3_4_6|2_3_4_7|2_3_4_9|2_3_5_6|2_3_5_7|2_3_5_9|2_3_6_7|2_3_6_8|2_3_6_9|2_3_7_8|2_3_7_9|2_3_8_9|2_4_5_6|2_4_5_7|2_4_5_9|2_4_6_7|2_4_6_8|2_4_6_9|2_4_7_8|2_4_7_9|2_4_8_9|2_5_6_7|2_5_6_9|2_5_7_9|2_6_7_8|2_6_8_9|2_7_8_9|3_4_5_6|3_4_5_7|3_4_5_8|3_4_5_9|3_4_6_7|3_4_6_8|3_4_6_9|3_4_7_8|3_4_7_9|3_4_8_9|3_5_6_7|3_5_6_8|3_5_6_9|3_5_7_8|3_5_7_9|3_5_8_9|3_6_7_8|3_6_7_9|3_6_8_9|3_7_8_9|4_5_6_7|4_5_6_8|4_5_6_9|4_5_7_8|4_5_7_9|4_6_7_8|4_6_8_9|4_7_8_9|5_6_7_8|5_6_7_9|5_6_8_9|5_7_8_9|6_7_8_9
#>    green rounds
#> 1:   120      1
#> Input contains 183 sets with 10 different taxa. 
#> The largest set has 4 taxa.
#>     n   k  FWPP FixingTaxa
#> 1: 10 183 FALSE      FALSE
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 input
#> 1: 1_2_3_4|1_2_3_5|1_2_3_7|1_2_3_8|1_2_3_9|1_2_3_10|1_2_4_5|1_2_4_8|1_2_4_9|1_2_4_10|1_2_5_7|1_2_5_9|1_2_5_10|1_2_7_8|1_2_7_9|1_2_7_10|1_2_8_9|1_2_8_10|1_2_9_10|1_3_4_5|1_3_4_7|1_3_4_8|1_3_4_9|1_3_4_10|1_3_5_6|1_3_5_7|1_3_5_8|1_3_5_9|1_3_5_10|1_3_6_7|1_3_6_8|1_3_6_9|1_3_6_10|1_3_7_8|1_3_7_9|1_3_7_10|1_3_8_9|1_3_8_10|1_3_9_10|1_4_5_6|1_4_5_8|1_4_5_9|1_4_5_10|1_4_6_7|1_4_6_8|1_4_6_9|1_4_6_10|1_4_7_8|1_4_7_9|1_4_7_10|1_4_8_9|1_4_8_10|1_4_9_10|1_5_6_7|1_5_6_8|1_5_6_9|1_5_6_10|1_5_7_9|1_5_7_10|1_5_8_9|1_5_9_10|1_6_7_8|1_6_7_9|1_6_7_10|1_6_8_9|1_6_8_10|1_6_9_10|1_7_8_9|1_7_8_10|1_7_9_10|1_8_9_10|2_3_4_5|2_3_4_6|2_3_4_7|2_3_4_8|2_3_4_9|2_3_4_10|2_3_5_7|2_3_5_8|2_3_5_9|2_3_5_10|2_3_6_7|2_3_6_8|2_3_6_9|2_3_6_10|2_3_7_8|2_3_7_9|2_3_7_10|2_3_8_9|2_3_8_10|2_3_9_10|2_4_5_6|2_4_5_7|2_4_5_8|2_4_5_9|2_4_6_7|2_4_6_8|2_4_6_10|2_4_7_8|2_4_7_9|2_4_7_10|2_4_8_9|2_4_8_10|2_5_6_7|2_5_6_8|2_5_6_9|2_5_6_10|2_5_7_8|2_5_7_9|2_5_7_10|2_5_8_9|2_5_8_10|2_5_9_10|2_6_7_8|2_6_7_10|2_6_8_9|2_6_9_10|2_7_8_9|2_7_8_10|2_8_9_10|3_4_5_6|3_4_5_7|3_4_5_8|3_4_5_9|3_4_5_10|3_4_6_7|3_4_6_8|3_4_6_9|3_4_6_10|3_4_7_8|3_4_7_9|3_4_7_10|3_4_8_9|3_4_8_10|3_5_6_7|3_5_6_8|3_5_6_9|3_5_6_10|3_5_7_10|3_5_8_9|3_5_8_10|3_5_9_10|3_6_7_8|3_6_7_9|3_6_7_10|3_6_8_9|3_6_8_10|3_6_9_10|3_7_8_9|3_7_8_10|3_7_9_10|3_8_9_10|4_5_6_7|4_5_6_8|4_5_6_9|4_5_6_10|4_5_7_8|4_5_7_9|4_5_7_10|4_5_8_9|4_5_8_10|4_5_9_10|4_6_7_8|4_6_7_9|4_6_7_10|4_6_8_9|4_6_9_10|4_7_8_9|4_7_9_10|4_8_9_10|5_6_7_8|5_6_7_9|5_6_7_10|5_6_8_10|5_6_9_10|5_7_8_10|5_7_9_10|5_8_9_10|6_7_8_9|6_7_8_10|6_7_9_10|6_8_9_10|7_8_9_10
#>    green rounds
#> 1:   203      1

myResult = rbindlist(dumTab)
myResult[,c(1,2,6,7)]
#>     n   k green rounds
#> 1:  6  12    12      0
#> 2:  7  31    31      0
#> 3:  8  62    65      1
#> 4:  9 113   120      1
#> 5: 10 183   203      1
save(myResult,file = "../results/03_4_UpperBound.RData")
```

Okay, it worked :)

# Session Info

``` r
sessionInfo()
#> R version 4.2.2 (2022-10-31 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 19045)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=English_United Kingdom.utf8 
#> [2] LC_CTYPE=English_United Kingdom.utf8   
#> [3] LC_MONETARY=English_United Kingdom.utf8
#> [4] LC_NUMERIC=C                           
#> [5] LC_TIME=English_United Kingdom.utf8    
#> 
#> attached base packages:
#> [1] grid      stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] cowplot_1.1.1           gtable_0.3.3            ggplot2_3.4.1          
#> [4] FixingTaxonTraceR_0.0.1 foreach_1.5.2           data.table_1.14.8      
#> 
#> loaded via a namespace (and not attached):
#>  [1] rstudioapi_0.14  knitr_1.42       magrittr_2.0.3   tidyselect_1.2.0
#>  [5] munsell_0.5.0    colorspace_2.1-0 R6_2.5.1         rlang_1.1.0     
#>  [9] fastmap_1.1.1    fansi_1.0.4      dplyr_1.1.0      tools_4.2.2     
#> [13] xfun_0.37        utf8_1.2.3       cli_3.6.0        withr_2.5.0     
#> [17] htmltools_0.5.4  iterators_1.0.14 yaml_2.3.7       digest_0.6.31   
#> [21] tibble_3.2.0     lifecycle_1.0.3  vctrs_0.5.2      codetools_0.2-18
#> [25] glue_1.6.2       evaluate_0.20    rmarkdown_2.20   compiler_4.2.2  
#> [29] pillar_1.9.0     generics_0.1.3   scales_1.2.1     pkgconfig_2.0.3
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
#> 
#> TOTAL TIME : 0.64 minutes
```
