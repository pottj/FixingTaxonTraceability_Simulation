
<!-- 02_Simulation_n07.md is generated from 02_Simulation_n07.Rmd. Please edit that file -->

# Introduction

<!-- badges: start -->
<!-- badges: end -->

I want to run a simulation for n=7 taxa. There are 35 possible
quadruples with n=7 taxa, and I want to test 100000 possible
combinations for each k out of 35 quadruples for phylogenetic
decisiveness.

For each random set of quadruples, I test both the
4-way-partition-property and the Fixing Taxon Traceability algorithm
(**FTT**, test), implemented in the R package **FixingTaxonTraceR**.

# Initialize

I use a file named *SourceFile.R* that contains all relevant R packages
and user-/server-specific path to the R library. If using this code, you
must make all the necessary changes within the template source file.

``` r
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
source("../helperFunctions/SimulationFunction.R")
```

# Get input data

In *myTab_n7*, there are all possible Four-way-partitions (4WPP) given
n=7 taxa.

``` r
# Load initial input data
test1 = FTT_createInput(fn="../data/S7_Decisive.txt",sepSym = "_",c = 4)
#> Input contains 14 sets with 7 different taxa. 
#> The largest set has 4 taxa.

# Sanity check of FTT package
test_FTT = FTT_algorithmRed(data = test1$data,verbose = T, c=4, n=7)
#> Using 14 of 35 4-tuples as input for algorithm (7 unique taxa). 
#>  This leaves 21 4-tuples unsolved.
#> In round #1, 0 4-tuples could be resolved ...
#> NOT RESOLVABLE VIA THIS ALGORITHM! 
#>  There are 21 remaining cross 4-tuples.
test_NRC = FTT_findNRC(data = test1)
#> THERE IS ONLY RAINBOW COLORING POSSIBLE OR 3-WAY COLORING
#>  It follows that the set is phylogenetically decisive

load("../results/01_partitions/partitions_n07.RData")
```

# Test-Loop with 10 combinations per k

To test less $k$, I use as lower bound the minimal triple covering,
$\frac{1}{4}\binom{n}{3}$. As upper bound I use
$\binom{n}{4} - (n-4) -1$, as above this $k$ all sets are phylogenetic
decisive.

``` r
test1$data[,status:=NA]
n=7
LowerBound = choose(n,3)/4
UpperBound = choose(n,4) - n + 3

dumTab = foreach(j=c(LowerBound:UpperBound))%do%{
  # j=20
  message("\nWorking on n=7, k=",j)
  time1 = Sys.time()
  myTest = SimulationFunction(number_taxa = n, 
                                    number_quads = j,
                                    repeats = 10,
                                    data1 = test1,
                                    data2 = myTab_n07,
                                    verbose = T)
  time2 = Sys.time()
  x0 = as.numeric(round(difftime(time2,time1,units = "mins"),3))
  message("       Total time for n=7, k=",j," & rep = 10: " ,round(difftime(time2,time1,units = "mins"),3)," minutes")
  
  outfn = paste0("../temp/02_SimulationsData_n07/SimulationResults_n7_k",j,".RData")
  save(myTest,file = outfn)
  
  tab = table(myTest$FWPP,myTest$FTT)
  x2 = myTest[FWPP == "PhyloDec" & FTT == "FTT",.N]
  x3 = myTest[FWPP == "NOT PhyloDec",.N]
  x4 = myTest[FWPP == "PhyloDec",.N]
  
  message("       There were ",x2," of ",x4," sets identified by Fischers algorith as decisive (",round(x2/x4,4)*100,"%)")
  
  res = data.table(k = j,
                   time = x0,
                   NR_NotPhyloDec = x3,
                   NR_PhyloDec = x4,
                   NR_FTT = x2)
  res
}
#> 
#> Working on n=7, k=8.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (8.75) muss ganze Zahl sein, gerundet auf 9
#>        Working on 10 repeats of 70607460 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=8.75 & rep = 10: 0.077 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=7, k=9.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (9.75) muss ganze Zahl sein, gerundet auf 10
#>        Working on 10 repeats of 183579396 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=9.75 & rep = 10: 0.075 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=7, k=10.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (10.75) muss ganze Zahl sein, gerundet auf 11
#>        Working on 10 repeats of 417225900 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=10.75 & rep = 10: 0.076 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=7, k=11.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (11.75) muss ganze Zahl sein, gerundet auf 12
#>        Working on 10 repeats of 834451800 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=11.75 & rep = 10: 0.083 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=7, k=12.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (12.75) muss ganze Zahl sein, gerundet auf 13
#>        Working on 10 repeats of 1476337800 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=12.75 & rep = 10: 0.073 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=7, k=13.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (13.75) muss ganze Zahl sein, gerundet auf 14
#>        Working on 10 repeats of 2319959400 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=13.75 & rep = 10: 0.106 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=7, k=14.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (14.75) muss ganze Zahl sein, gerundet auf 15
#>        Working on 10 repeats of 3247943160 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=14.75 & rep = 10: 0.118 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=7, k=15.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (15.75) muss ganze Zahl sein, gerundet auf 16
#>        Working on 10 repeats of 4059928950 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=15.75 & rep = 10: 0.114 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=7, k=16.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (16.75) muss ganze Zahl sein, gerundet auf 17
#>        Working on 10 repeats of 4537567650 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=16.75 & rep = 10: 0.108 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=7, k=17.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (17.75) muss ganze Zahl sein, gerundet auf 18
#>        Working on 10 repeats of 4537567650 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=17.75 & rep = 10: 0.126 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=7, k=18.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (18.75) muss ganze Zahl sein, gerundet auf 19
#>        Working on 10 repeats of 4059928950 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=18.75 & rep = 10: 0.106 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=7, k=19.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (19.75) muss ganze Zahl sein, gerundet auf 20
#>        Working on 10 repeats of 3247943160 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=19.75 & rep = 10: 0.12 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=7, k=20.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (20.75) muss ganze Zahl sein, gerundet auf 21
#>        Working on 10 repeats of 2319959400 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=20.75 & rep = 10: 0.095 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=7, k=21.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (21.75) muss ganze Zahl sein, gerundet auf 22
#>        Working on 10 repeats of 1476337800 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=21.75 & rep = 10: 0.093 minutes
#>        There were 2 of 3 sets identified by Fischers algorith as decisive (66.67%)
#> 
#> Working on n=7, k=22.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (22.75) muss ganze Zahl sein, gerundet auf 23
#>        Working on 10 repeats of 834451800 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=22.75 & rep = 10: 0.087 minutes
#>        There were 4 of 5 sets identified by Fischers algorith as decisive (80%)
#> 
#> Working on n=7, k=23.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (23.75) muss ganze Zahl sein, gerundet auf 24
#>        Working on 10 repeats of 417225900 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=23.75 & rep = 10: 0.069 minutes
#>        There were 6 of 6 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=7, k=24.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (24.75) muss ganze Zahl sein, gerundet auf 25
#>        Working on 10 repeats of 183579396 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=24.75 & rep = 10: 0.059 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=7, k=25.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (25.75) muss ganze Zahl sein, gerundet auf 26
#>        Working on 10 repeats of 70607460 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=25.75 & rep = 10: 0.058 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=7, k=26.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (26.75) muss ganze Zahl sein, gerundet auf 27
#>        Working on 10 repeats of 23535820 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=26.75 & rep = 10: 0.045 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=7, k=27.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (27.75) muss ganze Zahl sein, gerundet auf 28
#>        Working on 10 repeats of 6724520 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=27.75 & rep = 10: 0.044 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=7, k=28.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (28.75) muss ganze Zahl sein, gerundet auf 29
#>        Working on 10 repeats of 1623160 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=28.75 & rep = 10: 0.037 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=7, k=29.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (29.75) muss ganze Zahl sein, gerundet auf 30
#>        Working on 10 repeats of 324632 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=29.75 & rep = 10: 0.034 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=7, k=30.75
#>        Input data matches to given taxa size
#> Warning in choose(m, k): 'k' (30.75) muss ganze Zahl sein, gerundet auf 31
#>        Working on 10 repeats of 52360 combinations
#>               Working on combination 1
#>               Working on combination 2
#>               Working on combination 3
#>               Working on combination 4
#>               Working on combination 5
#>               Working on combination 6
#>               Working on combination 7
#>               Working on combination 8
#>               Working on combination 9
#>               Working on combination 10
#>        Total time for n=7, k=30.75 & rep = 10: 0.029 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
SimulationResults_n07_test = rbindlist(dumTab)

SimulationResults_n07_test[,posRate := NR_FTT/NR_PhyloDec]
SimulationResults_n07_test[NR_PhyloDec==0,posRate := NA]
SimulationResults_n07_test[NR_PhyloDec>0,]
#>         k  time NR_NotPhyloDec NR_PhyloDec NR_FTT   posRate
#>  1: 20.75 0.095              9           1      1 1.0000000
#>  2: 21.75 0.093              7           3      2 0.6666667
#>  3: 22.75 0.087              5           5      4 0.8000000
#>  4: 23.75 0.069              4           6      6 1.0000000
#>  5: 24.75 0.059              3           7      7 1.0000000
#>  6: 25.75 0.058              3           7      7 1.0000000
#>  7: 26.75 0.045              1           9      9 1.0000000
#>  8: 27.75 0.044              2           8      8 1.0000000
#>  9: 28.75 0.037              1           9      9 1.0000000
#> 10: 29.75 0.034              0          10     10 1.0000000
#> 11: 30.75 0.029              0          10     10 1.0000000
```

# Loop with all combinations

``` r
dumTab = foreach(j=c(LowerBound:UpperBound))%do%{
  # j=5
  message("\nWorking on n=7, k=",j)
  time1 = Sys.time()
  myTest = SimulationFunction(number_taxa = n, 
                              number_quads = j,
                              repeats = 10000,
                              data1 = test1,
                              data2 = myTab_n07,
                              verbose = T)
  time2 = Sys.time()
  x0 = as.numeric(round(difftime(time2,time1,units = "mins"),3))
  message("       Total time for n=7, k=",j," & rep = 100: " ,round(difftime(time2,time1,units = "mins"),3)," minutes")
  
  outfn = paste0("../results/02_SimulationsData_n07/SimulationResults_n7_k",j,"_rep10000.RData")
  save(myTest,file = outfn)

  tab = table(myTest$FWPP,myTest$FTT)
  x2 = myTest[FWPP == "PhyloDec" & FTT == "FTT",.N]
  x3 = myTest[FWPP == "NOT PhyloDec",.N]
  x4 = myTest[FWPP == "PhyloDec",.N]
  
  message("       There were ",x2," of ",x4," sets identified by Fischers algorith as decisive (",round(x2/x4,4)*100,"%)")
  
  res = data.table(k = j,
                   time = x0,
                   NR_NotPhyloDec = x3,
                   NR_PhyloDec = x4,
                   NR_FTT = x2)
  res
}
SimulationResults_n07 = rbindlist(dumTab)

SimulationResults_n07[,posRate := NR_FTT/NR_PhyloDec]
SimulationResults_n07[NR_PhyloDec==0,posRate := NA]
save(SimulationResults_n07, file = "../results/02_SimulationResults_n07.RData")

SimulationResults_n07[NR_PhyloDec>0,]
```

# Session Info

``` r
sessionInfo()
#> R version 4.2.2 (2022-10-31 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 22621)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8   
#> [3] LC_MONETARY=German_Germany.utf8 LC_NUMERIC=C                   
#> [5] LC_TIME=German_Germany.utf8    
#> 
#> attached base packages:
#> [1] grid      stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] cowplot_1.1.1           gtable_0.3.1            ggplot2_3.4.1          
#> [4] FixingTaxonTraceR_0.0.1 foreach_1.5.2           data.table_1.14.8      
#> 
#> loaded via a namespace (and not attached):
#>  [1] rstudioapi_0.14  knitr_1.42       magrittr_2.0.3   munsell_0.5.0   
#>  [5] colorspace_2.1-0 R6_2.5.1         rlang_1.0.6      fastmap_1.1.1   
#>  [9] fansi_1.0.4      tools_4.2.2      xfun_0.37        utf8_1.2.3      
#> [13] cli_3.6.0        withr_2.5.0      htmltools_0.5.4  iterators_1.0.14
#> [17] yaml_2.3.7       digest_0.6.31    tibble_3.2.0     lifecycle_1.0.3 
#> [21] vctrs_0.5.2      codetools_0.2-18 glue_1.6.2       evaluate_0.20   
#> [25] rmarkdown_2.20   compiler_4.2.2   pillar_1.8.1     scales_1.2.1    
#> [29] pkgconfig_2.0.3
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
#> 
#> TOTAL TIME : 2.625 minutes
```
