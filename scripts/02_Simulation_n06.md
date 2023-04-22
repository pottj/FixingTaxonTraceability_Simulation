
<!-- 02_Simulation_n06.md is generated from 02_Simulation_n06.Rmd. Please edit that file -->

# Introduction

<!-- badges: start -->
<!-- badges: end -->

I want to run a simulation for n=6 taxa. There are 15 possible
quadruples with n=6 taxa, and I want to test all possible combinations
for k out of 15 quadruples for phylogenetic decisiveness.

For each random set of quadruples, I test both the
4-way-partition-property using the No-Rainbow-Coloring algorithm
(**NRC**, truth) and the Fixing Taxon Traceability algorithm (**FTT**,
test), both implemented in the R package **FixingTaxonTraceR**.

# Initialize

I use a file named *SourceFile.R* that contains all relevant R packages
and user-/server-specific path to the R library. If using this code, you
must make all the necessary changes within the template source file.

``` r
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
source("../helperFunctions/SimulationFunction_exact.R")
```

# Get input data

In *myTab_n6*, there are all possible Four-way-partitions (4WP) given
n=6 taxa.

In *quadruple_data*, there are all possible quadruples given n=6 taxa.

``` r
# Load initial input data
test1 = FTT_createInput(fn="../data/S6_Decisive.txt",sepSym = "_",c = 4)
#> Input contains 10 sets with 6 different taxa. 
#> The largest set has 4 taxa.

# Sanity check of FTT package
test_FTT = FTT_algorithmRed(data = test1$data,verbose = T, c=4, n=6)
#> Using 10 of 15 4-tuples as input for algorithm (6 unique taxa). 
#>  This leaves 5 4-tuples unsolved.
#> In round #1, 2 4-tuples could be resolved ...
#> In round #2, 2 4-tuples could be resolved ...
#> In round #3, 1 4-tuples could be resolved ...
#> FIXING TAXON TRACEABLE
#>  It follows that the set is phylogenetically decisive
test_NRC = FTT_findNRC(data = test1)
#> THERE IS ONLY RAINBOW COLORING POSSIBLE OR 3-WAY COLORING
#>  It follows that the set is phylogenetically decisive

load("../results/01_partitions/partitions_n06.RData")
```

# Test-Loop with 10 combinations per k

To test less $k$, I use as lower bound the minimal triple covering,
$\frac{1}{4}\binom{n}{3}$. As upper bound I use
$\binom{n}{4} - (n-4) -1$, as above this $k$ all sets are phylogenetic
decisive.

``` r
test1$data[,status:=NA]
n=6
LowerBound = choose(n,3)/4
UpperBound = choose(n,4) - n + 3

dumTab = foreach(j=c(LowerBound:UpperBound))%do%{
  # j=5
  message("\nWorking on n=6, k=",j)
  time1 = Sys.time()
  myTest = SimulationFunction_exact(number_taxa = 6, 
                                    number_quads = j,
                                    repeats = 10,
                                    data1 = test1,
                                    data2 = myTab_n06,
                                    verbose = T)
  time2 = Sys.time()
  x0 = as.numeric(round(difftime(time2,time1,units = "mins"),3))
  message("       Total time for n=6, k=",j," & rep = 10: " ,round(difftime(time2,time1,units = "mins"),3)," minutes")
  
  outfn = paste0("../temp/02_SimulationsData_n06/SimulationResults_n6_k",j,".RData")
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
#> Working on n=6, k=5
#>        Input data matches to given taxa size
#>        Working on 10 repeats of 3003 combinations
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
#>        Total time for n=6, k=5 & rep = 10: 0.047 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=6, k=6
#>        Input data matches to given taxa size
#>        Working on 10 repeats of 5005 combinations
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
#>        Total time for n=6, k=6 & rep = 10: 0.044 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=6, k=7
#>        Input data matches to given taxa size
#>        Working on 10 repeats of 6435 combinations
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
#>        Total time for n=6, k=7 & rep = 10: 0.037 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=6, k=8
#>        Input data matches to given taxa size
#>        Working on 10 repeats of 6435 combinations
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
#>        Total time for n=6, k=8 & rep = 10: 0.038 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=6, k=9
#>        Input data matches to given taxa size
#>        Working on 10 repeats of 5005 combinations
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
#>        Total time for n=6, k=9 & rep = 10: 0.042 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=6, k=10
#>        Input data matches to given taxa size
#>        Working on 10 repeats of 3003 combinations
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
#>        Total time for n=6, k=10 & rep = 10: 0.036 minutes
#>        There were 5 of 5 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=6, k=11
#>        Input data matches to given taxa size
#>        Working on 10 repeats of 1365 combinations
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
#>        Total time for n=6, k=11 & rep = 10: 0.031 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=6, k=12
#>        Input data matches to given taxa size
#>        Working on 10 repeats of 455 combinations
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
#>        Total time for n=6, k=12 & rep = 10: 0.022 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
SimulationResults_n06_test = rbindlist(dumTab)

SimulationResults_n06_test[,posRate := NR_FTT/NR_PhyloDec]
SimulationResults_n06_test[NR_PhyloDec==0,posRate := NA]
SimulationResults_n06_test[NR_PhyloDec>0,]
#>     k  time NR_NotPhyloDec NR_PhyloDec NR_FTT posRate
#> 1: 10 0.036              5           5      5       1
#> 2: 11 0.031              0          10     10       1
#> 3: 12 0.022              0          10     10       1
```

# Loop with all combinations

I use for n=6 still the exact test, as the maximal combination number is
6435 \< 10000. I test only within a reasonable range

``` r
dumTab = foreach(j=c(LowerBound:UpperBound))%do%{
  # j=5
  message("\nWorking on n=6, k=",j)
  time1 = Sys.time()
  myTest = SimulationFunction_exact(number_taxa = 6, 
                                    number_quads = j,
                                    repeats = 10000,
                                    data1 = test1,
                                    data2 = myTab_n06,
                                    verbose = T)
  time2 = Sys.time()
  x0 = as.numeric(round(difftime(time2,time1,units = "mins"),3))
  message("       Total time for n=6, k=",j," & rep = 100: " ,round(difftime(time2,time1,units = "mins"),3)," minutes")
  
  outfn = paste0("../results/02_SimulationsData_n06/SimulationResults_n6_k",j,"_rep10000.RData")
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
SimulationResults_n06 = rbindlist(dumTab)

SimulationResults_n06[,posRate := NR_FTT/NR_PhyloDec]
SimulationResults_n06[NR_PhyloDec==0,posRate := NA]
save(SimulationResults_n06, file = "../results/02_SimulationResults_n06.RData")

SimulationResults_n06[NR_PhyloDec>0,]
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
#> TOTAL TIME : 0.457 minutes
```
