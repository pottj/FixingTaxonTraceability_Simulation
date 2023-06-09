
<!-- 02_Simulation_n08.md is generated from 02_Simulation_n08.Rmd. Please edit that file -->

# Introduction

<!-- badges: start -->
<!-- badges: end -->

I want to run a simulation for n=8 taxa. There are 70 possible
quadruples with n=8 taxa, and I want to test 100000 possible
combinations for each k out of 70 quadruples for phylogenetic
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

In *myTab_n8*, there are all possible Four-way-partitions (4WPP) given
n=8 taxa.

``` r
# Load initial input data
test1 = FTT_createInput(fn="../data/S8_Decisive.txt",sepSym = "_",c = 4)
#> Input contains 28 sets with 8 different taxa. 
#> The largest set has 4 taxa.

# Sanity check of FTT package
test_FTT = FTT_algorithmRed(data = test1$data,verbose = T, c=4, n=8)
#> Using 28 of 70 4-tuples as input for algorithm (8 unique taxa). 
#>  This leaves 42 4-tuples unsolved.
#> In round #1, 0 4-tuples could be resolved ...
#> NOT RESOLVABLE VIA THIS ALGORITHM! 
#>  There are 42 remaining cross 4-tuples.
test_NRC = FTT_findNRC(data = test1)
#> THERE IS ONLY RAINBOW COLORING POSSIBLE OR 3-WAY COLORING
#>  It follows that the set is phylogenetically decisive

load("../results/01_partitions/partitions_n08.RData")
```

# Test-Loop with 10 combinations per k

To test less $k$, I use as lower bound the minimal triple covering,
$\frac{1}{4}\binom{n}{3}$. As upper bound I use
$\binom{n}{4} - (n-4) -1$, as above this $k$ all sets are phylogenetic
decisive.

``` r
test1$data[,status:=NA]
n=8
LowerBound = choose(n,3)/4
UpperBound = choose(n,4) - n + 3

dumTab = foreach(j=c(LowerBound:UpperBound))%do%{
  # j=20
  message("\nWorking on n=8, k=",j)
  time1 = Sys.time()
  myTest = SimulationFunction(number_taxa = n, 
                                    number_quads = j,
                                    repeats = 10,
                                    data1 = test1,
                                    data2 = myTab_n08,
                                    verbose = F)
  time2 = Sys.time()
  x0 = as.numeric(round(difftime(time2,time1,units = "mins"),3))
  message("       Total time for n=8, k=",j," & rep = 10: " ,round(difftime(time2,time1,units = "mins"),3)," minutes")
  
  outfn = paste0("../temp/02_SimulationsData_n08/SimulationResults_n8_k",j,".RData")
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
#> Working on n=8, k=14
#>        Working on 10 repeats of 193253756909160 combinations
#>        Total time for n=8, k=14 & rep = 10: 0.189 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=15
#>        Working on 10 repeats of 721480692460864 combinations
#>        Total time for n=8, k=15 & rep = 10: 0.183 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=16
#>        Working on 10 repeats of 2480089880334220 combinations
#>        Total time for n=8, k=16 & rep = 10: 0.195 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=17
#>        Working on 10 repeats of 7877932561061640 combinations
#>        Total time for n=8, k=17 & rep = 10: 0.238 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=18
#>        Working on 10 repeats of 23196134763125940 combinations
#>        Total time for n=8, k=18 & rep = 10: 0.194 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=19
#>        Working on 10 repeats of 63484158299081520 combinations
#>        Total time for n=8, k=19 & rep = 10: 0.196 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=20
#>        Working on 10 repeats of 161884603662657856 combinations
#>        Total time for n=8, k=20 & rep = 10: 0.203 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=21
#>        Working on 10 repeats of 385439532530137728 combinations
#>        Total time for n=8, k=21 & rep = 10: 0.247 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=22
#>        Working on 10 repeats of 858478958817124864 combinations
#>        Total time for n=8, k=22 & rep = 10: 0.238 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=23
#>        Working on 10 repeats of 1791608261879217152 combinations
#>        Total time for n=8, k=23 & rep = 10: 0.231 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=24
#>        Working on 10 repeats of 3508566179513466880 combinations
#>        Total time for n=8, k=24 & rep = 10: 0.212 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=25
#>        Working on 10 repeats of 6455761770304779264 combinations
#>        Total time for n=8, k=25 & rep = 10: 0.253 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=26
#>        Working on 10 repeats of 11173433833219811328 combinations
#>        Total time for n=8, k=26 & rep = 10: 0.264 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=27
#>        Working on 10 repeats of 18208558839321174016 combinations
#>        Total time for n=8, k=27 & rep = 10: 0.346 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=28
#>        Working on 10 repeats of 27963143931814662144 combinations
#>        Total time for n=8, k=28 & rep = 10: 0.291 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=29
#>        Working on 10 repeats of 40498346384007438336 combinations
#>        Total time for n=8, k=29 & rep = 10: 0.297 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=30
#>        Working on 10 repeats of 55347740058143899658 combinations
#>        Total time for n=8, k=30 & rep = 10: 0.322 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=31
#>        Working on 10 repeats of 71416438784701202432 combinations
#>        Total time for n=8, k=31 & rep = 10: 0.309 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=32
#>        Working on 10 repeats of 87038784768855392256 combinations
#>        Total time for n=8, k=32 & rep = 10: 0.329 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=33
#>        Working on 10 repeats of 1.00226479430802e+20 combinations
#>        Total time for n=8, k=33 & rep = 10: 0.271 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=34
#>        Working on 10 repeats of 1.09069992321756e+20 combinations
#>        Total time for n=8, k=34 & rep = 10: 0.305 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=8, k=35
#>        Working on 10 repeats of 1.12186277816662e+20 combinations
#>        Total time for n=8, k=35 & rep = 10: 0.295 minutes
#>        There were 1 of 2 sets identified by Fischers algorith as decisive (50%)
#> 
#> Working on n=8, k=36
#>        Working on 10 repeats of 1.09069992321756e+20 combinations
#>        Total time for n=8, k=36 & rep = 10: 0.294 minutes
#>        There were 0 of 1 sets identified by Fischers algorith as decisive (0%)
#> 
#> Working on n=8, k=37
#>        Working on 10 repeats of 1.00226479430802e+20 combinations
#>        Total time for n=8, k=37 & rep = 10: 0.255 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=38
#>        Working on 10 repeats of 87038784768855392256 combinations
#>        Total time for n=8, k=38 & rep = 10: 0.232 minutes
#>        There were 2 of 2 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=39
#>        Working on 10 repeats of 71416438784701202432 combinations
#>        Total time for n=8, k=39 & rep = 10: 0.209 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=40
#>        Working on 10 repeats of 55347740058143899658 combinations
#>        Total time for n=8, k=40 & rep = 10: 0.184 minutes
#>        There were 4 of 4 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=41
#>        Working on 10 repeats of 40498346384007438336 combinations
#>        Total time for n=8, k=41 & rep = 10: 0.164 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=42
#>        Working on 10 repeats of 27963143931814662144 combinations
#>        Total time for n=8, k=42 & rep = 10: 0.204 minutes
#>        There were 5 of 6 sets identified by Fischers algorith as decisive (83.33%)
#> 
#> Working on n=8, k=43
#>        Working on 10 repeats of 18208558839321174016 combinations
#>        Total time for n=8, k=43 & rep = 10: 0.148 minutes
#>        There were 5 of 5 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=44
#>        Working on 10 repeats of 11173433833219811328 combinations
#>        Total time for n=8, k=44 & rep = 10: 0.133 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=45
#>        Working on 10 repeats of 6455761770304779264 combinations
#>        Total time for n=8, k=45 & rep = 10: 0.127 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=46
#>        Working on 10 repeats of 3508566179513466880 combinations
#>        Total time for n=8, k=46 & rep = 10: 0.12 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=47
#>        Working on 10 repeats of 1791608261879217152 combinations
#>        Total time for n=8, k=47 & rep = 10: 0.109 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=48
#>        Working on 10 repeats of 858478958817124864 combinations
#>        Total time for n=8, k=48 & rep = 10: 0.1 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=49
#>        Working on 10 repeats of 385439532530137728 combinations
#>        Total time for n=8, k=49 & rep = 10: 0.099 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=50
#>        Working on 10 repeats of 161884603662657856 combinations
#>        Total time for n=8, k=50 & rep = 10: 0.094 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=51
#>        Working on 10 repeats of 63484158299081520 combinations
#>        Total time for n=8, k=51 & rep = 10: 0.098 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=52
#>        Working on 10 repeats of 23196134763125940 combinations
#>        Total time for n=8, k=52 & rep = 10: 0.09 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=53
#>        Working on 10 repeats of 7877932561061640 combinations
#>        Total time for n=8, k=53 & rep = 10: 0.091 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=54
#>        Working on 10 repeats of 2480089880334220 combinations
#>        Total time for n=8, k=54 & rep = 10: 0.085 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=55
#>        Working on 10 repeats of 721480692460864 combinations
#>        Total time for n=8, k=55 & rep = 10: 0.082 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=56
#>        Working on 10 repeats of 193253756909160 combinations
#>        Total time for n=8, k=56 & rep = 10: 0.075 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=57
#>        Working on 10 repeats of 47465835030320 combinations
#>        Total time for n=8, k=57 & rep = 10: 0.08 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=58
#>        Working on 10 repeats of 10638894058520 combinations
#>        Total time for n=8, k=58 & rep = 10: 0.072 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=59
#>        Working on 10 repeats of 2163842859360 combinations
#>        Total time for n=8, k=59 & rep = 10: 0.072 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=60
#>        Working on 10 repeats of 396704524216 combinations
#>        Total time for n=8, k=60 & rep = 10: 0.067 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=61
#>        Working on 10 repeats of 65033528560 combinations
#>        Total time for n=8, k=61 & rep = 10: 0.062 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=62
#>        Working on 10 repeats of 9440350920 combinations
#>        Total time for n=8, k=62 & rep = 10: 0.061 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=63
#>        Working on 10 repeats of 1198774720 combinations
#>        Total time for n=8, k=63 & rep = 10: 0.06 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=64
#>        Working on 10 repeats of 131115985 combinations
#>        Total time for n=8, k=64 & rep = 10: 0.058 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=8, k=65
#>        Working on 10 repeats of 12103014 combinations
#>        Total time for n=8, k=65 & rep = 10: 0.061 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
SimulationResults_n08_test = rbindlist(dumTab)

SimulationResults_n08_test[,posRate := NR_FTT/NR_PhyloDec]
SimulationResults_n08_test[NR_PhyloDec==0,posRate := NA]
SimulationResults_n08_test[NR_PhyloDec>0,]
#>      k  time NR_NotPhyloDec NR_PhyloDec NR_FTT   posRate
#>  1: 35 0.295              8           2      1 0.5000000
#>  2: 36 0.294              9           1      0 0.0000000
#>  3: 37 0.255              9           1      1 1.0000000
#>  4: 38 0.232              8           2      2 1.0000000
#>  5: 39 0.209              9           1      1 1.0000000
#>  6: 40 0.184              6           4      4 1.0000000
#>  7: 41 0.164              2           8      8 1.0000000
#>  8: 42 0.204              4           6      5 0.8333333
#>  9: 43 0.148              5           5      5 1.0000000
#> 10: 44 0.133              3           7      7 1.0000000
#> 11: 45 0.127              2           8      8 1.0000000
#> 12: 46 0.120              1           9      9 1.0000000
#> 13: 47 0.109              0          10     10 1.0000000
#> 14: 48 0.100              0          10     10 1.0000000
#> 15: 49 0.099              1           9      9 1.0000000
#> 16: 50 0.094              0          10     10 1.0000000
#> 17: 51 0.098              0          10     10 1.0000000
#> 18: 52 0.090              0          10     10 1.0000000
#> 19: 53 0.091              0          10     10 1.0000000
#> 20: 54 0.085              0          10     10 1.0000000
#> 21: 55 0.082              0          10     10 1.0000000
#> 22: 56 0.075              0          10     10 1.0000000
#> 23: 57 0.080              0          10     10 1.0000000
#> 24: 58 0.072              0          10     10 1.0000000
#> 25: 59 0.072              0          10     10 1.0000000
#> 26: 60 0.067              0          10     10 1.0000000
#> 27: 61 0.062              0          10     10 1.0000000
#> 28: 62 0.061              0          10     10 1.0000000
#> 29: 63 0.060              0          10     10 1.0000000
#> 30: 64 0.058              0          10     10 1.0000000
#> 31: 65 0.061              0          10     10 1.0000000
#>      k  time NR_NotPhyloDec NR_PhyloDec NR_FTT   posRate
```

# Loop with all combinations

``` r
dumTab = foreach(j=c(LowerBound:UpperBound))%do%{
  # j=5
  message("\nWorking on n=8, k=",j)
  time1 = Sys.time()
  myTest = SimulationFunction(number_taxa = n, 
                              number_quads = j,
                              repeats = 10000,
                              data1 = test1,
                              data2 = myTab_n08,
                              verbose = T)
  time2 = Sys.time()
  x0 = as.numeric(round(difftime(time2,time1,units = "mins"),3))
  message("       Total time for n=8, k=",j," & rep = 100: " ,round(difftime(time2,time1,units = "mins"),3)," minutes")
  
  outfn = paste0("../results/02_SimulationsData_n08/SimulationResults_n8_k",j,"_rep10000.RData")
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
SimulationResults_n08 = rbindlist(dumTab)

SimulationResults_n08[,posRate := NR_FTT/NR_PhyloDec]
SimulationResults_n08[NR_PhyloDec==0,posRate := NA]
save(SimulationResults_n08, file = "../results/02_SimulationResults_n08.RData")

SimulationResults_n08[NR_PhyloDec>0,]
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
#> TOTAL TIME : 13.711 minutes
```
