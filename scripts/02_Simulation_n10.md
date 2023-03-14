
<!-- 02_Simulation_n10.md is generated from 02_Simulation_n10.Rmd. Please edit that file -->

# Introduction

<!-- badges: start -->
<!-- badges: end -->

I want to run a simulation for n=10 taxa. There are 70 possible
quadruples with n=10 taxa, and I want to test 100000 possible
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

In *myTab_n10*, there are all possible Four-way-partitions (4WPP) given
n=10 taxa.

``` r
# Load initial input data
test1 = FTT_createInput(fn="../data/S10_Decisive.txt",sepSym = "_",c = 4)
#> Input contains 7 sets with 10 different taxa. 
#> The largest set has 9 taxa.

# Sanity check of FTT package
test_FTT = FTT_algorithmRed(data = test1$data,verbose = T, c=4, n=10)
#> Using 196 of 210 4-tuples as input for algorithm (10 unique taxa). 
#>  This leaves 14 4-tuples unsolved.
#> In round #1, 14 4-tuples could be resolved ...
#> FIXING TAXON TRACEABLE
#>  It follows that the set is phylogenetically decisive
#test_NRC = FTT_findNRC(data = test1)

load("../results/01_partitions/partitions_n10.RData")
```

# Test-Loop with 10 combinations per k

To test less $k$, I use the more strict lower bound for touple numbers
(see FTT R Package for more explanation).

``` r
test1$data[,status:=NA]
n=10
LowerBound = (1/6) * 4 * choose(n,2)
UpperBound = choose(n,4) - n + 3

dumTab = foreach(j=c(LowerBound:UpperBound))%do%{
  # j=40
  message("\nWorking on n=10, k=",j)
  time1 = Sys.time()
  myTest = SimulationFunction(number_taxa = n, 
                              number_quads = j,
                              repeats = 10,
                              data1 = test1,
                              data2 = myTab_n10,
                              verbose = F,
                              FFT_only_if_PhyloDec = T)
  time2 = Sys.time()
  x0 = as.numeric(round(difftime(time2,time1,units = "mins"),3))
  message("       Total time for n=10, k=",j," & rep = 10: " ,round(difftime(time2,time1,units = "mins"),3)," minutes")
  
  outfn = paste0("../temp/02_SimulationsData_n10/SimulationResults_n10_k",j,".RData")
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
#> Working on n=10, k=30
#>        Working on 10 repeats of 1.98587137236679e+36 combinations
#>        Total time for n=10, k=30 & rep = 10: 0.499 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=31
#>        Working on 10 repeats of 1.15308660330975e+37 combinations
#>        Total time for n=10, k=31 & rep = 10: 0.51 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=32
#>        Working on 10 repeats of 6.4500781872639e+37 combinations
#>        Total time for n=10, k=32 & rep = 10: 0.53 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=33
#>        Working on 10 repeats of 3.47913308282726e+38 combinations
#>        Total time for n=10, k=33 & rep = 10: 0.563 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=34
#>        Working on 10 repeats of 1.81119575194243e+39 combinations
#>        Total time for n=10, k=34 & rep = 10: 0.575 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=35
#>        Working on 10 repeats of 9.10772720976741e+39 combinations
#>        Total time for n=10, k=35 & rep = 10: 0.586 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=36
#>        Working on 10 repeats of 4.427367393637e+40 combinations
#>        Total time for n=10, k=36 & rep = 10: 0.595 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=37
#>        Working on 10 repeats of 2.08205926079142e+41 combinations
#>        Total time for n=10, k=37 & rep = 10: 0.628 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=38
#>        Working on 10 repeats of 9.47884873991883e+41 combinations
#>        Total time for n=10, k=38 & rep = 10: 0.638 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=39
#>        Working on 10 repeats of 4.18041534170788e+42 combinations
#>        Total time for n=10, k=39 & rep = 10: 0.656 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=40
#>        Working on 10 repeats of 1.78712755858008e+43 combinations
#>        Total time for n=10, k=40 & rep = 10: 0.672 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=41
#>        Working on 10 repeats of 7.41004109655167e+43 combinations
#>        Total time for n=10, k=41 & rep = 10: 0.691 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=42
#>        Working on 10 repeats of 2.98165939361247e+44 combinations
#>        Total time for n=10, k=42 & rep = 10: 0.721 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=43
#>        Working on 10 repeats of 1.16492739099279e+45 combinations
#>        Total time for n=10, k=43 & rep = 10: 0.713 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=44
#>        Working on 10 repeats of 4.42142896126796e+45 combinations
#>        Total time for n=10, k=44 & rep = 10: 0.717 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=45
#>        Working on 10 repeats of 1.6310160168233e+46 combinations
#>        Total time for n=10, k=45 & rep = 10: 0.732 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=46
#>        Working on 10 repeats of 5.85038353860539e+46 combinations
#>        Total time for n=10, k=46 & rep = 10: 0.805 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=47
#>        Working on 10 repeats of 2.04141042623676e+47 combinations
#>        Total time for n=10, k=47 & rep = 10: 0.769 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=48
#>        Working on 10 repeats of 6.93228957242899e+47 combinations
#>        Total time for n=10, k=48 & rep = 10: 0.814 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=49
#>        Working on 10 repeats of 2.29189981782347e+48 combinations
#>        Total time for n=10, k=49 & rep = 10: 0.815 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=50
#>        Working on 10 repeats of 7.37991741339163e+48 combinations
#>        Total time for n=10, k=50 & rep = 10: 0.823 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=51
#>        Working on 10 repeats of 2.31526820812283e+49 combinations
#>        Total time for n=10, k=51 & rep = 10: 0.839 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=52
#>        Working on 10 repeats of 7.0793777902217e+49 combinations
#>        Total time for n=10, k=52 & rep = 10: 0.859 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=53
#>        Working on 10 repeats of 2.11045602048119e+50 combinations
#>        Total time for n=10, k=53 & rep = 10: 0.87 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=54
#>        Working on 10 repeats of 6.13595546695463e+50 combinations
#>        Total time for n=10, k=54 & rep = 10: 0.888 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=55
#>        Working on 10 repeats of 1.7403800960817e+51 combinations
#>        Total time for n=10, k=55 & rep = 10: 0.914 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=56
#>        Working on 10 repeats of 4.81712348022608e+51 combinations
#>        Total time for n=10, k=56 & rep = 10: 0.935 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=57
#>        Working on 10 repeats of 1.30146844904354e+52 combinations
#>        Total time for n=10, k=57 & rep = 10: 0.958 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=58
#>        Working on 10 repeats of 3.43318401213212e+52 combinations
#>        Total time for n=10, k=58 & rep = 10: 0.973 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=59
#>        Working on 10 repeats of 8.8448130482046e+52 combinations
#>        Total time for n=10, k=59 & rep = 10: 0.999 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=60
#>        Working on 10 repeats of 2.22594461713154e+53 combinations
#>        Total time for n=10, k=60 & rep = 10: 1.001 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=61
#>        Working on 10 repeats of 5.47363430442175e+53 combinations
#>        Total time for n=10, k=61 & rep = 10: 1.007 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=62
#>        Working on 10 repeats of 1.31543792154652e+54 combinations
#>        Total time for n=10, k=62 & rep = 10: 1.034 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=63
#>        Working on 10 repeats of 3.09023511728388e+54 combinations
#>        Total time for n=10, k=63 & rep = 10: 1.051 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=64
#>        Working on 10 repeats of 7.09788378501144e+54 combinations
#>        Total time for n=10, k=64 & rep = 10: 1.079 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=65
#>        Working on 10 repeats of 1.59429389632563e+55 combinations
#>        Total time for n=10, k=65 & rep = 10: 1.106 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=66
#>        Working on 10 repeats of 3.50261537829118e+55 combinations
#>        Total time for n=10, k=66 & rep = 10: 1.101 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=67
#>        Working on 10 repeats of 7.52800917125258e+55 combinations
#>        Total time for n=10, k=67 & rep = 10: 1.128 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=68
#>        Working on 10 repeats of 1.58309604630754e+56 combinations
#>        Total time for n=10, k=68 & rep = 10: 1.127 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=69
#>        Working on 10 repeats of 3.25796577645906e+56 combinations
#>        Total time for n=10, k=69 & rep = 10: 1.142 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=70
#>        Working on 10 repeats of 6.56247392115308e+56 combinations
#>        Total time for n=10, k=70 & rep = 10: 1.21 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=71
#>        Working on 10 repeats of 1.29400894219921e+57 combinations
#>        Total time for n=10, k=71 & rep = 10: 1.177 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=72
#>        Working on 10 repeats of 2.49815615230128e+57 combinations
#>        Total time for n=10, k=72 & rep = 10: 1.214 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=73
#>        Working on 10 repeats of 4.72254176736412e+57 combinations
#>        Total time for n=10, k=73 & rep = 10: 1.223 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=74
#>        Working on 10 repeats of 8.74308408282277e+57 combinations
#>        Total time for n=10, k=74 & rep = 10: 1.211 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=75
#>        Working on 10 repeats of 1.58541258035182e+58 combinations
#>        Total time for n=10, k=75 & rep = 10: 1.226 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=76
#>        Working on 10 repeats of 2.81619339930916e+58 combinations
#>        Total time for n=10, k=76 & rep = 10: 1.242 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=77
#>        Working on 10 repeats of 4.90090799360302e+58 combinations
#>        Total time for n=10, k=77 & rep = 10: 1.293 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=78
#>        Working on 10 repeats of 8.35667645063076e+58 combinations
#>        Total time for n=10, k=78 & rep = 10: 1.332 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=79
#>        Working on 10 repeats of 1.39630543225731e+59 combinations
#>        Total time for n=10, k=79 & rep = 10: 1.328 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=80
#>        Working on 10 repeats of 2.28645014532135e+59 combinations
#>        Total time for n=10, k=80 & rep = 10: 1.326 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=81
#>        Working on 10 repeats of 3.66961134434288e+59 combinations
#>        Total time for n=10, k=81 & rep = 10: 1.344 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=82
#>        Working on 10 repeats of 5.77292516366116e+59 combinations
#>        Total time for n=10, k=82 & rep = 10: 1.363 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=83
#>        Working on 10 repeats of 8.90282434877879e+59 combinations
#>        Total time for n=10, k=83 & rep = 10: 1.377 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=84
#>        Working on 10 repeats of 1.34602225273205e+60 combinations
#>        Total time for n=10, k=84 & rep = 10: 1.386 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=85
#>        Working on 10 repeats of 1.99528004522631e+60 combinations
#>        Total time for n=10, k=85 & rep = 10: 1.401 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=86
#>        Working on 10 repeats of 2.90011634480568e+60 combinations
#>        Total time for n=10, k=86 & rep = 10: 1.419 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=87
#>        Working on 10 repeats of 4.13349915811388e+60 combinations
#>        Total time for n=10, k=87 & rep = 10: 1.429 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=88
#>        Working on 10 repeats of 5.77750450509094e+60 combinations
#>        Total time for n=10, k=88 & rep = 10: 1.446 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=89
#>        Working on 10 repeats of 7.91972527664159e+60 combinations
#>        Total time for n=10, k=89 & rep = 10: 1.45 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=90
#>        Working on 10 repeats of 1.06476306497071e+61 combinations
#>        Total time for n=10, k=90 & rep = 10: 1.482 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=91
#>        Working on 10 repeats of 1.40408316259872e+61 combinations
#>        Total time for n=10, k=91 & rep = 10: 1.485 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=92
#>        Working on 10 repeats of 1.81615104727443e+61 combinations
#>        Total time for n=10, k=92 & rep = 10: 1.564 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=93
#>        Working on 10 repeats of 2.3043636943912e+61 combinations
#>        Total time for n=10, k=93 & rep = 10: 1.734 minutes
#>        There were 2 of 2 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=94
#>        Working on 10 repeats of 2.86819736429546e+61 combinations
#>        Total time for n=10, k=94 & rep = 10: 1.639 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=95
#>        Working on 10 repeats of 3.50221993956075e+61 combinations
#>        Total time for n=10, k=95 & rep = 10: 1.634 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=96
#>        Working on 10 repeats of 4.1953676359321e+61 combinations
#>        Total time for n=10, k=96 & rep = 10: 1.678 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=97
#>        Working on 10 repeats of 4.93063825253891e+61 combinations
#>        Total time for n=10, k=97 & rep = 10: 1.627 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=98
#>        Working on 10 repeats of 5.68532778098881e+61 combinations
#>        Total time for n=10, k=98 & rep = 10: 1.658 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=99
#>        Working on 10 repeats of 6.43188597445179e+61 combinations
#>        Total time for n=10, k=99 & rep = 10: 1.622 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=100
#>        Working on 10 repeats of 7.13939343164141e+61 combinations
#>        Total time for n=10, k=100 & rep = 10: 1.671 minutes
#>        There were 0 of 1 sets identified by Fischers algorith as decisive (0%)
#> 
#> Working on n=10, k=101
#>        Working on 10 repeats of 7.77557700475797e+61 combinations
#>        Total time for n=10, k=101 & rep = 10: 1.704 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=102
#>        Working on 10 repeats of 8.30919503449644e+61 combinations
#>        Total time for n=10, k=102 & rep = 10: 1.79 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=103
#>        Working on 10 repeats of 8.71255401675361e+61 combinations
#>        Total time for n=10, k=103 & rep = 10: 1.769 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=104
#>        Working on 10 repeats of 8.96387769031386e+61 combinations
#>        Total time for n=10, k=104 & rep = 10: 1.812 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=105
#>        Working on 10 repeats of 9.04924795403077e+61 combinations
#>        Total time for n=10, k=105 & rep = 10: 1.885 minutes
#>        There were 2 of 2 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=106
#>        Working on 10 repeats of 8.96387769031386e+61 combinations
#>        Total time for n=10, k=106 & rep = 10: 1.89 minutes
#>        There were 2 of 2 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=107
#>        Working on 10 repeats of 8.71255401675361e+61 combinations
#>        Total time for n=10, k=107 & rep = 10: 1.948 minutes
#>        There were 3 of 3 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=108
#>        Working on 10 repeats of 8.30919503449644e+61 combinations
#>        Total time for n=10, k=108 & rep = 10: 1.869 minutes
#>        There were 3 of 3 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=109
#>        Working on 10 repeats of 7.77557700475797e+61 combinations
#>        Total time for n=10, k=109 & rep = 10: 1.918 minutes
#>        There were 4 of 4 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=110
#>        Working on 10 repeats of 7.13939343164141e+61 combinations
#>        Total time for n=10, k=110 & rep = 10: 1.988 minutes
#>        There were 5 of 5 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=111
#>        Working on 10 repeats of 6.43188597445179e+61 combinations
#>        Total time for n=10, k=111 & rep = 10: 2.115 minutes
#>        There were 6 of 6 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=112
#>        Working on 10 repeats of 5.68532778098881e+61 combinations
#>        Total time for n=10, k=112 & rep = 10: 2.043 minutes
#>        There were 5 of 5 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=113
#>        Working on 10 repeats of 4.93063825253891e+61 combinations
#>        Total time for n=10, k=113 & rep = 10: 2.214 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=114
#>        Working on 10 repeats of 4.1953676359321e+61 combinations
#>        Total time for n=10, k=114 & rep = 10: 2.14 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=115
#>        Working on 10 repeats of 3.50221993956075e+61 combinations
#>        Total time for n=10, k=115 & rep = 10: 2.15 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=116
#>        Working on 10 repeats of 2.86819736429546e+61 combinations
#>        Total time for n=10, k=116 & rep = 10: 2.09 minutes
#>        There were 5 of 5 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=117
#>        Working on 10 repeats of 2.3043636943912e+61 combinations
#>        Total time for n=10, k=117 & rep = 10: 2.123 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=118
#>        Working on 10 repeats of 1.81615104727443e+61 combinations
#>        Total time for n=10, k=118 & rep = 10: 2.157 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=119
#>        Working on 10 repeats of 1.40408316259872e+61 combinations
#>        Total time for n=10, k=119 & rep = 10: 2.143 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=120
#>        Working on 10 repeats of 1.06476306497071e+61 combinations
#>        Total time for n=10, k=120 & rep = 10: 2.178 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=121
#>        Working on 10 repeats of 7.91972527664159e+60 combinations
#>        Total time for n=10, k=121 & rep = 10: 2.189 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=122
#>        Working on 10 repeats of 5.77750450509094e+60 combinations
#>        Total time for n=10, k=122 & rep = 10: 2.199 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=123
#>        Working on 10 repeats of 4.13349915811388e+60 combinations
#>        Total time for n=10, k=123 & rep = 10: 2.176 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=124
#>        Working on 10 repeats of 2.90011634480568e+60 combinations
#>        Total time for n=10, k=124 & rep = 10: 2.269 minutes
#>        There were 6 of 6 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=125
#>        Working on 10 repeats of 1.99528004522631e+60 combinations
#>        Total time for n=10, k=125 & rep = 10: 2.229 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=126
#>        Working on 10 repeats of 1.34602225273205e+60 combinations
#>        Total time for n=10, k=126 & rep = 10: 2.288 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=127
#>        Working on 10 repeats of 8.90282434877879e+59 combinations
#>        Total time for n=10, k=127 & rep = 10: 2.226 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=128
#>        Working on 10 repeats of 5.77292516366116e+59 combinations
#>        Total time for n=10, k=128 & rep = 10: 2.357 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=129
#>        Working on 10 repeats of 3.66961134434288e+59 combinations
#>        Total time for n=10, k=129 & rep = 10: 2.289 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=130
#>        Working on 10 repeats of 2.28645014532135e+59 combinations
#>        Total time for n=10, k=130 & rep = 10: 2.293 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=131
#>        Working on 10 repeats of 1.39630543225731e+59 combinations
#>        Total time for n=10, k=131 & rep = 10: 2.329 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=132
#>        Working on 10 repeats of 8.35667645063076e+58 combinations
#>        Total time for n=10, k=132 & rep = 10: 2.386 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=133
#>        Working on 10 repeats of 4.90090799360302e+58 combinations
#>        Total time for n=10, k=133 & rep = 10: 2.376 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=134
#>        Working on 10 repeats of 2.81619339930916e+58 combinations
#>        Total time for n=10, k=134 & rep = 10: 2.378 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=135
#>        Working on 10 repeats of 1.58541258035182e+58 combinations
#>        Total time for n=10, k=135 & rep = 10: 2.373 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=136
#>        Working on 10 repeats of 8.74308408282277e+57 combinations
#>        Total time for n=10, k=136 & rep = 10: 2.403 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=137
#>        Working on 10 repeats of 4.72254176736412e+57 combinations
#>        Total time for n=10, k=137 & rep = 10: 2.455 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=138
#>        Working on 10 repeats of 2.49815615230128e+57 combinations
#>        Total time for n=10, k=138 & rep = 10: 2.468 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=139
#>        Working on 10 repeats of 1.29400894219921e+57 combinations
#>        Total time for n=10, k=139 & rep = 10: 3.111 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=140
#>        Working on 10 repeats of 6.56247392115308e+56 combinations
#>        Total time for n=10, k=140 & rep = 10: 5.154 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=141
#>        Working on 10 repeats of 3.25796577645906e+56 combinations
#>        Total time for n=10, k=141 & rep = 10: 5.297 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=142
#>        Working on 10 repeats of 1.58309604630754e+56 combinations
#>        Total time for n=10, k=142 & rep = 10: 5.283 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=143
#>        Working on 10 repeats of 7.52800917125258e+55 combinations
#>        Total time for n=10, k=143 & rep = 10: 5.225 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=144
#>        Working on 10 repeats of 3.50261537829118e+55 combinations
#>        Total time for n=10, k=144 & rep = 10: 5.26 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=145
#>        Working on 10 repeats of 1.59429389632563e+55 combinations
#>        Total time for n=10, k=145 & rep = 10: 5.372 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=146
#>        Working on 10 repeats of 7.09788378501144e+54 combinations
#>        Total time for n=10, k=146 & rep = 10: 5.483 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=147
#>        Working on 10 repeats of 3.09023511728388e+54 combinations
#>        Total time for n=10, k=147 & rep = 10: 5.568 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=148
#>        Working on 10 repeats of 1.31543792154652e+54 combinations
#>        Total time for n=10, k=148 & rep = 10: 5.476 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=149
#>        Working on 10 repeats of 5.47363430442175e+53 combinations
#>        Total time for n=10, k=149 & rep = 10: 5.508 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=150
#>        Working on 10 repeats of 2.22594461713154e+53 combinations
#>        Total time for n=10, k=150 & rep = 10: 5.547 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=151
#>        Working on 10 repeats of 8.8448130482046e+52 combinations
#>        Total time for n=10, k=151 & rep = 10: 5.568 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=152
#>        Working on 10 repeats of 3.43318401213212e+52 combinations
#>        Total time for n=10, k=152 & rep = 10: 5.455 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=153
#>        Working on 10 repeats of 1.30146844904354e+52 combinations
#>        Total time for n=10, k=153 & rep = 10: 5.55 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=154
#>        Working on 10 repeats of 4.81712348022608e+51 combinations
#>        Total time for n=10, k=154 & rep = 10: 5.685 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=155
#>        Working on 10 repeats of 1.7403800960817e+51 combinations
#>        Total time for n=10, k=155 & rep = 10: 3.428 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=156
#>        Working on 10 repeats of 6.13595546695463e+50 combinations
#>        Total time for n=10, k=156 & rep = 10: 2.726 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=157
#>        Working on 10 repeats of 2.11045602048119e+50 combinations
#>        Total time for n=10, k=157 & rep = 10: 2.996 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=158
#>        Working on 10 repeats of 7.0793777902217e+49 combinations
#>        Total time for n=10, k=158 & rep = 10: 2.822 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=159
#>        Working on 10 repeats of 2.31526820812283e+49 combinations
#>        Total time for n=10, k=159 & rep = 10: 2.813 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=160
#>        Working on 10 repeats of 7.37991741339163e+48 combinations
#>        Total time for n=10, k=160 & rep = 10: 2.821 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=161
#>        Working on 10 repeats of 2.29189981782347e+48 combinations
#>        Total time for n=10, k=161 & rep = 10: 2.826 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=162
#>        Working on 10 repeats of 6.93228957242899e+47 combinations
#>        Total time for n=10, k=162 & rep = 10: 2.851 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=163
#>        Working on 10 repeats of 2.04141042623676e+47 combinations
#>        Total time for n=10, k=163 & rep = 10: 2.873 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=164
#>        Working on 10 repeats of 5.85038353860539e+46 combinations
#>        Total time for n=10, k=164 & rep = 10: 2.782 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=165
#>        Working on 10 repeats of 1.6310160168233e+46 combinations
#>        Total time for n=10, k=165 & rep = 10: 2.806 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=166
#>        Working on 10 repeats of 4.42142896126796e+45 combinations
#>        Total time for n=10, k=166 & rep = 10: 2.792 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=167
#>        Working on 10 repeats of 1.16492739099279e+45 combinations
#>        Total time for n=10, k=167 & rep = 10: 2.798 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=168
#>        Working on 10 repeats of 2.98165939361247e+44 combinations
#>        Total time for n=10, k=168 & rep = 10: 2.815 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=169
#>        Working on 10 repeats of 7.41004109655167e+43 combinations
#>        Total time for n=10, k=169 & rep = 10: 2.96 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=170
#>        Working on 10 repeats of 1.78712755858008e+43 combinations
#>        Total time for n=10, k=170 & rep = 10: 2.845 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=171
#>        Working on 10 repeats of 4.18041534170788e+42 combinations
#>        Total time for n=10, k=171 & rep = 10: 2.915 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=172
#>        Working on 10 repeats of 9.47884873991883e+41 combinations
#>        Total time for n=10, k=172 & rep = 10: 2.887 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=173
#>        Working on 10 repeats of 2.08205926079142e+41 combinations
#>        Total time for n=10, k=173 & rep = 10: 2.89 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=174
#>        Working on 10 repeats of 4.427367393637e+40 combinations
#>        Total time for n=10, k=174 & rep = 10: 2.903 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=175
#>        Working on 10 repeats of 9.10772720976741e+39 combinations
#>        Total time for n=10, k=175 & rep = 10: 2.909 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=176
#>        Working on 10 repeats of 1.81119575194243e+39 combinations
#>        Total time for n=10, k=176 & rep = 10: 2.963 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=177
#>        Working on 10 repeats of 3.47913308282726e+38 combinations
#>        Total time for n=10, k=177 & rep = 10: 2.939 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=178
#>        Working on 10 repeats of 6.4500781872639e+37 combinations
#>        Total time for n=10, k=178 & rep = 10: 2.958 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=179
#>        Working on 10 repeats of 1.15308660330975e+37 combinations
#>        Total time for n=10, k=179 & rep = 10: 2.97 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=180
#>        Working on 10 repeats of 1.98587137236679e+36 combinations
#>        Total time for n=10, k=180 & rep = 10: 3.076 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=181
#>        Working on 10 repeats of 3.29149951221018e+35 combinations
#>        Total time for n=10, k=181 & rep = 10: 2.998 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=182
#>        Working on 10 repeats of 5.24469702495028e+34 combinations
#>        Total time for n=10, k=182 & rep = 10: 3.004 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=183
#>        Working on 10 repeats of 8.02467304363978e+33 combinations
#>        Total time for n=10, k=183 & rep = 10: 3.018 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=184
#>        Working on 10 repeats of 1.17753354444714e+33 combinations
#>        Total time for n=10, k=184 & rep = 10: 3.035 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=185
#>        Working on 10 repeats of 1.6549120084122e+32 combinations
#>        Total time for n=10, k=185 & rep = 10: 3.042 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=186
#>        Working on 10 repeats of 2.22434409732822e+31 combinations
#>        Total time for n=10, k=186 & rep = 10: 3.099 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=187
#>        Working on 10 repeats of 2.85477317304157e+30 combinations
#>        Total time for n=10, k=187 & rep = 10: 3.097 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=188
#>        Working on 10 repeats of 3.49254164787001e+29 combinations
#>        Total time for n=10, k=188 & rep = 10: 3.096 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=189
#>        Working on 10 repeats of 4.06539239434604e+28 combinations
#>        Total time for n=10, k=189 & rep = 10: 3.119 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=190
#>        Working on 10 repeats of 4.49332843585615e+27 combinations
#>        Total time for n=10, k=190 & rep = 10: 3.128 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=191
#>        Working on 10 repeats of 4.70505595377607e+26 combinations
#>        Total time for n=10, k=191 & rep = 10: 3.146 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=192
#>        Working on 10 repeats of 4.65604495425757e+25 combinations
#>        Total time for n=10, k=192 & rep = 10: 3.183 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=193
#>        Working on 10 repeats of 4.3424253459397e+24 combinations
#>        Total time for n=10, k=193 & rep = 10: 3.168 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=194
#>        Working on 10 repeats of 3.80521808664819e+23 combinations
#>        Total time for n=10, k=194 & rep = 10: 3.185 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=195
#>        Working on 10 repeats of 3.12223022494211e+22 combinations
#>        Total time for n=10, k=195 & rep = 10: 3.184 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=196
#>        Working on 10 repeats of 2.38946190684345e+21 combinations
#>        Total time for n=10, k=196 & rep = 10: 3.225 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=197
#>        Working on 10 repeats of 1.6980947561324e+20 combinations
#>        Total time for n=10, k=197 & rep = 10: 3.219 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=198
#>        Working on 10 repeats of 11149106984707682304 combinations
#>        Total time for n=10, k=198 & rep = 10: 3.237 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=199
#>        Working on 10 repeats of 672307958876845184 combinations
#>        Total time for n=10, k=199 & rep = 10: 3.241 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=200
#>        Working on 10 repeats of 36976937738226480 combinations
#>        Total time for n=10, k=200 & rep = 10: 3.272 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=201
#>        Working on 10 repeats of 1839648643692860 combinations
#>        Total time for n=10, k=201 & rep = 10: 3.388 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=202
#>        Working on 10 repeats of 81964543530870 combinations
#>        Total time for n=10, k=202 & rep = 10: 3.301 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=203
#>        Working on 10 repeats of 3230129794320 combinations
#>        Total time for n=10, k=203 & rep = 10: 3.315 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
SimulationResults_n10_test = rbindlist(dumTab)

SimulationResults_n10_test[,posRate := NR_FTT/NR_PhyloDec]
SimulationResults_n10_test[NR_PhyloDec==0,posRate := NA]
SimulationResults_n10_test[NR_PhyloDec>0,]
#>        k  time NR_NotPhyloDec NR_PhyloDec NR_FTT posRate
#>   1:  92 1.564              9           1      1       1
#>   2:  93 1.734              8           2      2       1
#>   3:  96 1.678              9           1      1       1
#>   4:  97 1.627              9           1      1       1
#>   5:  98 1.658              9           1      1       1
#>  ---                                                    
#> 103: 199 3.241              0          10     10       1
#> 104: 200 3.272              0          10     10       1
#> 105: 201 3.388              0          10     10       1
#> 106: 202 3.301              0          10     10       1
#> 107: 203 3.315              0          10     10       1
```

# Loop with all combinations

``` r
dumTab = foreach(j=c(LowerBound:UpperBound))%do%{
  # j=5
  message("\nWorking on n=10, k=",j)
  time1 = Sys.time()
  myTest = SimulationFunction(number_taxa = n, 
                              number_quads = j,
                              repeats = 10000,
                              data1 = test1,
                              data2 = myTab_n10,
                              verbose = T,
                              FFT_only_if_PhyloDec = T)
  time2 = Sys.time()
  x0 = as.numeric(round(difftime(time2,time1,units = "mins"),3))
  message("       Total time for n=10, k=",j," & rep = 100: " ,round(difftime(time2,time1,units = "mins"),3)," minutes")
  
  outfn = paste0("../results/02_SimulationsData_n10/SimulationResults_n10_k",k,"_rep10000.RData")
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
SimulationResults_n10 = rbindlist(dumTab)

SimulationResults_n10[,posRate := NR_FTT/NR_PhyloDec]
SimulationResults_n10[NR_PhyloDec==0,posRate := NA]
save(SimulationResults_n10, file = "../results/02_SimulationResults_n10.RData")

SimulationResults_n10[NR_PhyloDec>0,]
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
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] FixingTaxonTraceR_0.0.1 foreach_1.5.2           data.table_1.14.8      
#> 
#> loaded via a namespace (and not attached):
#>  [1] codetools_0.2-18 digest_0.6.31    magrittr_2.0.3   evaluate_0.20   
#>  [5] rlang_1.0.6      cli_3.6.0        rstudioapi_0.14  rmarkdown_2.20  
#>  [9] iterators_1.0.14 tools_4.2.2      xfun_0.37        yaml_2.3.7      
#> [13] fastmap_1.1.1    compiler_4.2.2   htmltools_0.5.4  knitr_1.42
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
#> 
#> TOTAL TIME : 391.556 minutes
```
