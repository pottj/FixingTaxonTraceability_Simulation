
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

To test less $k$, I use as lower bound the minimal triple covering,
$\frac{1}{4}\binom{n}{3}$. As upper bound I use
$\binom{n}{4} - (n-4) -1$, as above this $k$ all sets are phylogenetic
decisive.

``` r
test1$data[,status:=NA]
n=10
LowerBound = choose(n,3)/4
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
#>        Total time for n=10, k=30 & rep = 10: 0.549 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=31
#>        Working on 10 repeats of 1.15308660330975e+37 combinations
#>        Total time for n=10, k=31 & rep = 10: 0.552 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=32
#>        Working on 10 repeats of 6.4500781872639e+37 combinations
#>        Total time for n=10, k=32 & rep = 10: 0.63 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=33
#>        Working on 10 repeats of 3.47913308282726e+38 combinations
#>        Total time for n=10, k=33 & rep = 10: 0.607 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=34
#>        Working on 10 repeats of 1.81119575194243e+39 combinations
#>        Total time for n=10, k=34 & rep = 10: 0.6 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=35
#>        Working on 10 repeats of 9.10772720976741e+39 combinations
#>        Total time for n=10, k=35 & rep = 10: 0.613 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=36
#>        Working on 10 repeats of 4.427367393637e+40 combinations
#>        Total time for n=10, k=36 & rep = 10: 0.629 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=37
#>        Working on 10 repeats of 2.08205926079142e+41 combinations
#>        Total time for n=10, k=37 & rep = 10: 0.654 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=38
#>        Working on 10 repeats of 9.47884873991883e+41 combinations
#>        Total time for n=10, k=38 & rep = 10: 0.674 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=39
#>        Working on 10 repeats of 4.18041534170788e+42 combinations
#>        Total time for n=10, k=39 & rep = 10: 0.688 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=40
#>        Working on 10 repeats of 1.78712755858008e+43 combinations
#>        Total time for n=10, k=40 & rep = 10: 0.717 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=41
#>        Working on 10 repeats of 7.41004109655167e+43 combinations
#>        Total time for n=10, k=41 & rep = 10: 0.732 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=42
#>        Working on 10 repeats of 2.98165939361247e+44 combinations
#>        Total time for n=10, k=42 & rep = 10: 0.751 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=43
#>        Working on 10 repeats of 1.16492739099279e+45 combinations
#>        Total time for n=10, k=43 & rep = 10: 0.763 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=44
#>        Working on 10 repeats of 4.42142896126796e+45 combinations
#>        Total time for n=10, k=44 & rep = 10: 0.783 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=45
#>        Working on 10 repeats of 1.6310160168233e+46 combinations
#>        Total time for n=10, k=45 & rep = 10: 0.793 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=46
#>        Working on 10 repeats of 5.85038353860539e+46 combinations
#>        Total time for n=10, k=46 & rep = 10: 0.817 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=47
#>        Working on 10 repeats of 2.04141042623676e+47 combinations
#>        Total time for n=10, k=47 & rep = 10: 0.834 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=48
#>        Working on 10 repeats of 6.93228957242899e+47 combinations
#>        Total time for n=10, k=48 & rep = 10: 0.847 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=49
#>        Working on 10 repeats of 2.29189981782347e+48 combinations
#>        Total time for n=10, k=49 & rep = 10: 0.872 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=50
#>        Working on 10 repeats of 7.37991741339163e+48 combinations
#>        Total time for n=10, k=50 & rep = 10: 0.898 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=51
#>        Working on 10 repeats of 2.31526820812283e+49 combinations
#>        Total time for n=10, k=51 & rep = 10: 0.905 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=52
#>        Working on 10 repeats of 7.0793777902217e+49 combinations
#>        Total time for n=10, k=52 & rep = 10: 0.923 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=53
#>        Working on 10 repeats of 2.11045602048119e+50 combinations
#>        Total time for n=10, k=53 & rep = 10: 0.945 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=54
#>        Working on 10 repeats of 6.13595546695463e+50 combinations
#>        Total time for n=10, k=54 & rep = 10: 0.955 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=55
#>        Working on 10 repeats of 1.7403800960817e+51 combinations
#>        Total time for n=10, k=55 & rep = 10: 0.972 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=56
#>        Working on 10 repeats of 4.81712348022608e+51 combinations
#>        Total time for n=10, k=56 & rep = 10: 0.988 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=57
#>        Working on 10 repeats of 1.30146844904354e+52 combinations
#>        Total time for n=10, k=57 & rep = 10: 1.005 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=58
#>        Working on 10 repeats of 3.43318401213212e+52 combinations
#>        Total time for n=10, k=58 & rep = 10: 1.034 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=59
#>        Working on 10 repeats of 8.8448130482046e+52 combinations
#>        Total time for n=10, k=59 & rep = 10: 1.049 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=60
#>        Working on 10 repeats of 2.22594461713154e+53 combinations
#>        Total time for n=10, k=60 & rep = 10: 1.062 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=61
#>        Working on 10 repeats of 5.47363430442175e+53 combinations
#>        Total time for n=10, k=61 & rep = 10: 1.077 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=62
#>        Working on 10 repeats of 1.31543792154652e+54 combinations
#>        Total time for n=10, k=62 & rep = 10: 1.095 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=63
#>        Working on 10 repeats of 3.09023511728388e+54 combinations
#>        Total time for n=10, k=63 & rep = 10: 1.103 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=64
#>        Working on 10 repeats of 7.09788378501144e+54 combinations
#>        Total time for n=10, k=64 & rep = 10: 1.13 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=65
#>        Working on 10 repeats of 1.59429389632563e+55 combinations
#>        Total time for n=10, k=65 & rep = 10: 1.147 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=66
#>        Working on 10 repeats of 3.50261537829118e+55 combinations
#>        Total time for n=10, k=66 & rep = 10: 1.167 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=67
#>        Working on 10 repeats of 7.52800917125258e+55 combinations
#>        Total time for n=10, k=67 & rep = 10: 1.176 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=68
#>        Working on 10 repeats of 1.58309604630754e+56 combinations
#>        Total time for n=10, k=68 & rep = 10: 1.193 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=69
#>        Working on 10 repeats of 3.25796577645906e+56 combinations
#>        Total time for n=10, k=69 & rep = 10: 1.258 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=70
#>        Working on 10 repeats of 6.56247392115308e+56 combinations
#>        Total time for n=10, k=70 & rep = 10: 1.255 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=71
#>        Working on 10 repeats of 1.29400894219921e+57 combinations
#>        Total time for n=10, k=71 & rep = 10: 1.222 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=72
#>        Working on 10 repeats of 2.49815615230128e+57 combinations
#>        Total time for n=10, k=72 & rep = 10: 1.212 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=73
#>        Working on 10 repeats of 4.72254176736412e+57 combinations
#>        Total time for n=10, k=73 & rep = 10: 1.255 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=74
#>        Working on 10 repeats of 8.74308408282277e+57 combinations
#>        Total time for n=10, k=74 & rep = 10: 1.247 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=75
#>        Working on 10 repeats of 1.58541258035182e+58 combinations
#>        Total time for n=10, k=75 & rep = 10: 1.25 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=76
#>        Working on 10 repeats of 2.81619339930916e+58 combinations
#>        Total time for n=10, k=76 & rep = 10: 1.269 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=77
#>        Working on 10 repeats of 4.90090799360302e+58 combinations
#>        Total time for n=10, k=77 & rep = 10: 1.288 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=78
#>        Working on 10 repeats of 8.35667645063076e+58 combinations
#>        Total time for n=10, k=78 & rep = 10: 1.306 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=79
#>        Working on 10 repeats of 1.39630543225731e+59 combinations
#>        Total time for n=10, k=79 & rep = 10: 1.314 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=80
#>        Working on 10 repeats of 2.28645014532135e+59 combinations
#>        Total time for n=10, k=80 & rep = 10: 1.337 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=81
#>        Working on 10 repeats of 3.66961134434288e+59 combinations
#>        Total time for n=10, k=81 & rep = 10: 1.35 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=82
#>        Working on 10 repeats of 5.77292516366116e+59 combinations
#>        Total time for n=10, k=82 & rep = 10: 1.355 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=83
#>        Working on 10 repeats of 8.90282434877879e+59 combinations
#>        Total time for n=10, k=83 & rep = 10: 1.408 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=84
#>        Working on 10 repeats of 1.34602225273205e+60 combinations
#>        Total time for n=10, k=84 & rep = 10: 1.406 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=85
#>        Working on 10 repeats of 1.99528004522631e+60 combinations
#>        Total time for n=10, k=85 & rep = 10: 1.414 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=86
#>        Working on 10 repeats of 2.90011634480568e+60 combinations
#>        Total time for n=10, k=86 & rep = 10: 1.425 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=87
#>        Working on 10 repeats of 4.13349915811388e+60 combinations
#>        Total time for n=10, k=87 & rep = 10: 1.447 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=88
#>        Working on 10 repeats of 5.77750450509094e+60 combinations
#>        Total time for n=10, k=88 & rep = 10: 1.465 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=89
#>        Working on 10 repeats of 7.91972527664159e+60 combinations
#>        Total time for n=10, k=89 & rep = 10: 1.481 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=90
#>        Working on 10 repeats of 1.06476306497071e+61 combinations
#>        Total time for n=10, k=90 & rep = 10: 1.495 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=91
#>        Working on 10 repeats of 1.40408316259872e+61 combinations
#>        Total time for n=10, k=91 & rep = 10: 1.512 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=92
#>        Working on 10 repeats of 1.81615104727443e+61 combinations
#>        Total time for n=10, k=92 & rep = 10: 1.606 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=93
#>        Working on 10 repeats of 2.3043636943912e+61 combinations
#>        Total time for n=10, k=93 & rep = 10: 1.824 minutes
#>        There were 2 of 2 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=94
#>        Working on 10 repeats of 2.86819736429546e+61 combinations
#>        Total time for n=10, k=94 & rep = 10: 1.565 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=95
#>        Working on 10 repeats of 3.50221993956075e+61 combinations
#>        Total time for n=10, k=95 & rep = 10: 1.568 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=96
#>        Working on 10 repeats of 4.1953676359321e+61 combinations
#>        Total time for n=10, k=96 & rep = 10: 1.66 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=97
#>        Working on 10 repeats of 4.93063825253891e+61 combinations
#>        Total time for n=10, k=97 & rep = 10: 1.679 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=98
#>        Working on 10 repeats of 5.68532778098881e+61 combinations
#>        Total time for n=10, k=98 & rep = 10: 1.84 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=99
#>        Working on 10 repeats of 6.43188597445179e+61 combinations
#>        Total time for n=10, k=99 & rep = 10: 1.876 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=100
#>        Working on 10 repeats of 7.13939343164141e+61 combinations
#>        Total time for n=10, k=100 & rep = 10: 1.967 minutes
#>        There were 0 of 1 sets identified by Fischers algorith as decisive (0%)
#> 
#> Working on n=10, k=101
#>        Working on 10 repeats of 7.77557700475797e+61 combinations
#>        Total time for n=10, k=101 & rep = 10: 1.953 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=102
#>        Working on 10 repeats of 8.30919503449644e+61 combinations
#>        Total time for n=10, k=102 & rep = 10: 1.98 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=103
#>        Working on 10 repeats of 8.71255401675361e+61 combinations
#>        Total time for n=10, k=103 & rep = 10: 1.873 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=10, k=104
#>        Working on 10 repeats of 8.96387769031386e+61 combinations
#>        Total time for n=10, k=104 & rep = 10: 1.975 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=105
#>        Working on 10 repeats of 9.04924795403077e+61 combinations
#>        Total time for n=10, k=105 & rep = 10: 2.155 minutes
#>        There were 2 of 2 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=106
#>        Working on 10 repeats of 8.96387769031386e+61 combinations
#>        Total time for n=10, k=106 & rep = 10: 2.067 minutes
#>        There were 2 of 2 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=107
#>        Working on 10 repeats of 8.71255401675361e+61 combinations
#>        Total time for n=10, k=107 & rep = 10: 2.081 minutes
#>        There were 3 of 3 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=108
#>        Working on 10 repeats of 8.30919503449644e+61 combinations
#>        Total time for n=10, k=108 & rep = 10: 2.038 minutes
#>        There were 3 of 3 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=109
#>        Working on 10 repeats of 7.77557700475797e+61 combinations
#>        Total time for n=10, k=109 & rep = 10: 2.108 minutes
#>        There were 4 of 4 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=110
#>        Working on 10 repeats of 7.13939343164141e+61 combinations
#>        Total time for n=10, k=110 & rep = 10: 2.145 minutes
#>        There were 5 of 5 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=111
#>        Working on 10 repeats of 6.43188597445179e+61 combinations
#>        Total time for n=10, k=111 & rep = 10: 2.266 minutes
#>        There were 6 of 6 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=112
#>        Working on 10 repeats of 5.68532778098881e+61 combinations
#>        Total time for n=10, k=112 & rep = 10: 2.194 minutes
#>        There were 5 of 5 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=113
#>        Working on 10 repeats of 4.93063825253891e+61 combinations
#>        Total time for n=10, k=113 & rep = 10: 2.331 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=114
#>        Working on 10 repeats of 4.1953676359321e+61 combinations
#>        Total time for n=10, k=114 & rep = 10: 2.425 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=115
#>        Working on 10 repeats of 3.50221993956075e+61 combinations
#>        Total time for n=10, k=115 & rep = 10: 2.459 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=116
#>        Working on 10 repeats of 2.86819736429546e+61 combinations
#>        Total time for n=10, k=116 & rep = 10: 2.088 minutes
#>        There were 5 of 5 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=117
#>        Working on 10 repeats of 2.3043636943912e+61 combinations
#>        Total time for n=10, k=117 & rep = 10: 2.141 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=118
#>        Working on 10 repeats of 1.81615104727443e+61 combinations
#>        Total time for n=10, k=118 & rep = 10: 2.185 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=119
#>        Working on 10 repeats of 1.40408316259872e+61 combinations
#>        Total time for n=10, k=119 & rep = 10: 2.142 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=120
#>        Working on 10 repeats of 1.06476306497071e+61 combinations
#>        Total time for n=10, k=120 & rep = 10: 2.21 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=121
#>        Working on 10 repeats of 7.91972527664159e+60 combinations
#>        Total time for n=10, k=121 & rep = 10: 2.223 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=122
#>        Working on 10 repeats of 5.77750450509094e+60 combinations
#>        Total time for n=10, k=122 & rep = 10: 2.24 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=123
#>        Working on 10 repeats of 4.13349915811388e+60 combinations
#>        Total time for n=10, k=123 & rep = 10: 2.2 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=124
#>        Working on 10 repeats of 2.90011634480568e+60 combinations
#>        Total time for n=10, k=124 & rep = 10: 2.223 minutes
#>        There were 6 of 6 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=125
#>        Working on 10 repeats of 1.99528004522631e+60 combinations
#>        Total time for n=10, k=125 & rep = 10: 2.371 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=126
#>        Working on 10 repeats of 1.34602225273205e+60 combinations
#>        Total time for n=10, k=126 & rep = 10: 2.284 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=127
#>        Working on 10 repeats of 8.90282434877879e+59 combinations
#>        Total time for n=10, k=127 & rep = 10: 2.308 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=128
#>        Working on 10 repeats of 5.77292516366116e+59 combinations
#>        Total time for n=10, k=128 & rep = 10: 2.327 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=129
#>        Working on 10 repeats of 3.66961134434288e+59 combinations
#>        Total time for n=10, k=129 & rep = 10: 2.316 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=130
#>        Working on 10 repeats of 2.28645014532135e+59 combinations
#>        Total time for n=10, k=130 & rep = 10: 2.433 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=131
#>        Working on 10 repeats of 1.39630543225731e+59 combinations
#>        Total time for n=10, k=131 & rep = 10: 2.425 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=132
#>        Working on 10 repeats of 8.35667645063076e+58 combinations
#>        Total time for n=10, k=132 & rep = 10: 2.416 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=133
#>        Working on 10 repeats of 4.90090799360302e+58 combinations
#>        Total time for n=10, k=133 & rep = 10: 2.483 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=134
#>        Working on 10 repeats of 2.81619339930916e+58 combinations
#>        Total time for n=10, k=134 & rep = 10: 2.469 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=135
#>        Working on 10 repeats of 1.58541258035182e+58 combinations
#>        Total time for n=10, k=135 & rep = 10: 2.478 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=136
#>        Working on 10 repeats of 8.74308408282277e+57 combinations
#>        Total time for n=10, k=136 & rep = 10: 2.529 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=137
#>        Working on 10 repeats of 4.72254176736412e+57 combinations
#>        Total time for n=10, k=137 & rep = 10: 2.531 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=138
#>        Working on 10 repeats of 2.49815615230128e+57 combinations
#>        Total time for n=10, k=138 & rep = 10: 2.515 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=139
#>        Working on 10 repeats of 1.29400894219921e+57 combinations
#>        Total time for n=10, k=139 & rep = 10: 2.558 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=140
#>        Working on 10 repeats of 6.56247392115308e+56 combinations
#>        Total time for n=10, k=140 & rep = 10: 2.553 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=141
#>        Working on 10 repeats of 3.25796577645906e+56 combinations
#>        Total time for n=10, k=141 & rep = 10: 2.569 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=142
#>        Working on 10 repeats of 1.58309604630754e+56 combinations
#>        Total time for n=10, k=142 & rep = 10: 2.59 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=143
#>        Working on 10 repeats of 7.52800917125258e+55 combinations
#>        Total time for n=10, k=143 & rep = 10: 2.616 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=144
#>        Working on 10 repeats of 3.50261537829118e+55 combinations
#>        Total time for n=10, k=144 & rep = 10: 2.666 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=145
#>        Working on 10 repeats of 1.59429389632563e+55 combinations
#>        Total time for n=10, k=145 & rep = 10: 2.608 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=146
#>        Working on 10 repeats of 7.09788378501144e+54 combinations
#>        Total time for n=10, k=146 & rep = 10: 2.625 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=147
#>        Working on 10 repeats of 3.09023511728388e+54 combinations
#>        Total time for n=10, k=147 & rep = 10: 2.63 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=148
#>        Working on 10 repeats of 1.31543792154652e+54 combinations
#>        Total time for n=10, k=148 & rep = 10: 2.645 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=149
#>        Working on 10 repeats of 5.47363430442175e+53 combinations
#>        Total time for n=10, k=149 & rep = 10: 2.663 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=150
#>        Working on 10 repeats of 2.22594461713154e+53 combinations
#>        Total time for n=10, k=150 & rep = 10: 2.687 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=151
#>        Working on 10 repeats of 8.8448130482046e+52 combinations
#>        Total time for n=10, k=151 & rep = 10: 2.696 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=152
#>        Working on 10 repeats of 3.43318401213212e+52 combinations
#>        Total time for n=10, k=152 & rep = 10: 2.738 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=153
#>        Working on 10 repeats of 1.30146844904354e+52 combinations
#>        Total time for n=10, k=153 & rep = 10: 2.725 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=154
#>        Working on 10 repeats of 4.81712348022608e+51 combinations
#>        Total time for n=10, k=154 & rep = 10: 2.733 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=155
#>        Working on 10 repeats of 1.7403800960817e+51 combinations
#>        Total time for n=10, k=155 & rep = 10: 2.755 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=156
#>        Working on 10 repeats of 6.13595546695463e+50 combinations
#>        Total time for n=10, k=156 & rep = 10: 2.724 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=157
#>        Working on 10 repeats of 2.11045602048119e+50 combinations
#>        Total time for n=10, k=157 & rep = 10: 2.661 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=158
#>        Working on 10 repeats of 7.0793777902217e+49 combinations
#>        Total time for n=10, k=158 & rep = 10: 2.728 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=159
#>        Working on 10 repeats of 2.31526820812283e+49 combinations
#>        Total time for n=10, k=159 & rep = 10: 2.683 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=160
#>        Working on 10 repeats of 7.37991741339163e+48 combinations
#>        Total time for n=10, k=160 & rep = 10: 2.693 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=161
#>        Working on 10 repeats of 2.29189981782347e+48 combinations
#>        Total time for n=10, k=161 & rep = 10: 2.733 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=162
#>        Working on 10 repeats of 6.93228957242899e+47 combinations
#>        Total time for n=10, k=162 & rep = 10: 2.721 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=163
#>        Working on 10 repeats of 2.04141042623676e+47 combinations
#>        Total time for n=10, k=163 & rep = 10: 2.753 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=164
#>        Working on 10 repeats of 5.85038353860539e+46 combinations
#>        Total time for n=10, k=164 & rep = 10: 2.751 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=165
#>        Working on 10 repeats of 1.6310160168233e+46 combinations
#>        Total time for n=10, k=165 & rep = 10: 2.771 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=166
#>        Working on 10 repeats of 4.42142896126796e+45 combinations
#>        Total time for n=10, k=166 & rep = 10: 2.769 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=167
#>        Working on 10 repeats of 1.16492739099279e+45 combinations
#>        Total time for n=10, k=167 & rep = 10: 2.781 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=168
#>        Working on 10 repeats of 2.98165939361247e+44 combinations
#>        Total time for n=10, k=168 & rep = 10: 2.824 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=169
#>        Working on 10 repeats of 7.41004109655167e+43 combinations
#>        Total time for n=10, k=169 & rep = 10: 2.823 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=170
#>        Working on 10 repeats of 1.78712755858008e+43 combinations
#>        Total time for n=10, k=170 & rep = 10: 2.832 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=171
#>        Working on 10 repeats of 4.18041534170788e+42 combinations
#>        Total time for n=10, k=171 & rep = 10: 2.855 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=172
#>        Working on 10 repeats of 9.47884873991883e+41 combinations
#>        Total time for n=10, k=172 & rep = 10: 508.407 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=173
#>        Working on 10 repeats of 2.08205926079142e+41 combinations
#>        Total time for n=10, k=173 & rep = 10: 2.913 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=174
#>        Working on 10 repeats of 4.427367393637e+40 combinations
#>        Total time for n=10, k=174 & rep = 10: 2.959 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=175
#>        Working on 10 repeats of 9.10772720976741e+39 combinations
#>        Total time for n=10, k=175 & rep = 10: 2.907 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=176
#>        Working on 10 repeats of 1.81119575194243e+39 combinations
#>        Total time for n=10, k=176 & rep = 10: 2.918 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=177
#>        Working on 10 repeats of 3.47913308282726e+38 combinations
#>        Total time for n=10, k=177 & rep = 10: 2.939 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=178
#>        Working on 10 repeats of 6.4500781872639e+37 combinations
#>        Total time for n=10, k=178 & rep = 10: 2.948 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=179
#>        Working on 10 repeats of 1.15308660330975e+37 combinations
#>        Total time for n=10, k=179 & rep = 10: 2.96 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=180
#>        Working on 10 repeats of 1.98587137236679e+36 combinations
#>        Total time for n=10, k=180 & rep = 10: 2.973 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=181
#>        Working on 10 repeats of 3.29149951221018e+35 combinations
#>        Total time for n=10, k=181 & rep = 10: 2.984 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=182
#>        Working on 10 repeats of 5.24469702495028e+34 combinations
#>        Total time for n=10, k=182 & rep = 10: 3.001 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=183
#>        Working on 10 repeats of 8.02467304363978e+33 combinations
#>        Total time for n=10, k=183 & rep = 10: 3.01 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=184
#>        Working on 10 repeats of 1.17753354444714e+33 combinations
#>        Total time for n=10, k=184 & rep = 10: 3.022 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=185
#>        Working on 10 repeats of 1.6549120084122e+32 combinations
#>        Total time for n=10, k=185 & rep = 10: 3.046 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=186
#>        Working on 10 repeats of 2.22434409732822e+31 combinations
#>        Total time for n=10, k=186 & rep = 10: 3.052 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=187
#>        Working on 10 repeats of 2.85477317304157e+30 combinations
#>        Total time for n=10, k=187 & rep = 10: 3.063 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=188
#>        Working on 10 repeats of 3.49254164787001e+29 combinations
#>        Total time for n=10, k=188 & rep = 10: 3.081 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=189
#>        Working on 10 repeats of 4.06539239434604e+28 combinations
#>        Total time for n=10, k=189 & rep = 10: 3.097 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=190
#>        Working on 10 repeats of 4.49332843585615e+27 combinations
#>        Total time for n=10, k=190 & rep = 10: 3.105 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=191
#>        Working on 10 repeats of 4.70505595377607e+26 combinations
#>        Total time for n=10, k=191 & rep = 10: 3.121 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=192
#>        Working on 10 repeats of 4.65604495425757e+25 combinations
#>        Total time for n=10, k=192 & rep = 10: 49.678 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=193
#>        Working on 10 repeats of 4.3424253459397e+24 combinations
#>        Total time for n=10, k=193 & rep = 10: 3.166 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=194
#>        Working on 10 repeats of 3.80521808664819e+23 combinations
#>        Total time for n=10, k=194 & rep = 10: 3.173 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=195
#>        Working on 10 repeats of 3.12223022494211e+22 combinations
#>        Total time for n=10, k=195 & rep = 10: 3.285 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=196
#>        Working on 10 repeats of 2.38946190684345e+21 combinations
#>        Total time for n=10, k=196 & rep = 10: 3.258 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=197
#>        Working on 10 repeats of 1.6980947561324e+20 combinations
#>        Total time for n=10, k=197 & rep = 10: 3.319 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=198
#>        Working on 10 repeats of 11149106984707682304 combinations
#>        Total time for n=10, k=198 & rep = 10: 3.508 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=199
#>        Working on 10 repeats of 672307958876845184 combinations
#>        Total time for n=10, k=199 & rep = 10: 3.444 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=200
#>        Working on 10 repeats of 36976937738226480 combinations
#>        Total time for n=10, k=200 & rep = 10: 3.503 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=201
#>        Working on 10 repeats of 1839648643692860 combinations
#>        Total time for n=10, k=201 & rep = 10: 3.681 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=202
#>        Working on 10 repeats of 81964543530870 combinations
#>        Total time for n=10, k=202 & rep = 10: 3.394 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=10, k=203
#>        Working on 10 repeats of 3230129794320 combinations
#>        Total time for n=10, k=203 & rep = 10: 3.538 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
SimulationResults_n10_test = rbindlist(dumTab)

SimulationResults_n10_test[,posRate := NR_FTT/NR_PhyloDec]
SimulationResults_n10_test[NR_PhyloDec==0,posRate := NA]
SimulationResults_n10_test[NR_PhyloDec>0,]
#>        k  time NR_NotPhyloDec NR_PhyloDec NR_FTT posRate
#>   1:  92 1.606              9           1      1       1
#>   2:  93 1.824              8           2      2       1
#>   3:  96 1.660              9           1      1       1
#>   4:  97 1.679              9           1      1       1
#>   5:  98 1.840              9           1      1       1
#>  ---                                                    
#> 103: 199 3.444              0          10     10       1
#> 104: 200 3.503              0          10     10       1
#> 105: 201 3.681              0          10     10       1
#> 106: 202 3.394              0          10     10       1
#> 107: 203 3.538              0          10     10       1
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
  
  outfn = paste0("../results/02_SimulationsData_n10/SimulationResults_n10_k",j,"_rep10000.RData")
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
#> TOTAL TIME : 907.801 minutes
```
