
<!-- 02_Simulation_n09.md is generated from 02_Simulation_n09.Rmd. Please edit that file -->

# Introduction

<!-- badges: start -->
<!-- badges: end -->

I want to run a simulation for n=9 taxa. There are 70 possible
quadruples with n=9 taxa, and I want to test 100000 possible
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

In *myTab_n9*, there are all possible Four-way-partitions (4WPP) given
n=9 taxa.

``` r
# Load initial input data
test1 = FTT_createInput(fn="../data/S9_Decisive.txt",sepSym = "_",c = 4)
#> Input contains 7 sets with 9 different taxa. 
#> The largest set has 8 taxa.

# Sanity check of FTT package
test_FTT = FTT_algorithmRed(data = test1$data,verbose = T, c=4, n=9)
#> Using 112 of 126 4-tuples as input for algorithm (9 unique taxa). 
#>  This leaves 14 4-tuples unsolved.
#> In round #1, 10 4-tuples could be resolved ...
#> In round #2, 4 4-tuples could be resolved ...
#> FIXING TAXON TRACEABLE
#>  It follows that the set is phylogenetically decisive
#test_NRC = FTT_findNRC(data = test1)

load("../results/01_partitions/partitions_n09.RData")
```

# Test-Loop with 10 combinations per k

To test less $k$, I use the more strict lower bound for touple numbers
(see FTT R Package for more explanation).

``` r
test1$data[,status:=NA]
n=9
LowerBound = (1/6) * 4 * choose(n,2)
UpperBound = choose(n,4) - n + 3

dumTab = foreach(j=c(LowerBound:UpperBound))%do%{
  # j=24
  message("\nWorking on n=9, k=",j)
  time1 = Sys.time()
  myTest = SimulationFunction(number_taxa = n, 
                              number_quads = j,
                              repeats = 10,
                              data1 = test1,
                              data2 = myTab_n09,
                              verbose = F,
                              FFT_only_if_PhyloDec = T)
  time2 = Sys.time()
  x0 = as.numeric(round(difftime(time2,time1,units = "mins"),3))
  message("       Total time for n=9, k=",j," & rep = 10: " ,round(difftime(time2,time1,units = "mins"),3)," minutes")
  
  outfn = paste0("../temp/02_SimulationsData_n09/SimulationResults_n9_k",j,".RData")
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
#> Working on n=9, k=24
#>        Working on 10 repeats of 3.97663318623775e+25 combinations
#>        Total time for n=9, k=24 & rep = 10: 0.065 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=25
#>        Working on 10 repeats of 1.622466339985e+26 combinations
#>        Total time for n=9, k=25 & rep = 10: 0.066 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=26
#>        Working on 10 repeats of 6.30265770532635e+26 combinations
#>        Total time for n=9, k=26 & rep = 10: 0.069 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=27
#>        Working on 10 repeats of 2.33431766863939e+27 combinations
#>        Total time for n=9, k=27 & rep = 10: 0.072 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=28
#>        Working on 10 repeats of 8.25348032840355e+27 combinations
#>        Total time for n=9, k=28 & rep = 10: 0.076 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=29
#>        Working on 10 repeats of 2.78910714546051e+28 combinations
#>        Total time for n=9, k=29 & rep = 10: 0.085 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=30
#>        Working on 10 repeats of 9.01811310365565e+28 combinations
#>        Total time for n=9, k=30 & rep = 10: 0.084 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=31
#>        Working on 10 repeats of 2.79270599339011e+29 combinations
#>        Total time for n=9, k=31 & rep = 10: 0.092 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=32
#>        Working on 10 repeats of 8.29084591787702e+29 combinations
#>        Total time for n=9, k=32 & rep = 10: 0.084 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=33
#>        Working on 10 repeats of 2.36163489781951e+30 combinations
#>        Total time for n=9, k=33 & rep = 10: 0.088 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=34
#>        Working on 10 repeats of 6.4597660440357e+30 combinations
#>        Total time for n=9, k=34 & rep = 10: 0.089 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=35
#>        Working on 10 repeats of 1.69799564586081e+31 combinations
#>        Total time for n=9, k=35 & rep = 10: 0.093 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=36
#>        Working on 10 repeats of 4.29215566037039e+31 combinations
#>        Total time for n=9, k=36 & rep = 10: 0.094 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=37
#>        Working on 10 repeats of 1.04403786333334e+32 combinations
#>        Total time for n=9, k=37 & rep = 10: 0.098 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=38
#>        Working on 10 repeats of 2.44524657464913e+32 combinations
#>        Total time for n=9, k=38 & rep = 10: 0.098 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=39
#>        Working on 10 repeats of 5.51747945049032e+32 combinations
#>        Total time for n=9, k=39 & rep = 10: 0.102 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=40
#>        Working on 10 repeats of 1.20005178048167e+33 combinations
#>        Total time for n=9, k=40 & rep = 10: 0.11 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=41
#>        Working on 10 repeats of 2.51718178344935e+33 combinations
#>        Total time for n=9, k=41 & rep = 10: 0.108 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=42
#>        Working on 10 repeats of 5.09429646650455e+33 combinations
#>        Total time for n=9, k=42 & rep = 10: 0.109 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=43
#>        Working on 10 repeats of 9.95164891131125e+33 combinations
#>        Total time for n=9, k=43 & rep = 10: 0.111 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=44
#>        Working on 10 repeats of 1.87724286281556e+34 combinations
#>        Total time for n=9, k=44 & rep = 10: 0.116 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=45
#>        Working on 10 repeats of 3.42075366113052e+34 combinations
#>        Total time for n=9, k=45 & rep = 10: 0.118 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=46
#>        Working on 10 repeats of 6.02350101199075e+34 combinations
#>        Total time for n=9, k=46 & rep = 10: 0.12 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=47
#>        Working on 10 repeats of 1.02527676799842e+35 combinations
#>        Total time for n=9, k=47 & rep = 10: 0.122 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=48
#>        Working on 10 repeats of 1.68743468066407e+35 combinations
#>        Total time for n=9, k=48 & rep = 10: 0.138 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=49
#>        Working on 10 repeats of 2.68612051207747e+35 combinations
#>        Total time for n=9, k=49 & rep = 10: 0.127 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=50
#>        Working on 10 repeats of 4.13662558859935e+35 combinations
#>        Total time for n=9, k=50 & rep = 10: 0.13 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=51
#>        Working on 10 repeats of 6.16438323006961e+35 combinations
#>        Total time for n=9, k=51 & rep = 10: 0.135 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=52
#>        Working on 10 repeats of 8.890937351062e+35 combinations
#>        Total time for n=9, k=52 & rep = 10: 0.138 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=53
#>        Working on 10 repeats of 1.24137615845016e+36 combinations
#>        Total time for n=9, k=53 & rep = 10: 0.138 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=54
#>        Working on 10 repeats of 1.67815665864561e+36 combinations
#>        Total time for n=9, k=54 & rep = 10: 0.141 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=55
#>        Working on 10 repeats of 2.19685962586328e+36 combinations
#>        Total time for n=9, k=55 & rep = 10: 0.143 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=56
#>        Working on 10 repeats of 2.78530416850525e+36 combinations
#>        Total time for n=9, k=56 & rep = 10: 0.146 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=57
#>        Working on 10 repeats of 3.42054897886612e+36 combinations
#>        Total time for n=9, k=57 & rep = 10: 0.149 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=58
#>        Working on 10 repeats of 4.06927378520275e+36 combinations
#>        Total time for n=9, k=58 & rep = 10: 0.189 minutes
#>        There were 1 of 1 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=59
#>        Working on 10 repeats of 4.69001046430148e+36 combinations
#>        Total time for n=9, k=59 & rep = 10: 0.158 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=60
#>        Working on 10 repeats of 5.23717835180331e+36 combinations
#>        Total time for n=9, k=60 & rep = 10: 0.159 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=61
#>        Working on 10 repeats of 5.66645526588558e+36 combinations
#>        Total time for n=9, k=61 & rep = 10: 0.158 minutes
#>        There were 0 of 0 sets identified by Fischers algorith as decisive (NaN%)
#> 
#> Working on n=9, k=62
#>        Working on 10 repeats of 5.94063858520265e+36 combinations
#>        Total time for n=9, k=62 & rep = 10: 0.208 minutes
#>        There were 0 of 1 sets identified by Fischers algorith as decisive (0%)
#> 
#> Working on n=9, k=63
#>        Working on 10 repeats of 6.03493443576138e+36 combinations
#>        Total time for n=9, k=63 & rep = 10: 0.227 minutes
#>        There were 1 of 2 sets identified by Fischers algorith as decisive (50%)
#> 
#> Working on n=9, k=64
#>        Working on 10 repeats of 5.94063858520265e+36 combinations
#>        Total time for n=9, k=64 & rep = 10: 0.201 minutes
#>        There were 0 of 1 sets identified by Fischers algorith as decisive (0%)
#> 
#> Working on n=9, k=65
#>        Working on 10 repeats of 5.66645526588558e+36 combinations
#>        Total time for n=9, k=65 & rep = 10: 0.316 minutes
#>        There were 5 of 6 sets identified by Fischers algorith as decisive (83.33%)
#> 
#> Working on n=9, k=66
#>        Working on 10 repeats of 5.23717835180331e+36 combinations
#>        Total time for n=9, k=66 & rep = 10: 0.267 minutes
#>        There were 4 of 4 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=67
#>        Working on 10 repeats of 4.69001046430148e+36 combinations
#>        Total time for n=9, k=67 & rep = 10: 0.273 minutes
#>        There were 4 of 4 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=68
#>        Working on 10 repeats of 4.06927378520275e+36 combinations
#>        Total time for n=9, k=68 & rep = 10: 0.339 minutes
#>        There were 6 of 7 sets identified by Fischers algorith as decisive (85.71%)
#> 
#> Working on n=9, k=69
#>        Working on 10 repeats of 3.42054897886612e+36 combinations
#>        Total time for n=9, k=69 & rep = 10: 0.301 minutes
#>        There were 4 of 5 sets identified by Fischers algorith as decisive (80%)
#> 
#> Working on n=9, k=70
#>        Working on 10 repeats of 2.78530416850525e+36 combinations
#>        Total time for n=9, k=70 & rep = 10: 0.264 minutes
#>        There were 4 of 4 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=71
#>        Working on 10 repeats of 2.19685962586328e+36 combinations
#>        Total time for n=9, k=71 & rep = 10: 0.303 minutes
#>        There were 6 of 6 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=72
#>        Working on 10 repeats of 1.67815665864561e+36 combinations
#>        Total time for n=9, k=72 & rep = 10: 0.294 minutes
#>        There were 6 of 6 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=73
#>        Working on 10 repeats of 1.24137615845016e+36 combinations
#>        Total time for n=9, k=73 & rep = 10: 0.296 minutes
#>        There were 6 of 6 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=74
#>        Working on 10 repeats of 8.890937351062e+35 combinations
#>        Total time for n=9, k=74 & rep = 10: 0.295 minutes
#>        There were 6 of 6 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=75
#>        Working on 10 repeats of 6.16438323006961e+35 combinations
#>        Total time for n=9, k=75 & rep = 10: 0.348 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=76
#>        Working on 10 repeats of 4.13662558859935e+35 combinations
#>        Total time for n=9, k=76 & rep = 10: 0.323 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=77
#>        Working on 10 repeats of 2.68612051207747e+35 combinations
#>        Total time for n=9, k=77 & rep = 10: 0.322 minutes
#>        There were 7 of 7 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=78
#>        Working on 10 repeats of 1.68743468066407e+35 combinations
#>        Total time for n=9, k=78 & rep = 10: 0.323 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=79
#>        Working on 10 repeats of 1.02527676799842e+35 combinations
#>        Total time for n=9, k=79 & rep = 10: 0.336 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=80
#>        Working on 10 repeats of 6.02350101199075e+34 combinations
#>        Total time for n=9, k=80 & rep = 10: 0.318 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=81
#>        Working on 10 repeats of 3.42075366113052e+34 combinations
#>        Total time for n=9, k=81 & rep = 10: 0.353 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=82
#>        Working on 10 repeats of 1.87724286281556e+34 combinations
#>        Total time for n=9, k=82 & rep = 10: 0.361 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=83
#>        Working on 10 repeats of 9.95164891131125e+33 combinations
#>        Total time for n=9, k=83 & rep = 10: 0.345 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=84
#>        Working on 10 repeats of 5.09429646650455e+33 combinations
#>        Total time for n=9, k=84 & rep = 10: 0.329 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=85
#>        Working on 10 repeats of 2.51718178344935e+33 combinations
#>        Total time for n=9, k=85 & rep = 10: 0.307 minutes
#>        There were 8 of 8 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=86
#>        Working on 10 repeats of 1.20005178048167e+33 combinations
#>        Total time for n=9, k=86 & rep = 10: 0.319 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=87
#>        Working on 10 repeats of 5.51747945049032e+32 combinations
#>        Total time for n=9, k=87 & rep = 10: 0.326 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=88
#>        Working on 10 repeats of 2.44524657464913e+32 combinations
#>        Total time for n=9, k=88 & rep = 10: 0.317 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=89
#>        Working on 10 repeats of 1.04403786333334e+32 combinations
#>        Total time for n=9, k=89 & rep = 10: 0.319 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=90
#>        Working on 10 repeats of 4.29215566037039e+31 combinations
#>        Total time for n=9, k=90 & rep = 10: 0.329 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=91
#>        Working on 10 repeats of 1.69799564586081e+31 combinations
#>        Total time for n=9, k=91 & rep = 10: 0.318 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=92
#>        Working on 10 repeats of 6.4597660440357e+30 combinations
#>        Total time for n=9, k=92 & rep = 10: 0.314 minutes
#>        There were 9 of 9 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=93
#>        Working on 10 repeats of 2.36163489781951e+30 combinations
#>        Total time for n=9, k=93 & rep = 10: 0.316 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=94
#>        Working on 10 repeats of 8.29084591787702e+29 combinations
#>        Total time for n=9, k=94 & rep = 10: 0.317 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=95
#>        Working on 10 repeats of 2.79270599339011e+29 combinations
#>        Total time for n=9, k=95 & rep = 10: 0.333 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=96
#>        Working on 10 repeats of 9.01811310365565e+28 combinations
#>        Total time for n=9, k=96 & rep = 10: 0.316 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=97
#>        Working on 10 repeats of 2.78910714546051e+28 combinations
#>        Total time for n=9, k=97 & rep = 10: 0.318 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=98
#>        Working on 10 repeats of 8.25348032840355e+27 combinations
#>        Total time for n=9, k=98 & rep = 10: 0.322 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=99
#>        Working on 10 repeats of 2.33431766863939e+27 combinations
#>        Total time for n=9, k=99 & rep = 10: 0.321 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=100
#>        Working on 10 repeats of 6.30265770532635e+26 combinations
#>        Total time for n=9, k=100 & rep = 10: 0.319 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=101
#>        Working on 10 repeats of 1.622466339985e+26 combinations
#>        Total time for n=9, k=101 & rep = 10: 0.319 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=102
#>        Working on 10 repeats of 3.97663318623775e+25 combinations
#>        Total time for n=9, k=102 & rep = 10: 0.315 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=103
#>        Working on 10 repeats of 9.26594140482581e+24 combinations
#>        Total time for n=9, k=103 & rep = 10: 0.314 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=104
#>        Working on 10 repeats of 2.0491985799134e+24 combinations
#>        Total time for n=9, k=104 & rep = 10: 0.323 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=105
#>        Working on 10 repeats of 4.29355892934236e+23 combinations
#>        Total time for n=9, k=105 & rep = 10: 0.322 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=106
#>        Working on 10 repeats of 8.50610731284808e+22 combinations
#>        Total time for n=9, k=106 & rep = 10: 0.338 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=107
#>        Working on 10 repeats of 1.58992660053235e+22 combinations
#>        Total time for n=9, k=107 & rep = 10: 0.32 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=108
#>        Working on 10 repeats of 2.79709309352914e+21 combinations
#>        Total time for n=9, k=108 & rep = 10: 0.354 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=109
#>        Working on 10 repeats of 4.61905281500225e+20 combinations
#>        Total time for n=9, k=109 & rep = 10: 0.367 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=110
#>        Working on 10 repeats of 71385361686398353418 combinations
#>        Total time for n=9, k=110 & rep = 10: 0.335 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=111
#>        Working on 10 repeats of 10289781864706068480 combinations
#>        Total time for n=9, k=111 & rep = 10: 0.347 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=112
#>        Working on 10 repeats of 1378095785451705600 combinations
#>        Total time for n=9, k=112 & rep = 10: 0.337 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=113
#>        Working on 10 repeats of 170737530940919296 combinations
#>        Total time for n=9, k=113 & rep = 10: 0.359 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=114
#>        Working on 10 repeats of 19470069317824128 combinations
#>        Total time for n=9, k=114 & rep = 10: 0.346 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=115
#>        Working on 10 repeats of 2031659407077300 combinations
#>        Total time for n=9, k=115 & rep = 10: 0.371 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=116
#>        Working on 10 repeats of 192657357567675 combinations
#>        Total time for n=9, k=116 & rep = 10: 0.321 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=117
#>        Working on 10 repeats of 16466440817750 combinations
#>        Total time for n=9, k=117 & rep = 10: 0.325 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=118
#>        Working on 10 repeats of 1255914977625 combinations
#>        Total time for n=9, k=118 & rep = 10: 0.322 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=119
#>        Working on 10 repeats of 84431259000 combinations
#>        Total time for n=9, k=119 & rep = 10: 0.325 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
#> 
#> Working on n=9, k=120
#>        Working on 10 repeats of 4925156775 combinations
#>        Total time for n=9, k=120 & rep = 10: 0.328 minutes
#>        There were 10 of 10 sets identified by Fischers algorith as decisive (100%)
SimulationResults_n09_test = rbindlist(dumTab)

SimulationResults_n09_test[,posRate := NR_FTT/NR_PhyloDec]
SimulationResults_n09_test[NR_PhyloDec==0,posRate := NA]
SimulationResults_n09_test[NR_PhyloDec>0,]
#>       k  time NR_NotPhyloDec NR_PhyloDec NR_FTT   posRate
#>  1:  58 0.189              9           1      1 1.0000000
#>  2:  62 0.208              9           1      0 0.0000000
#>  3:  63 0.227              8           2      1 0.5000000
#>  4:  64 0.201              9           1      0 0.0000000
#>  5:  65 0.316              4           6      5 0.8333333
#>  6:  66 0.267              6           4      4 1.0000000
#>  7:  67 0.273              6           4      4 1.0000000
#>  8:  68 0.339              3           7      6 0.8571429
#>  9:  69 0.301              5           5      4 0.8000000
#> 10:  70 0.264              6           4      4 1.0000000
#> 11:  71 0.303              4           6      6 1.0000000
#> 12:  72 0.294              4           6      6 1.0000000
#> 13:  73 0.296              4           6      6 1.0000000
#> 14:  74 0.295              4           6      6 1.0000000
#> 15:  75 0.348              2           8      8 1.0000000
#> 16:  76 0.323              2           8      8 1.0000000
#> 17:  77 0.322              3           7      7 1.0000000
#> 18:  78 0.323              2           8      8 1.0000000
#> 19:  79 0.336              1           9      9 1.0000000
#> 20:  80 0.318              2           8      8 1.0000000
#> 21:  81 0.353              0          10     10 1.0000000
#> 22:  82 0.361              0          10     10 1.0000000
#> 23:  83 0.345              0          10     10 1.0000000
#> 24:  84 0.329              0          10     10 1.0000000
#> 25:  85 0.307              2           8      8 1.0000000
#> 26:  86 0.319              1           9      9 1.0000000
#> 27:  87 0.326              0          10     10 1.0000000
#> 28:  88 0.317              0          10     10 1.0000000
#> 29:  89 0.319              0          10     10 1.0000000
#> 30:  90 0.329              0          10     10 1.0000000
#> 31:  91 0.318              0          10     10 1.0000000
#> 32:  92 0.314              1           9      9 1.0000000
#> 33:  93 0.316              0          10     10 1.0000000
#> 34:  94 0.317              0          10     10 1.0000000
#> 35:  95 0.333              0          10     10 1.0000000
#> 36:  96 0.316              0          10     10 1.0000000
#> 37:  97 0.318              0          10     10 1.0000000
#> 38:  98 0.322              0          10     10 1.0000000
#> 39:  99 0.321              0          10     10 1.0000000
#> 40: 100 0.319              0          10     10 1.0000000
#> 41: 101 0.319              0          10     10 1.0000000
#> 42: 102 0.315              0          10     10 1.0000000
#> 43: 103 0.314              0          10     10 1.0000000
#> 44: 104 0.323              0          10     10 1.0000000
#> 45: 105 0.322              0          10     10 1.0000000
#> 46: 106 0.338              0          10     10 1.0000000
#> 47: 107 0.320              0          10     10 1.0000000
#> 48: 108 0.354              0          10     10 1.0000000
#> 49: 109 0.367              0          10     10 1.0000000
#> 50: 110 0.335              0          10     10 1.0000000
#> 51: 111 0.347              0          10     10 1.0000000
#> 52: 112 0.337              0          10     10 1.0000000
#> 53: 113 0.359              0          10     10 1.0000000
#> 54: 114 0.346              0          10     10 1.0000000
#> 55: 115 0.371              0          10     10 1.0000000
#> 56: 116 0.321              0          10     10 1.0000000
#> 57: 117 0.325              0          10     10 1.0000000
#> 58: 118 0.322              0          10     10 1.0000000
#> 59: 119 0.325              0          10     10 1.0000000
#> 60: 120 0.328              0          10     10 1.0000000
#>       k  time NR_NotPhyloDec NR_PhyloDec NR_FTT   posRate
```

# Loop with all combinations

``` r
dumTab = foreach(j=c(LowerBound:UpperBound))%do%{
  # j=5
  message("\nWorking on n=9, k=",j)
  time1 = Sys.time()
  myTest = SimulationFunction(number_taxa = n, 
                              number_quads = j,
                              repeats = 10000,
                              data1 = test1,
                              data2 = myTab_n09,
                              verbose = T,
                              FFT_only_if_PhyloDec = T)
  time2 = Sys.time()
  x0 = as.numeric(round(difftime(time2,time1,units = "mins"),3))
  message("       Total time for n=9, k=",j," & rep = 100: " ,round(difftime(time2,time1,units = "mins"),3)," minutes")
  
  outfn = paste0("../results/02_SimulationsData_n09/SimulationResults_n9_k",k,"_rep10000.RData")
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
SimulationResults_n09 = rbindlist(dumTab)

SimulationResults_n09[,posRate := NR_FTT/NR_PhyloDec]
SimulationResults_n09[NR_PhyloDec==0,posRate := NA]
save(SimulationResults_n09, file = "../results/02_SimulationResults_n09.RData")

SimulationResults_n09[NR_PhyloDec>0,]
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
#> TOTAL TIME : 23.115 minutes
```
