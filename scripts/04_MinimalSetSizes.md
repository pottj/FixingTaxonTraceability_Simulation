Minimal decisive sets
================
Janne Pott
27/04/2023

<!-- 04_MinimalSetSizes.md is generated from 04_MinimalSetSizes.Rmd. Please edit that file -->

# Introduction

<!-- badges: start -->
<!-- badges: end -->

In our publication, we give lower bounds for each simulated $k_n$ for
which phylogenetic decisive sets are known.

Here, I load example sets with the respective $k_n$ and test for both
decisiveness using the *findNRC* algorithm of Parvini et al. and our
*fixingTaxonTraceR* algorithm.

# Initialize

I use a file names *SourceFile.R* that contains all relevant R packages
and user-/server-specific path to the R library. If using this code, you
must make all the necessary changes within the template source file.

``` r
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
```

# Load data

The example data is stored within the simulation repository. For $n=6$,
$n=7$ and $n=8$, the examples are explicitly also given in the
publication (Examples 3.10 and 3.23). The set for $n=9$ was taken from
the simulation with $k_9=46$, were exactly one set was decisive. We
tested if one additional input sample could be removed and indeed,
$\{2, 3, 5, 9\}$ was not necessary for covering all 4-way partitions.
For $n=10$, we used the Steiner Quadruple System as baseline and added
step-wise cross-quadruples which covered the most of the un-covered
partitions.

``` r
# load example sets
data_n06 = FTT_createInput(fn = "../data/S6_Decisive_minimal.txt",
                           sepSym = "_")
#> Input contains 9 sets with 6 different taxa. 
#> The largest set has 4 taxa.
data_n07 = FTT_createInput(fn = "../data/S7_Decisive_minimal.txt",
                           sepSym = "_")
#> Input contains 14 sets with 7 different taxa. 
#> The largest set has 4 taxa.
data_n08 = FTT_createInput(fn = "../data/S8_Decisive_minimal.txt",
                           sepSym = "_")
#> Input contains 26 sets with 8 different taxa. 
#> The largest set has 4 taxa.
data_n09 = FTT_createInput(fn = "../data/S9_Decisive_minimal.txt",
                           sepSym = "_")
#> Input contains 45 sets with 9 different taxa. 
#> The largest set has 4 taxa.
data_n10 = FTT_createInput(fn = "../data/S10_Decisive_minimal.txt",
                           sepSym = "_")
#> Input contains 58 sets with 10 different taxa. 
#> The largest set has 4 taxa.

# load partitions
load("../results/01_partitions/partitions_n06.RData")
load("../results/01_partitions/partitions_n07.RData")
load("../results/01_partitions/partitions_n08.RData")
load("../results/01_partitions/partitions_n09.RData")
load("../results/01_partitions/partitions_n10.RData")
```

# Test partitions

To test for phylogenetic decisiveness, one can simply test if all
possible 4-way partitions are covered. All possible partitions for
$n=6, ..., 10$ were created in the simulation, and are used to test the
5 minimal sets.

``` r
time3 = Sys.time()

dummy06 = data_n06$data[status=="input",]
dummy07 = data_n07$data[status=="input",]
dummy08 = data_n08$data[status=="input",]
dummy09 = data_n09$data[status=="input",]
dummy10 = data_n10$data[status=="input",]

myTab_n06[,count :=0]
myTab_n07[,count :=0]
myTab_n08[,count :=0]
myTab_n09[,count :=0]
myTab_n10[,count :=0]

for(i in 1:dim(dummy06)[1]){
  #i=1
  myInput = dummy06[i,ctuple]
  filt = grepl(myInput,myTab_n06$allQuads)
  myTab_n06[filt==T,status := "covered"]
  myTab_n06[filt==T,count := count + 1]
  myTab_n06
}
for(i in 1:dim(dummy07)[1]){
  #i=1
  myInput = dummy07[i,ctuple]
  filt = grepl(myInput,myTab_n07$allQuads)
  myTab_n07[filt==T,status := "covered"]
  myTab_n07[filt==T,count := count + 1]
  myTab_n07
}
for(i in 1:dim(dummy08)[1]){
  #i=1
  myInput = dummy08[i,ctuple]
  filt = grepl(myInput,myTab_n08$allQuads)
  myTab_n08[filt==T,status := "covered"]
  myTab_n08[filt==T,count := count + 1]
  myTab_n08
}
for(i in 1:dim(dummy09)[1]){
  #i=1
  myInput = dummy09[i,ctuple]
  filt = grepl(myInput,myTab_n09$allQuads)
  myTab_n09[filt==T,status := "covered"]
  myTab_n09[filt==T,count := count + 1]
  myTab_n09
}
for(i in 1:dim(dummy10)[1]){
  #i=1
  myInput = dummy10[i,ctuple]
  filt = grepl(myInput,myTab_n10$allQuads)
  myTab_n10[filt==T,status := "covered"]
  myTab_n10[filt==T,count := count + 1]
  myTab_n10
}

myTab_n06[,table(count)]
#> count
#>  1  2  3  4 
#> 12 30 20  3
myTab_n07[,table(count)]
#> count
#>   1   2   3   4   5 
#>  56 112 126  42  14
myTab_n08[,table(count)]
#> count
#>   1   2   3   4   5   6   8  10 
#>  44 268 264 693 156 247  20   9
myTab_n09[,table(count)]
#> count
#>    1    2    3    4    5    6    7    8    9   10   11   12   13   14   16 
#>   76  312  701 1114 1279 1309 1123  825  505  313  148   44   16    4    1
myTab_n10[,table(count)]
#> count
#>    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#>  144  696 1677 3316 4347 5380 5261 4404 3325 2306 1508  955  491  208   60   24 
#>   17 
#>    3

message("\nTOTAL TIME for partitions: " ,
        round(difftime(Sys.time(),time3,units = "mins"),3),
        " minutes")
#> 
#> TOTAL TIME for partitions: 0.11 minutes
```

For all 5 sets, the minimal count per partition is $1>0$, e.g. all
partitions are covered at least once, which gives us phylogenetic
decisiveness.

# Test FTT

All set sizes are below the lower bound of the FTT algorithm. Hence this
section is just a sanity check, as none of them can be resolved with
fixing taxa.

``` r
time1 = Sys.time()

test1_06_FTT = FTT_algorithmRed(data = data_n06$data,n=6,verbose = T)
#> Using 9 of 15 4-tuples as input for algorithm (6 unique taxa). 
#>  This leaves 6 4-tuples unsolved.
#> In round #1, 0 4-tuples could be resolved ...
#> NOT RESOLVABLE VIA THIS ALGORITHM! 
#>  There are 6 remaining cross 4-tuples.
test1_07_FTT = FTT_algorithmRed(data = data_n07$data,n=7,verbose = T)
#> Using 14 of 35 4-tuples as input for algorithm (7 unique taxa). 
#>  This leaves 21 4-tuples unsolved.
#> In round #1, 0 4-tuples could be resolved ...
#> NOT RESOLVABLE VIA THIS ALGORITHM! 
#>  There are 21 remaining cross 4-tuples.
test1_08_FTT = FTT_algorithmRed(data = data_n08$data,n=8,verbose = T)
#> Using 26 of 70 4-tuples as input for algorithm (8 unique taxa). 
#>  This leaves 44 4-tuples unsolved.
#> In round #1, 0 4-tuples could be resolved ...
#> NOT RESOLVABLE VIA THIS ALGORITHM! 
#>  There are 44 remaining cross 4-tuples.
test1_09_FTT = FTT_algorithmRed(data = data_n09$data,n=9,verbose = T)
#> Using 45 of 126 4-tuples as input for algorithm (9 unique taxa). 
#>  This leaves 81 4-tuples unsolved.
#> In round #1, 2 4-tuples could be resolved ...
#> In round #2, 0 4-tuples could be resolved ...
#> NOT RESOLVABLE VIA THIS ALGORITHM! 
#>  There are 79 remaining cross 4-tuples.
test1_10_FTT = FTT_algorithmRed(data = data_n10$data,n=10,verbose = T)
#> Using 58 of 210 4-tuples as input for algorithm (10 unique taxa). 
#>  This leaves 152 4-tuples unsolved.
#> In round #1, 0 4-tuples could be resolved ...
#> NOT RESOLVABLE VIA THIS ALGORITHM! 
#>  There are 152 remaining cross 4-tuples.

message("\nTOTAL TIME for FTT: " ,
        round(difftime(Sys.time(),time1,units = "mins"),3),
        " minutes")
#> 
#> TOTAL TIME for FTT: 0.074 minutes
```

# Test NRC

As all sets are decisive, there should only be rainbow-coloring
possible, e.g. regardless how you try to color the $n$ taxa with 4
colors, there is at least one quadruple with 4 colors. This is a direct
result of the 4-way partition property: each partition is one possible
coloring, and all partitions are covered at least once.

Please note: the *FTT_findNRC* algorithm is an R adaption of Parvini et
al. *findNRC* algorithm implemented in Python. The run time is not
optimized.

``` r
time2 = Sys.time()

test2_06_NRC = FTT_findNRC(data = data_n06)
#> THERE IS ONLY RAINBOW COLORING POSSIBLE OR 3-WAY COLORING
#>  It follows that the set is phylogenetically decisive
test2_07_NRC = FTT_findNRC(data = data_n07)
#> THERE IS ONLY RAINBOW COLORING POSSIBLE OR 3-WAY COLORING
#>  It follows that the set is phylogenetically decisive
test2_08_NRC = FTT_findNRC(data = data_n08)
#> THERE IS ONLY RAINBOW COLORING POSSIBLE OR 3-WAY COLORING
#>  It follows that the set is phylogenetically decisive
test2_09_NRC = FTT_findNRC(data = data_n09)
#> THERE IS ONLY RAINBOW COLORING POSSIBLE OR 3-WAY COLORING
#>  It follows that the set is phylogenetically decisive
test2_10_NRC = FTT_findNRC(data = data_n10)
#> THERE IS ONLY RAINBOW COLORING POSSIBLE OR 3-WAY COLORING
#>  It follows that the set is phylogenetically decisive

message("\nTOTAL TIME for NRC: " ,
        round(difftime(Sys.time(),time2,units = "mins"),3),
        " minutes")
#> 
#> TOTAL TIME for NRC: 36.603 minutes
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
message("\nTOTAL TIME : " ,
        round(difftime(Sys.time(),time0,units = "mins"),3),
        " minutes")
#> 
#> TOTAL TIME : 36.904 minutes
```
