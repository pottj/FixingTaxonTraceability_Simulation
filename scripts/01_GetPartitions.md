
<!-- 01_GetPartitions.md is generated from 01_GetPartitions.Rmd. Please edit that file -->

# Introduction

<!-- badges: start -->
<!-- badges: end -->

I want to create data tables containing all possible four-way partitions
with $n=6, ..., 10$ taxa. The data table should contain the following
columns:

-   set_i: taxa included in partition i of 4, separated by “\|”
-   allQuads: all quadruples possible with the given partition,
    separated by “\|”

According to Stirling numbers of the second kind, there are 34,105
different partitions for 10 taxa.

There are nine types of partitions:

1.  a \| b \| c \| defghij –> 6, 7, 8, 9, 10
2.  a \| b \| cd \| efghij –> 6, 7, 8, 9, 10
3.  a \| b \| cde \| fghij –> 8, 9, 10
4.  a \| b \| cdef \| ghij –> 10
5.  a \| bc \| de \| fghij –> 7, 8, 9, 10
6.  a \| bc \| def \| ghij –> 9, 10
7.  a \| bcd \| efg \| hij –> 10
8.  ab \| cd \| ef \| ghij –> 8, 9, 10
9.  ab \| cd \| efg \| hij –> 10

I create them separately and combine them at the end.

In addition, I test the table with some example data sets (not decisive:
Mareikes Mathematica example, decisive: some trees added to the previous
example)

# Initialize

``` r
rm(list = ls())
time0<-Sys.time()

source("../helperFunctions/StirlingFunction.R")
source("../helperFunctions/PartitioningFunction.R")
source("../helperFunctions/PartitioningFunction_HelpR.R")
source("../SourceFile.R")
```

# Check & general parameters

``` r
StirlingFunction(n=6,k=4)
#> [1] 65
StirlingFunction(n=7,k=4)
#> [1] 350
StirlingFunction(n=8,k=4)
#> [1] 1701
StirlingFunction(n=9,k=4)
#> [1] 7770
StirlingFunction(n=10,k=4)
#> [1] 34105
```

# Set 1 (all n)

a \| b \| c \| def == 1 + 1 + 1 + 3 a \| b \| c \| defg == 1 + 1 + 1 + 4
a \| b \| c \| defgh == 1 + 1 + 1 + 5 a \| b \| c \| defghi == 1 + 1 + 1
+ 6 a \| b \| c \| defghij == 1 + 1 + 1 + 7

``` r
tab1_n6 = PartitioningFunction(1,1,1,3)
tab1_n7 = PartitioningFunction(1,1,1,4)
tab1_n8 = PartitioningFunction(1,1,1,5)
tab1_n9 = PartitioningFunction(1,1,1,6)
tab1_n10 = PartitioningFunction(1,1,1,7)
```

# Set 2 (all n)

a \| b \| cd \| ef == 1 + 1 + 2 + 2 a \| b \| cd \| efg == 1 + 1 + 2 + 3
a \| b \| cd \| efgh == 1 + 1 + 2 + 4 a \| b \| cd \| efghi == 1 + 1 + 2
+ 5 a \| b \| cd \| efghij == 1 + 1 + 2 + 6

``` r
tab2_n6 = PartitioningFunction(1,1,2,2)
tab2_n7 = PartitioningFunction(1,1,2,3)
tab2_n8 = PartitioningFunction(1,1,2,4)
tab2_n9 = PartitioningFunction(1,1,2,5)
tab2_n10 = PartitioningFunction(1,1,2,6)
```

# Set 3 (8, 9, 10)

a \| b \| cde \| fgh == 1 + 1 + 3 + 3 a \| b \| cde \| fghi == 1 + 1 + 3
+ 4 a \| b \| cde \| fghij == 1 + 1 + 3 + 5

``` r
tab3_n8 = PartitioningFunction(1,1,3,3)
tab3_n9 = PartitioningFunction(1,1,3,4)
tab3_n10 = PartitioningFunction(1,1,3,5)
```

# Set 4 (10)

a \| b \| cdef \| ghij == 1 + 1 + 4 + 4

``` r
tab4_n10 = PartitioningFunction(1,1,4,4)
```

# Set 5 (7, 8, 9, 10)

a \| bc \| de \| fg == 1 + 2 + 2 + 2 a \| bc \| de \| fgh == 1 + 2 + 2 +
3 a \| bc \| de \| fghi == 1 + 2 + 2 + 4 a \| bc \| de \| fghij == 1 + 2
+ 2 + 5

``` r
tab5_n7 = PartitioningFunction(1,2,2,2)
tab5_n8 = PartitioningFunction(1,2,2,3)
tab5_n9 = PartitioningFunction(1,2,2,4)
tab5_n10 = PartitioningFunction(1,2,2,5)
```

# Set 6 (9, 10)

a \| bc \| def \| ghi == 1 + 2 + 3 + 3 a \| bc \| def \| ghij == 1 + 2 +
3 + 4

``` r
tab6_n9 = PartitioningFunction(1,2,3,3)
tab6_n10 = PartitioningFunction(1,2,3,4)
```

# Set 7 (10)

a \| bcd \| efg \| hij == 1 + 3 + 3 + 3

``` r
tab7_n10 = PartitioningFunction(1,3,3,3)
```

# Set 8 (8,9,10)

ab \| cd \| ef \| gh == 2 + 2 + 2 + 2 ab \| cd \| ef \| ghi == 2 + 2 + 2
+ 3 ab \| cd \| ef \| ghij == 2 + 2 + 2 + 4

``` r
tab8_n8 = PartitioningFunction(2,2,2,2)
tab8_n9 = PartitioningFunction(2,2,2,3)
tab8_n10 = PartitioningFunction(2,2,2,4)
```

# Set 9

ab \| bc \| efg \| hij == 2 + 2 + 3 + 3

``` r
tab9_n10 = PartitioningFunction(2,2,3,3)
```

# Combine & Check

``` r
myTab_n06 = rbind(tab1_n6,tab2_n6)
stopifnot(dim(myTab_n06)[1]==StirlingFunction(n=6,k=4))

myTab_n07 = rbind(tab1_n7,tab2_n7,tab5_n7)
stopifnot(dim(myTab_n07)[1]==StirlingFunction(n=7,k=4))

myTab_n08 = rbind(tab1_n8,tab2_n8,tab3_n8,tab5_n8,tab8_n8)
stopifnot(dim(myTab_n08)[1]==StirlingFunction(n=8,k=4))

myTab_n09 = rbind(tab1_n9,tab2_n9,tab3_n9,tab5_n9,tab6_n9,tab8_n9)
stopifnot(dim(myTab_n09)[1]==StirlingFunction(n=9,k=4))

myTab_n10 = rbind(tab1_n10,tab2_n10,tab3_n10,tab4_n10,
                  tab5_n10,tab6_n10,tab7_n10,tab8_n10,tab9_n10)
stopifnot(dim(myTab_n10)[1]==StirlingFunction(n=10,k=4))
```

# Save

``` r
save(myTab_n06,file = "../results/01_partitions/partitions_n06.RData")
save(myTab_n07,file = "../results/01_partitions/partitions_n07.RData")
save(myTab_n08,file = "../results/01_partitions/partitions_n08.RData")
save(myTab_n09,file = "../results/01_partitions/partitions_n09.RData")
save(myTab_n10,file = "../results/01_partitions/partitions_n10.RData")
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
#>  [5] rlang_1.0.6      cli_3.6.0        rmarkdown_2.20   iterators_1.0.14
#>  [9] tools_4.2.2      xfun_0.37        yaml_2.3.7       fastmap_1.1.1   
#> [13] compiler_4.2.2   htmltools_0.5.4  knitr_1.42
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
#> 
#> TOTAL TIME : 5.82 minutes
```
