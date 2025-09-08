# Tree_shrew_scRNA_code
Data Analysis Process of "Tree shrew immune cell atlas identifies NR1H3⁺ tissue macrophages with conserved anti-inflammatory function"

Raw sequencing data is available at SRA databank under the accession number PRJNA784168.

T. belangeri Cell Atlas website https://stmm-lab.shinyapps.io/tbca/  offering interactive search feature based on tissue, cell type and gene through ShinyCell tool(PMID: 33774659).

My environment: (not all are neccessary)

```txt
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 18.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Asia/Shanghai
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] patchwork_1.2.0     tibble_3.2.1        purrr_1.0.2         metap_1.4           multtest_2.60.0     Biobase_2.64.0     
 [7] BiocGenerics_0.50.0 dplyr_1.1.4         Matrix_1.7-0        cowplot_1.1.3       ggplot2_3.5.1       Seurat_5.1.0       
[13] SeuratObject_5.0.2  sp_2.1-4           

loaded via a namespace (and not attached):
  [1] mathjaxr_1.4-0         RColorBrewer_1.1-3     rstudioapi_0.16.0      jsonlite_1.8.8         magrittr_2.0.3        
  [6] TH.data_1.0-10         spatstat.utils_3.0-5   vctrs_0.6.5            ROCR_1.0-11            spatstat.explore_3.3-1
 [11] htmltools_0.5.8.1      plotrix_3.8-1          sctransform_0.4.1      parallelly_1.38.0      KernSmooth_2.23-24    
 [16] htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9             sandwich_3.0-1         plotly_4.10.4         
 [21] zoo_1.8-12             igraph_2.0.3           mime_0.12              lifecycle_1.0.4        pkgconfig_2.0.3       
 [26] R6_2.5.1               fastmap_1.2.0          rbibutils_2.2          fitdistrplus_1.2-1     future_1.33.2         
 [31] shiny_1.8.1.1          digest_0.6.36          numDeriv_2016.8-1.1    colorspace_2.1-1       tensor_1.5            
 [36] RSpectra_0.16-2        irlba_2.3.5.1          progressr_0.14.0       fansi_1.0.6            spatstat.sparse_3.1-0 
 [41] httr_1.4.7             TFisher_0.2.0          polyclip_1.10-7        abind_1.4-5            compiler_4.4.1        
 [46] withr_3.0.0            mutoss_0.1-12          fastDummies_1.7.3      MASS_7.3-61            tools_4.4.1           
 [51] lmtest_0.9-40          httpuv_1.6.15          future.apply_1.11.2    goftest_1.2-3          glue_1.7.0            
 [56] nlme_3.1-165           promises_1.3.0         grid_4.4.1             Rtsne_0.17             cluster_2.1.6         
 [61] reshape2_1.4.4         generics_0.1.3         gtable_0.3.5           spatstat.data_3.1-2    tidyr_1.3.1           
 [66] sn_2.0.0               data.table_1.15.4      utf8_1.2.4             spatstat.geom_3.3-2    RcppAnnoy_0.0.22      
 [71] ggrepel_0.9.5          RANN_2.6.1             pillar_1.9.0           stringr_1.5.1          spam_2.10-0           
 [76] RcppHNSW_0.6.0         later_1.3.2            splines_4.4.1          lattice_0.22-6         survival_3.7-0        
 [81] deldir_2.0-4           tidyselect_1.2.1       miniUI_0.1.1.1         pbapply_1.7-2          gridExtra_2.3         
 [86] scattermore_1.2        stats4_4.4.1           matrixStats_1.3.0      stringi_1.8.4          lazyeval_0.2.2        
 [91] codetools_0.2-20       cli_3.6.3              uwot_0.2.2             xtable_1.8-4           reticulate_1.38.0     
 [96] Rdpack_2.1.2           munsell_0.5.1          Rcpp_1.0.13-1          globals_0.16.3         spatstat.random_3.3-1 
[101] png_0.1-8              spatstat.univar_3.0-0  parallel_4.4.1         dotCall64_1.1-1        listenv_0.9.1         
[106] viridisLite_0.4.2      mvtnorm_1.1-2          scales_1.3.0           ggridges_0.5.6         leiden_0.4.3.1        
[111] tmvnsim_1.0-2          rlang_1.1.4            multcomp_1.4-17        mnormt_2.0.2          
```

