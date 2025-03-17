#!/bin/sh -x
#PBS -l nodes=1:ppn=20
#PBS -N null3
#PBS -l walltime=30000:00:00       
#PBS -m abe
#PBS -M 1294794885@qq.com

cd /data/home/shpli3/R_projects/BSS_exclude_tree_raw/code/fit/fit_plot_top50_ages1_35_model_comparison/null
date
hostname
/data/apps/R-4.2.2/lib64/R/bin/R CMD BATCH  3.fit.r 3.fit.rout
date
