#!/bin/sh -x
#PBS -l nodes=1:ppn=20
#PBS -N ricker7
#PBS -l walltime=30000:00:00       
#PBS -m abe
#PBS -M 1294794885@qq.com

cd /data/home/shpli3/R_projects/BSS_exclude_tree_raw/code/fit/fit_plot_top50_ages1_35_equal_interval_model_comparison/ricker
date
hostname
/data/apps/R-4.2.2/lib64/R/bin/R CMD BATCH  7.fit.r 7.fit.rout
date
