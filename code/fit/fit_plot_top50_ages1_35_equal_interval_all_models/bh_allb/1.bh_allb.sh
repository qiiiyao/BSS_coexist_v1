#!/bin/sh -x
#PBS -l nodes=1:ppn=20
#PBS -N bh_allb1
#PBS -l walltime=30000:00:00       
#PBS -m abe
#PBS -M 1294794885@qq.com

cd /data/home/shpli3/R_projects/BSS_exclude_tree_raw/code/fit/fit_plot_top50_ages1_35_equal_interval_model_comparison/bh_allb
date
hostname
/data/apps/R-4.2.2/lib64/R/bin/R CMD BATCH  1.fit.r 1.fit.rout
date
