#!/bin/sh -x
#PBS -l nodes=1:ppn=20
#PBS -N mod_com
#PBS -l walltime=30000:00:00       
#PBS -m abe
#PBS -M 1294794885@qq.com

cd '/data/home/shpli3/R_projects/BSS_exclude_tree_raw/code/results analysing'
date
hostname
/data/apps/R-4.2.2/lib64/R/bin/R CMD BATCH  fit_plot_ages1_35_top50_mod_equal_interval_mod_comparison.R fit_plot_ages1_35_top50_mod_equal_interval_mod_comparison.Rout
date
