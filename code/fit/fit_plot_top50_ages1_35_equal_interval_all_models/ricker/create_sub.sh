#!/bin/bash
cd '/data/home/shpli3/R_projects/BSS_exclude_tree_raw/code/fit/fit_plot_top50_ages1_35_equal_interval_model_comparison/ricker'

for i in {2..10}; do
    touch "${i}.fit.r"
done

# 源文件名，包含要复制的代码
src_file="1.fit.r"

# 循环从2到10，将源文件的内容复制到每个目标文件中，并修改 i 的值
for i in {2..10}; do
    # 创建目标文件名
    dest_file="${i}.fit.r"
    
    # 使用 sed 复制并替换 i 的值，然后将结果输出到目标文件
    sed "s/plot_list\[,1\]/plot_list\[,${i}\]/g" "$src_file" > "$dest_file"
done

# Base script content
base_script='#!/bin/sh -x
#PBS -l nodes=1:ppn=20
#PBS -N rickerNUMBER
#PBS -l walltime=30000:00:00       
#PBS -m abe
#PBS -M 1294794885@qq.com

cd '/data/home/shpli3/R_projects/BSS_exclude_tree_raw/code/fit/fit_plot_top50_ages1_35_equal_interval_model_comparison/ricker'
date
hostname
/data/apps/R-4.2.2/lib64/R/bin/R CMD BATCH  NUMBER.fit.r NUMBER.fit.rout
date'

# Loop to create scripts
for i in $(seq 1 10)
do
    script_content="${base_script//NUMBER/$i}"
    echo "$script_content" > "${i}.ricker.sh"
done

for i in $(seq 1 10)
do
    qsub "${i}.ricker.sh"
done


for i in $(seq 38041 38047)
do
    bkill "${i}"
done

cd '/data/home/shpli3/R_projects/BSS_exclude_tree_raw/fit_results/plot_ages1_35_top50_equal_interval_model_comparison/ricker'
find summary -type f -name "*.rdata" | wc -l
ls -l -tr 'posterior'
ls -l -tr 'summary'

cd '/data/home/shpli3/R_projects/BSS_exclude_tree_raw/code/fit/fit_plot_top50_ages1_35_equal_interval_model_comparison/ricker'
file 8.fit.rout
tail -n 1000 8.fit.rout
tail -n 1000 8.fit.rout
grep "hello" '/data/home/shpli3/R_projects/BSS_exclude_tree_raw/code/fit/fit_plot_top50_ages1_35_equal_interval_model_comparison/ricker/9.ricker.sh'