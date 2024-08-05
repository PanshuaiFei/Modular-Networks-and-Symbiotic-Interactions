#!/bin/bash

# 设置起始目录和结束目录编号

    cd 0

    # 创建目录
    yes | mkdir bootstrap_counts bootstrap_correlation

    # 运行 fastspar
    yes | fastspar --otu_table 0.tsv --correlation median_correlation.tsv --covariance median_covariance.tsv --threads 20

    # 运行 fastspar_bootstrap
    yes | fastspar_bootstrap --otu_table 0.tsv --number 1000 --prefix bootstrap_counts/fake_data --threads 20

    # 并行运行 fastspar
    yes | parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*

    # 运行 fastspar_pvalues
    yes | fastspar_pvalues --otu_table 0.tsv --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_fake_data_ --permutations 1000 --outfile pvalues.tsv --threads 20

    #删除过程性文件
    rm -rf bootstrap_correlation bootstrap_counts
    # 切换回上级目录
    cd ..




# 设置起始目录和结束目录编号

    cd healthy

    # 创建目录
    yes | mkdir bootstrap_counts bootstrap_correlation

    # 运行 fastspar
    yes | fastspar --otu_table healthy.tsv --correlation median_correlation.tsv --covariance median_covariance.tsv --threads 20

    # 运行 fastspar_bootstrap
    yes | fastspar_bootstrap --otu_table healthy.tsv --number 1000 --prefix bootstrap_counts/fake_data --threads 20

    # 并行运行 fastspar
    yes | parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*

    # 运行 fastspar_pvalues
    yes | fastspar_pvalues --otu_table healthy.tsv --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_fake_data_ --permutations 1000 --outfile pvalues.tsv --threads 20

    #删除过程性文件
    rm -rf bootstrap_correlation bootstrap_counts
    # 切换回上级目录
    cd ..


# 切换到目录
    cd HS

    # 创建目录
    yes | mkdir bootstrap_counts bootstrap_correlation

    # 运行 fastspar
    yes | fastspar --otu_table HS.tsv --correlation median_correlation.tsv --covariance median_covariance.tsv --threads 20

    # 运行 fastspar_bootstrap
    yes | fastspar_bootstrap --otu_table HS.tsv --number 1000 --prefix bootstrap_counts/fake_data --threads 20

    # 并行运行 fastspar
    yes | parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*

    # 运行 fastspar_pvalues
    yes | fastspar_pvalues --otu_table HS.tsv --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_fake_data_ --permutations 1000 --outfile pvalues.tsv --threads 20

    #删除过程性文件
    rm -rf bootstrap_correlation bootstrap_counts
    # 切换回上级目录
    cd ..



# 切换到目录
    cd I

    # 创建目录
    yes | mkdir bootstrap_counts bootstrap_correlation

    # 运行 fastspar
    yes | fastspar --otu_table I.tsv --correlation median_correlation.tsv --covariance median_covariance.tsv --threads 20

    # 运行 fastspar_bootstrap
    yes | fastspar_bootstrap --otu_table I.tsv --number 1000 --prefix bootstrap_counts/fake_data --threads 20

    # 并行运行 fastspar
    yes | parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*

    # 运行 fastspar_pvalues
    yes | fastspar_pvalues --otu_table I.tsv --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_fake_data_ --permutations 1000 --outfile pvalues.tsv --threads 20

    #删除过程性文件
    rm -rf bootstrap_correlation bootstrap_counts
    # 切换回上级目录
    cd ..



# 切换到目录
    cd II

    # 创建目录
    yes | mkdir bootstrap_counts bootstrap_correlation

    # 运行 fastspar
    yes | fastspar --otu_table II.tsv --correlation median_correlation.tsv --covariance median_covariance.tsv --threads 20

    # 运行 fastspar_bootstrap
    yes | fastspar_bootstrap --otu_table II.tsv --number 1000 --prefix bootstrap_counts/fake_data --threads 20

    # 并行运行 fastspar
    yes | parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*

    # 运行 fastspar_pvalues
    yes | fastspar_pvalues --otu_table II.tsv --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_fake_data_ --permutations 1000 --outfile pvalues.tsv --threads 20

    #删除过程性文件
    rm -rf bootstrap_correlation bootstrap_counts
    # 切换回上级目录
    cd ..


# 切换到目录
    cd III

    # 创建目录
    yes | mkdir bootstrap_counts bootstrap_correlation

    # 运行 fastspar
    yes | fastspar --otu_table III.tsv --correlation median_correlation.tsv --covariance median_covariance.tsv --threads 20

    # 运行 fastspar_bootstrap
    yes | fastspar_bootstrap --otu_table III.tsv --number 1000 --prefix bootstrap_counts/fake_data --threads 20

    # 并行运行 fastspar
    yes | parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*

    # 运行 fastspar_pvalues
    yes | fastspar_pvalues --otu_table III.tsv --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_fake_data_ --permutations 1000 --outfile pvalues.tsv --threads 20

    #删除过程性文件
    rm -rf bootstrap_correlation bootstrap_counts
    # 切换回上级目录
    cd ..




# 切换到目录
    cd IV

    # 创建目录
    yes | mkdir bootstrap_counts bootstrap_correlation

    # 运行 fastspar
    yes | fastspar --otu_table IV.tsv --correlation median_correlation.tsv --covariance median_covariance.tsv --threads 20

    # 运行 fastspar_bootstrap
    yes | fastspar_bootstrap --otu_table IV.tsv --number 1000 --prefix bootstrap_counts/fake_data --threads 20

    # 并行运行 fastspar
    yes | parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*

    # 运行 fastspar_pvalues
    yes | fastspar_pvalues --otu_table IV.tsv --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_fake_data_ --permutations 1000 --outfile pvalues.tsv --threads 20

    #删除过程性文件
    rm -rf bootstrap_correlation bootstrap_counts
    # 切换回上级目录
    cd ..




# 切换到目录
    cd MP

    # 创建目录
    yes | mkdir bootstrap_counts bootstrap_correlation

    # 运行 fastspar
    yes | fastspar --otu_table MP.tsv --correlation median_correlation.tsv --covariance median_covariance.tsv --threads 20

    # 运行 fastspar_bootstrap
    yes | fastspar_bootstrap --otu_table MP.tsv --number 1000 --prefix bootstrap_counts/fake_data --threads 20

    # 并行运行 fastspar
    yes | parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*

    # 运行 fastspar_pvalues
    yes | fastspar_pvalues --otu_table MP.tsv --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_fake_data_ --permutations 1000 --outfile pvalues.tsv --threads 20

    #删除过程性文件
    rm -rf bootstrap_correlation bootstrap_counts
    # 切换回上级目录
    cd ..


    cd total

    # 创建目录
    yes | mkdir bootstrap_counts bootstrap_correlation

    # 运行 fastspar
    yes | fastspar --otu_table total.tsv --correlation median_correlation.tsv --covariance median_covariance.tsv --threads 20

    # 运行 fastspar_bootstrap
    yes | fastspar_bootstrap --otu_table total.tsv --number 1000 --prefix bootstrap_counts/fake_data --threads 20

    # 并行运行 fastspar
    yes | parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*

    # 运行 fastspar_pvalues
    yes | fastspar_pvalues --otu_table total.tsv --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_fake_data_ --permutations 1000 --outfile pvalues.tsv --threads 20

    #删除过程性文件
    rm -rf bootstrap_correlation bootstrap_counts
    # 切换回上级目录
    cd ..



