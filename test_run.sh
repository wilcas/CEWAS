#!/bin/bash

python cewas.py --sum_stats="data/SCZ2019_refmt.txt" \
    --output="test_results.txt"\
    --mapping="data/test_data/ensemblHG19_chr22.tsv"\
    --cov_file_snp_fmt="data/test_data/ROSMAP_CEWAS_snp_methy_cov"\
    --cov_file_methy="data/test_data/ROSMAP_CEWAS_methy_expr_cov.txt.gz"\
    --predict_db_methy_fmt="data/test_data/ROSMAP_CEWAS_snp_methy"\
    --predict_db_expr="data/test_data/ROSMAP_CEWAS_methy_expr.db"\
    --or_column="OR"\
    --pvalue_column="P"\
    --se_column="SE"