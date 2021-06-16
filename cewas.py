import argparse
import glob
import joblib
import os
import subprocess
import pandas as pd

def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--output',
        help="Name of file for CEWAS results",
        required=True
    )
    parser.add_argument(
        '--sum_stats',
        help="Path to GWAS summary statistics",
        required=True
    )

    parser.add_argument(
        '--beta_column',
        help="Name of effect size column in GWAS summary statistics",
    )
    parser.add_argument(
        '--or_column',
        help="Name of odds ratio column in GWAS summary statistics"
    )
    parser.add_argument(
        '--se_column',
        help="Name of standard error column in GWAS summary statistics"
    )
    parser.add_argument(
        '--pvalue_column',
        help="Name of P value column in GWAS summary statistics"
    )
    parser.add_argument(
        '--zscore_column',
        help="Name of z-score column in GWAS summary statistics"
    )

    parser.add_argument(
        '--mapping',
        help="ENSEMBL mart file for re-mapping gene names",
        default="data/ensemblHG19.tsv"
    )
    parser.add_argument(
        '--reference_allele',
        help="Reference allele column in GWAS summary statistics",
        default="A1"
    )
    parser.add_argument(
        '--alternative_allele',
        help="Alternative allele column in GWAS summary statistics",
        default="A2"
    )
    parser.add_argument(
        '--snp_column',
        help="SNP id column in GWAS summary statistics",
        default="SNP"
    )
    parser.add_argument(
        '--cov_file_snp_fmt',
        help="Covariance file prefix for snp to methylation model",
        default="data/covfiles/ROSMAP_CEWAS_snp_methy_cov"
    )
    parser.add_argument(
        '--cov_file_methy',
        help="Covariance file for methylation to expression model",
        default="data/covfiles/ROSMAP_CEWAS_pred_methy_expr_cov.txt.gz"
    )
    parser.add_argument(
        '--predict_db_methy_fmt',
        help="SNP to methylation predictdb prefix",
        default="data/predictdb/ROSMAP_CEWAS_snp_methy"
    )
    parser.add_argument(
        '--predict_db_expr',
        help="Methylation to expression predictdb",
        default="data/predictdb/ROSMAP_CEWAS_methy_expr_pred.db"
    )
    parser.add_argument(
        '--ncores',
        help="Number of cores available for job",
        default=1,
        type=int
    )

    return parser


def format_cewas(input_f, mapping_f, output):
    result = pd.read_csv(input_f)
    mapping = pd.read_csv(mapping_f, sep="\s+")
    merged = pd.merge(result,mapping,how="left",left_on=['gene'], right_on=['name2'])
    final = merged[["gene","gene_name","chrom","txStart","zscore","pvalue"]]
    final.columns = ["ENSG", "geneSymbol", "chr", "tss", "z", "p"]
    final = final.drop_duplicates("ENSG")
    final.to_csv(output,index=False)


def run_snp_methy_by_chr(args,chrom):
    cur_predict_db= "{}_chr{}.db".format(args.predict_db_methy_fmt,chrom)
    cur_cov_file = "{}_chr{}.txt.gz".format(args.cov_file_snp_fmt,chrom)
    cur_output = "{}.methy.chr{}".format(args.output,chrom)
    if not (os.path.exists(cur_cov_file) and os.path.exists(cur_predict_db)):
        print("SNP to methylation data is missing for chromosome {}, skipping".format(chrom))
        return 0
    if os.path.exists(cur_output):
        print("{} already exists, please remove it to recalculate".format(cur_output))
        return 0
    if args.zscore_column is not None:
        with open(os.devnull,'w') as f:
            subprocess.run(
                """
                {executable} \\
                    --model_db_path {predict_db_methy} \\
                    --gwas_file {gwas_file} \\
                    --snp_column {snp_column} \\
                    --zscore_column {zscore_column}\\
                    --non_effect_allele_column {alternative_allele}\\
                    --effect_allele_column {reference_allele}\\
                    --covariance {cov_file}\\
                    --output_file {output}\\
                    --verbosity 10\\
                    --keep_non_rsid
                """.format(
                        executable=executable,
                        predict_db_methy=cur_predict_db,
                        gwas_file=args.sum_stats,
                        snp_column=args.snp_column,
                        zscore_column=args.zscore_column,
                        alternative_allele=args.alternative_allele,
                        reference_allele=args.reference_allele,
                        cov_file=cur_cov_file,
                        output=cur_output
                ),
                shell=True,
                check=True,
                stdout=f,
                stderr=f
            )
    else:
        effect = "--or_column {}\\".format(args.or_column) if args.or_column is not None else "--beta_column {}\\".format(args.beta_column)
        with open(os.devnull,'w') as f:
            subprocess.run(
                """
                {executable} \\
                    --model_db_path {predict_db_methy} \\
                    --gwas_file {gwas_file} \\
                    {effect}
                    --snp_column {snp_column} \\
                    --pvalue_column {pvalue_column}\\
                    --se_column {se_column}\\
                    --non_effect_allele_column {alternative_allele}\\
                    --effect_allele_column {reference_allele}\\
                    --covariance {cov_file}\\
                    --output_file {output}\\
                    --verbosity 10\\
                    --keep_non_rsid
                """.format(
                        executable=executable,
                        predict_db_methy=cur_predict_db,
                        gwas_file=args.sum_stats,
                        snp_column=args.snp_column,
                        pvalue_column=args.pvalue_column,
                        effect=effect,
                        se_column=args.se_column,
                        alternative_allele=args.alternative_allele,
                        reference_allele=args.reference_allele,
                        cov_file=cur_cov_file,
                        output=cur_output
                ),
                shell=True,
                check=True,
                stdout=f,
                stderr=f
            )


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    path = os.path.dirname(os.path.realpath(__file__))
    executable = path + "/MetaXcan/software/SPrediXcan.py"
    methy_result = "{}.methy".format(args.output)
    if os.path.exists(methy_result):
        print("{} already exists please delete to recalculate. Running additional associations".format(methy_result))
        result = pd.read_csv(methy_result,sep="\t")
        result['A']='A'
        result['C']='C'
    else:
        print("Starting SNP to methylation associations by chromosome")
        joblib.Parallel(n_jobs=args.ncores, backend="multiprocessing",verbose=100)(
            joblib.delayed(run_snp_methy_by_chr)(args,i) for i in range(1,23)
        )
        #merge results for each chromosome
        results = []
        for chrom in range(1,23):
            cur_result = "{}.methy.chr{}".format(args.output,chrom)
            if os.path.exists(cur_result):
                df = pd.read_csv(cur_result)
                results.append(df)
                os.remove(cur_result)
        if len(results) > 1:
            result=pd.concat(results)
        else:
            result = df
        result.dropna(thresh=2)
        result['A']='A'
        result['C']='C'

    result.to_csv(methy_result,index=False,sep="\t")
    with open(os.devnull,'w') as f:
        subprocess.run(
            """
            {executable} \\
                --model_db_path {predict_db_expr} \\
                --gwas_file {gwas_file} \\
                --snp_column {snp_column}\\
                --zscore_column {zscore_column}\\
                --non_effect_allele_column {alternative_allele}\\
                --effect_allele_column {reference_allele}\\
                --covariance {cov_file}\\
                --output_file {output}\\
                --verbosity 10\\
                --keep_non_rsid
            """.format(
                    executable=executable,
                    predict_db_expr=args.predict_db_expr,
                    gwas_file=methy_result,
                    snp_column= "gene",
                    zscore_column="zscore",
                    alternative_allele="A",
                    reference_allele="C",
                    cov_file=args.cov_file_methy,
                    output=args.output
            ),
            shell=True,
            check=True,
            stdout=f,
            stderr=f
        )
    result.drop(columns=['A','C']).to_csv(methy_result,index=False,sep="\t")
    format_cewas(args.output,args.mapping,args.output)
