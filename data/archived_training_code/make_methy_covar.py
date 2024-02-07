import gzip
import numpy as np
import os
import pandas as pd
import sqlalchemy as db


def read_methy(methy_f):
    with open(methy_f, "r") as f:
        line = f.readline().strip().split("\t")
        ncol = len(line)
        probenames = line[2:]
    data = np.loadtxt(methy_f, usecols=range(2, ncol), skiprows=1, delimiter="\t")
    return pd.DataFrame(data=data, columns=probenames)


def predicted_methy_covariance(methy_f, predict_db, outfile):
    df = read_methy(methy_f)
    engine = db.create_engine(f"sqlite:///{predict_db}")
    connection = engine.connect()
    metadata = db.MetaData()
    weights = db.Table("weights", metadata, autoload=True, autoload_with=engine)
    query = db.select([weights.columns.rsid, weights.columns.gene]).distinct()
    db_df = pd.read_sql(query, engine)
    for (grp, dat) in db_df.groupby("gene"):
        cur_df = df[dat.rsid[dat.rsid.isin(df.columns)]]
        if not cur_df.empty:
            cur_cov = cur_df.cov().stack().reset_index()
            cur_cov.insert(0, "gene", grp)
            cur_cov.columns = ["GENE", "RSID1", "RSID2", "VALUE"]
            cur_cov.to_csv(
                outfile,
                header=not os.path.isfile(outfile),
                index=False,
                sep=" ",
                mode="a",
                compression="gzip",
            )


if __name__ == "__main__":
    methy_f = "/zfs3/scratch/william.casazza/ROSMAP_matrixeqtl/ROSMAP_CEWAS_snp_methy_predicted_expression.txt"
    predict_db = "/zfs3/users/william.casazza/william.casazza/multimetaxcan/ROSMAP_CEWAS_methy_expr_nonzero.db"
    outfile = "/zfs3/users/william.casazza/william.casazza/multimetaxcan/cewas_replication/rosmap_imputed_methy_expr_model_par_nonzero_cov.txt.gz"
    predicted_methy_covariance(methy_f, predict_db, outfile)

    # methy_f = "/zfs3/scratch/william.casazza/CMC_alt_fmt/robust_across_all_values_methy_rosmap_nonan.txt"
    # predict_db = "/zfs3/users/william.casazza/william.casazza/multimetaxcan/cewas_replication/cmc_imputed_methy_expr_model_par_nonzero.db"
    # outfile = "/zfs3/users/william.casazza/william.casazza/multimetaxcan/cewas_replication/cmc_imputed_methy_expr_model_par_nonzero_cov.txt.gz"
    # predicted_methy_covariance(methy_f,predict_db,outfile)
