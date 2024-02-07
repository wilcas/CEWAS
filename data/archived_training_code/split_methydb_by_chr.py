import os
import pandas as pd
import sqlalchemy as db


def get_predictdb_tables(chrom, dbname):
    engine = db.create_engine(f"sqlite:///{dbname}")
    conn = engine.connect()
    meta = db.MetaData()
    weights = db.Table("weights", meta, autoload=True, autoload_with=engine)
    extra = db.Table("extra", meta, autoload=True, autoload_with=engine)
    gene_pos = db.Table("gene_pos", meta, autoload=True, autoload_with=engine)
    genes = db.select([gene_pos.columns.gene]).where(gene_pos.columns.chr == chrom)
    tmp_weights = db.select([weights]).where(weights.columns.gene.in_(genes))
    tmp_extra = db.select([extra]).where(extra.columns.gene.in_(genes))
    return pd.read_sql(tmp_extra, conn), pd.read_sql(tmp_weights, conn)


def write_predictdb(dbname, extra, weights):
    engine = db.create_engine(f"sqlite:///{dbname}")
    conn = engine.connect()
    meta = db.MetaData()
    extra.to_sql("extra", conn, index=False)
    weights.to_sql("weights", conn, index=False)
    meta.reflect(bind=engine)
    extra_table = meta.tables["extra"]
    weights_table = meta.tables["weights"]
    extra_gene = db.Index("extra_gene", extra_table.columns.gene)
    extra_gene.create(bind=engine)
    weights_gene = db.Index("weights_gene", weights_table.columns.gene)
    weights_gene.create(bind=engine)
    weights_rsid = db.Index("weights_rsid", weights_table.columns.rsid)
    weights_rsid.create(bind=engine)
    weights_rsid_gene = db.Index(
        "weights_rsid_gene", weights_table.columns.rsid, weights_table.columns.gene
    )
    weights_rsid_gene.create(bind=engine)
    return 0


if __name__ == "__main__":
    dbname = "/zfs3/users/william.casazza/william.casazza/multimetaxcan/cewas_replication/jaffe_mQTL.db"
    for i in range(1, 23):
        new_db = f"/zfs3/users/william.casazza/william.casazza/multimetaxcan/cewas_replication/jaffe_mQTL_by_chr/chr{i}.db"
        if not os.path.exists(new_db):
            tmp_extra, tmp_weights = get_predictdb_tables(i, dbname)
            write_predictdb(new_db, tmp_extra, tmp_weights)
