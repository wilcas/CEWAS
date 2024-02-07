import gzip
import pandas as pd

if __name__ == "__main__":
    coords = pd.read_csv("~/jaffe_mQTL_coord.csv").drop_duplicates()
    header = b"GENE RSID1 RSID2 VALUE\n"
    for i in range(1, 22):
        chr_i_cpg = coords["cpg"][coords["chr"] == str(i)].to_numpy()
        print(chr_i_cpg)
        with gzip.open("jaffe_mQTL.cov.txt.gz", "r") as f:
            for line in f:

                if line.split(b" ")[0].decode() in chr_i_cpg:
                    print("here")
                    with gzip.open(f"mQTL_cov_by_chr/chr{i}.cov.txt.gz", "ab") as out_f:
                        if out_f.tell() == 0:
                            out_f.write(header)
                        out_f.write(line)
