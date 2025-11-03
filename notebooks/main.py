import polars as pl
from pathlib import Path

def main():
    fastq_files = Path("/zata/zippy/ramirezc/atac-smk/.test").glob("*.fastq.gz")
    R1, R2 = [], []

    for file in fastq_files:
        if file.name.endswith("R1.fastq.gz"):
            R1.append(str(file))
        elif file.name.endswith("R2.fastq.gz"):
            R2.append(str(file))
    R1, R2 = sorted(R1), sorted(R2)

    df = pl.DataFrame({"R1": R1, "R2": R2})
    with_dataset = df.with_columns(pl.col("R1").str.split("_").list.get(1).alias("Sample"))
    print(with_dataset)
    print(with_dataset.columns)

    with_dataset.write_csv("/zata/zippy/ramirezc/atac-smk/.test/metadata.tsv", separator="\t")

if __name__ == "__main__":
    main()
