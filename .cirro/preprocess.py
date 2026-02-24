import json
from cirro.helpers.preprocess_dataset import PreprocessDataset


def extract_all_bams(ds):
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    # All BAMs (exclude .bam.bai)
    bams = sorted([
        f for f in df["file"].tolist()
        if f.endswith(".bam") and not f.endswith(".bam.bai")
    ])

    if not bams:
        raise ValueError("No BAMs found in dataset")

    return bams


def main():
    ds = PreprocessDataset.from_running()

    print("=== ds.files preview ===")
    print(ds.files.head(30).to_string(index=False))

    bams = extract_all_bams(ds)

    ds.add_param("bams", bams)

    print("\nFinal parameters:")
    print(json.dumps(ds.params, indent=2, default=str))


if __name__ == "__main__":
    main()