def modify_features_file(file_path: str) -> None:
    """Mofile the "<prefix>_features.tsv.gz" file.

    Modify the "<prefix>_features.tsv.gz" file to have a third column with the value
    "Gene Expression".

    Args:
        file_path: The path to the "<prefix>_features.tsv.gz" file.
    """
    with gzip.open(filename=file_path, mode="rt") as input_file:
        lines = input_file.readlines()
        if not lines[0].strip().endswith("Gene Expression"):
            with tempfile.NamedTemporaryFile() as temp_file:
                temp_file_path = temp_file.name
                with gzip.open(filename=temp_file_path, mode="wt") as output_file:
                    for line in lines:
                        line = line.strip() + "\tGene Expression\n"
                        output_file.write(line)
                shutil.copy2(src=temp_file_path, dst=file_path)


def download_and_extract_tar_file(url: str, dir_path: str) -> None:
    """Download and extract a TAR file.

    Args:
        url: The URL of the TAR file.
        dir_path: The directory path where the TAR file will be stored and extracted.
    """
    # Download the TAR file.
    tar_file_path = os.path.join(dir_path, "data.tar")
    download_file(url=url, path=tar_file_path)

    # Extract the TAR file.
    with tarfile.open(name=tar_file_path, mode="r") as tar_file:
        tar_file.extractall(path=dir_path, filter=None)
