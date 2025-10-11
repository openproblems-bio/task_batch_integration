import os
import tarfile
import tempfile
import zipfile

def unpack_directory(directory):
    """
    Unpack a directory to a temporary location (if needed)

    Args:
        directory (str): Path to a directory, .zip, or .tar.gz file.

    Returns:
        tuple: (unpacked_directory (str), temp_directory (TemporaryDirectory or None))
            unpacked_directory: Path to the unpacked directory.
            temp_directory: TemporaryDirectory object if a temp dir was created, else None.
    """
    print(f"Unpacking directory: '{directory}'", flush=True)

    if os.path.isdir(directory):
        print(f"Returning provided directory: '{directory}'", flush=True)
        temp_directory = None
        unpacked_directory = directory
    else:
        temp_directory = tempfile.TemporaryDirectory()
        unpacked_directory = temp_directory.name

        if zipfile.is_zipfile(directory):
            print("Extracting .zip...", flush=True)
            with zipfile.ZipFile(directory, "r") as zip_file:
                zip_file.extractall(unpacked_directory)
        elif tarfile.is_tarfile(directory) and directory.endswith(".tar.gz"):
            print("Extracting .tar.gz...", flush=True)
            with tarfile.open(directory, "r:gz") as tar_file:
                tar_file.extractall(unpacked_directory)
                unpacked_directory = os.path.join(unpacked_directory, os.listdir(unpacked_directory)[0])
        else:
            raise ValueError(
                "The 'directory' argument should be a directory, a .zip file or a .tar.gz file"
            )
        print(f"Extracted to '{unpacked_directory}'", flush=True)

    return (unpacked_directory, temp_directory)
