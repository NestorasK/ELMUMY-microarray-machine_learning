import glob, os


def find_files(directory_path, filename_pattern, verbose=True):
    """Find matching files to a pattern in a directory"""

    # Use glob to find files that match the specified pattern
    matching_files = glob.glob(
        os.path.join(directory_path, "**", filename_pattern), recursive=True
    )

    # Print the list of matching files
    if verbose:
        print("Matching Files:")
        for file_path in matching_files:
            print(file_path)

    return matching_files
