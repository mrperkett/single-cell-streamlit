"""
A script to download example data for the web application
"""

import pooch


def main():

    # download data to ~/.cache/scverse_tutorials
    base_url = "doi:10.6084/m9.figshare.22716739.v1/"
    dir_path = "~/.cache/scverse_tutorials"
    example_data = pooch.create(
        path=dir_path,
        base_url=base_url,
    )
    example_data.load_registry_from_doi()

    for filename in example_data.registry_files:
        print(f"Attempting to fetch '{filename}`")
        example_data.fetch(filename)


if __name__ == "__main__":
    main()
