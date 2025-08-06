#!/bin/bash

# Set main directory containing all the subfolders
BASEDIR="data/max_fit_data"

# List of subfolders inside $BASEDIR
subfolders=("genotype" "phenotype" "ind")

# Loop through each subfolder
for sub in "${subfolders[@]}"; do
  subdir="${BASEDIR}/${sub}"

  # Loop over all .zip files in this subfolder
  for zipfile in "$subdir"/*.zip; do
    [[ -e "$zipfile" ]] || continue  # skip if no zip files

    zipname=$(basename "$zipfile")
    foldername="${zipname%.zip}"
    targetdir="${subdir}/${foldername}"

    # Skip if the folder already exists
    if [[ -d "$targetdir" ]]; then
      echo "⏩ Skipping $zipname — folder already exists."
      continue
    fi

    echo "📦 Unzipping $zipname in $subdir..."

    # Create a temporary directory to unzip into
    tmpdir="${subdir}/__tmp_unzip__"
    rm -rf "$tmpdir"
    mkdir -p "$tmpdir"

    # Unzip into the temporary directory
    unzip -q "$zipfile" -d "$tmpdir"

    # Detect the top-level folder or file in the zip
    entries=("$tmpdir"/*)
    if [[ ${#entries[@]} -eq 1 && -d "${entries[0]}" ]]; then
      # Single folder — rename it
      mv "${entries[0]}" "$targetdir"
      echo "✔ Renamed extracted folder to ${foldername}"
    else
      # Multiple files or complex structure — move all into a folder
      mkdir -p "$targetdir"
      mv "$tmpdir"/* "$targetdir/"
      echo "⚠ Multiple items extracted — moved into folder ${foldername}"
    fi

    # Clean up temp directory
    rm -rf "$tmpdir"
  done
done
