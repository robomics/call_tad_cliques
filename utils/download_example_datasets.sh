#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

set -e
set -u

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path_to_output_folder/"
  exit 1
fi

if ! command -v curl &> /dev/null; then
  2>&1 echo "Unable to find curl in your PATH"
  exit 1
fi

if ! command -v sha256sum &> /dev/null; then
  2>&1 echo "Unable to find sha256sum in your PATH"
  exit 1
fi


urls=('https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/8b9d6836-0d6b-4c2f-9eaa-323f4fd7b6e4/4DNFI74YHN5W.mcool'
      'https://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz'
      'https://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/gap.txt.gz')

checksums=('a5e3dc9e9e6b68b577c086dcbf34e8a21094bf34fa27461440642e9273ee228b'
           '30513728ad2afc50fc75e0ba1d48e5602c491d5b7313505bb5c5e2a9c67bd2e5'
           'fc9da5605c3458e45925eb344c37e8584f13613b00ddf91ba11d86f5bac6ec98')

outdir="$1"
mkdir -p "$outdir"

function download_dataset {
  src="$1"
  dest="${outdir%/}/$(basename "$src")"
  checksum="$2"

  if [ -f "$dest" ]; then
    if sha256sum -c <(echo "$checksum  $dest") &> /dev/null; then
      2>&1 echo "File $(basename "$dest") already downloaded. SKIPPING!"
      return
    fi
  fi

  echo "curl -L \"$src\" -o \"$dest\""
  curl -L "$src" -o "$dest"
  sha256sum -c <(echo "$checksum  $dest")
}

for i in "${!urls[@]}"; do
  download_dataset "${urls[$i]}" "${checksums[$i]}"
done
