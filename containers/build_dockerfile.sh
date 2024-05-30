# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

if [ $# -ne 1 ]; then
  1>&2 echo "Usage: $0 Dockerfile"
  1>&2 echo "Example: $0 utils__v1.0.0.Dockerfile"
  exit 1
fi

git_root="$(git rev-parse --show-toplevel)"
img_prefix="$(git remote -v | grep push | sed -E 's/.*:(.*)\.git.*/\1/' | tr '/' '-')"

dockerfile="$1"
dockerfile_name="$(basename "${dockerfile[*]}" .Dockerfile)"
name="$(echo "$dockerfile_name" | sed -E 's/(.*)__.*/\1/')"
version="$(echo "$dockerfile_name" | sed -E 's/.*__v(.*)/\1/')"

sudo docker build \
  -t "$name:latest" \
  -t "$name:$version" \
  -f "$dockerfile" \
  --build-arg="CONTAINER_VERSION=$version" \
   "$git_root"

if ! command -v apptainer &> /dev/null; then
  2>&1 echo 'Unable to find apptainer in your PATH'
  2>&1 echo 'Skipping apptainer build step!'
  exit 0
fi

outdir="$git_root/containers/cache"
sif="$outdir/$img_prefix-$name-$version.img"
mkdir -p "$outdir"

sudo apptainer build -F "$sif" "docker-daemon://$name:$version"

uid="$(id -u)"
gid="$(id -g)"

sudo chown "$uid:$gid" "$sif"

set -x
apptainer exec "$sif" sh -c 'echo Hi!' > /dev/null
