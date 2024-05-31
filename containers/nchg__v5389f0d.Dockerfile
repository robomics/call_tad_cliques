# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

ARG CONTAINER_VERSION

FROM ghcr.io/paulsengroup/nchg:sha-5389f0d AS base

ARG CONTAINER_TITLE
ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN apt-get update \
&&  apt-get install -y pigz procps zstd \
&& rm -rf /var/lib/apt/lists/*


ENTRYPOINT []
WORKDIR /data

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/robomics/call_tad_cliques'
LABEL org.opencontainers.image.documentation='https://github.com/robomics/call_tad_cliques'
LABEL org.opencontainers.image.source='https://github.com/robomics/call_tad_cliques'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-nchg}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
