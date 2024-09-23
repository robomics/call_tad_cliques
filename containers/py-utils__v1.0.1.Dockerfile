# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.5.10-noble AS builder

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install -y \
        -c conda-forge \
        -c bioconda \
        'python>=3.12' \
        'bioframe>0.7' \
        'pandas>2' \
        procps-ng \
&& micromamba clean --all -y

WORKDIR /data

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/robomics/call_tad_cliques'
LABEL org.opencontainers.image.documentation='https://github.com/robomics/call_tad_cliques'
LABEL org.opencontainers.image.source='https://github.com/robomics/call_tad_cliques'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-py-utils}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
