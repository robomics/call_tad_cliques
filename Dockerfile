# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.4.0-jammy AS builder

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ARG NCHG_GIT_URL='https://github.com/Chrom3D/preprocess_scripts.git'
ARG NCHG_GIT_COMMIT='ca149e12cf929e90a74a779fd7e7f24a0101fb97'

RUN micromamba install -y \
        -c conda-forge \
        cxx-compiler \
        git \
        make \
        unzip

RUN git clone "$NCHG_GIT_URL" /tmp/nchg && cd /tmp/nchg \
&& git checkout "$NCHG_GIT_COMMIT" \
&& unzip /tmp/nchg/NCHG_hic.zip -d /tmp/ \
&& sed -i 's/g++/c++/g' '/tmp/NCHG_hic/Makefile' \
&& (cd /tmp/NCHG_hic && make -j $(nproc)) \
&& install -Dm755 /tmp/NCHG_hic/NCHG /tmp/NCHG \
&& /tmp/NCHG --help


FROM mambaorg/micromamba:1.4.0-jammy AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/

RUN micromamba install -y \
        -c conda-forge \
        -c bioconda \
        --file /tmp/env.yml \
&& micromamba clean --all -y \
&& rm /tmp/env.yml

COPY --from=builder --chown=root:root /tmp/NCHG /usr/local/bin/NCHG

USER root
RUN chmod 755 /usr/local/bin/NCHG
USER mambauser

WORKDIR /data

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/robomics/call_tad_cliques'
LABEL org.opencontainers.image.documentation='https://github.com/robomics/call_tad_cliques'
LABEL org.opencontainers.image.source='https://github.com/robomics/call_tad_cliques'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-call-tad-cliques}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
