# Single stage: Runtime image with pre-built binary
# Build the binary locally first with: cargo build --release
FROM ubuntu:22.04

LABEL maintainer="Starfish-Bioscience" \
      version="0.9.0" \
      description="CoverM with CRAM support - Full dependencies"

ENV DEBIAN_FRONTEND=noninteractive

# Install base dependencies from apt (only what's not available/better in conda)
RUN apt-get update && \
    apt-get install -y \
        wget \
        curl \
        ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Miniforge (no TOS issues)
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

ENV PATH=/opt/conda/bin:$PATH

# Configure conda and install bioconda packages
# Note: samtools must be from conda (apt version 1.13 has non-UTF-8 output that crashes coverm)
# Note: dashing excluded due to zlib dependency conflicts
RUN conda config --add channels bioconda && \
    conda config --set channel_priority strict && \
    conda install -y \
        samtools=1.23 \
        minimap2=2.30 \
        bwa \
        bwa-mem2=2.3 \
        skani \
        fastani \
        strobealign && \
    conda clean -a -y

# Copy pre-built binary (must be built locally first)
COPY target/release/coverm /usr/local/bin/

RUN chmod +x /usr/local/bin/coverm

# Verify installation
RUN coverm --version

WORKDIR /data

ENTRYPOINT ["coverm"]
CMD ["--help"]
