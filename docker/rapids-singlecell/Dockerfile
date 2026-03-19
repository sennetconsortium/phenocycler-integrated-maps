FROM ghcr.io/scverse/rapids-singlecell-cu12:v0.15.0rc4

ENV DEBIAN_FRONTEND=noninteractive

# Get CUDA toolkit
RUN apt-get update && apt-get install -y wget && rm -rf /var/lib/apt/lists/*
RUN wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-ubuntu2404.pin
RUN mv cuda-ubuntu2404.pin /etc/apt/preferences.d/cuda-repository-pin-600
RUN wget https://developer.download.nvidia.com/compute/cuda/12.9.0/local_installers/cuda-repo-ubuntu2404-12-9-local_12.9.0-575.51.03-1_amd64.deb
RUN dpkg -i cuda-repo-ubuntu2404-12-9-local_12.9.0-575.51.03-1_amd64.deb
RUN cp /var/cuda-repo-ubuntu2404-12-9-local/cuda-*-keyring.gpg /usr/share/keyrings/
RUN apt-get update \
    && apt-get -y install cuda-toolkit-12-9 \
    && rm -rf /var/cache/apt/*

COPY bin/ /opt

WORKDIR /opt

CMD ["/bin/bash"]