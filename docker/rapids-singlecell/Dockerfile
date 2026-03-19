FROM ghcr.io/scverse/rapids-singlecell-cu12:v0.15.0rc4

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    autoconf \
    automake \
    bison \
    fastqc \
    flex \
    git \
    libfreetype6-dev \
    libxml2-dev \
    libtool \
    m4 \
    zlib1g-dev \
    software-properties-common \
    pkg-config \
 && rm -rf /var/cache/apt/*

RUN add-apt-repository ppa:deadsnakes/ppa

RUN apt-get -y install \
    python3.9 \
    python3.9-dev \
    python3-pip \
    python3.9-venv \
 && rm -rf /var/cache/apt/*

## Python packages
RUN python3.9 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
COPY docker/rapids-singlecell/requirements.txt /opt/requirements.txt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

COPY bin/ /opt

WORKDIR /opt

CMD ["/bin/bash"]