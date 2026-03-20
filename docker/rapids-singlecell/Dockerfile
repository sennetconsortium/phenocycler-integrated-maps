FROM ghcr.io/scverse/rapids-singlecell-cu12:v0.15.0rc4

ENV DEBIAN_FRONTEND=noninteractive

RUN pip install leidenalg
RUN pip install igraph

COPY bin/ /opt

WORKDIR /opt

CMD ["/bin/bash"]