FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/bioinf-tools/:$PATH
ENV LANG=C.UTF-8
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

ARG MYKROBE_DIR=/mykrobe
COPY . $MYKROBE_DIR
RUN $MYKROBE_DIR/ci/install_mykrobe_linux.sh $MYKROBE_DIR
CMD mykrobe
