FROM opencadc/astropy:3.8-slim

RUN apt-get update -y
RUN apt-get install -y \
    build-essential \
    git

RUN pip install cadcdata \
    cadctap \
    caom2 \
    caom2repo \
    caom2utils \
    deprecated \
    ftputil \
    importlib-metadata \
    pytz \
    PyYAML \
    spherical-geometry \
    vos

WORKDIR /usr/src/app

ARG OPENCADC_BRANCH=master
ARG OPENCADC_REPO=opencadc

RUN git clone https://github.com/${OPENCADC_REPO}/caom2tools.git --branch ${OPENCADC_BRANCH} --single-branch && \
  pip install ./caom2tools/caom2utils

RUN git clone https://github.com/${OPENCADC_REPO}/caom2pipe.git && \
  pip install ./caom2pipe

RUN git clone https://github.com/${OPENCADC_REPO}/cfhtProc2caom2.git && \
  cp ./cfhtProc2caom2/scripts/config.yml / && \
  cp ./cfhtProc2caom2/scripts/docker-entrypoint.sh / && \
  pip install ./cfhtProc2caom2

RUN apt-get purge -y git

ENTRYPOINT ["/docker-entrypoint.sh"]

