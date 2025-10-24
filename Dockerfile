FROM opencadc/astropy:3.8-slim

RUN apt-get update -y && apt-get install -y \
    build-essential \
    git && \
    rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

RUN pip install cadcdata \
    cadctap \
    caom2 \
    caom2repo \
    ftputil \
    importlib-metadata \
    python-dateutil \
    PyYAML \
    spherical-geometry \
    vos

WORKDIR /usr/src/app

ARG COLLECTION_BRANCH=master
ARG COLLECTION_REPO=opencadc
ARG OPENCADC_BRANCH=master
ARG OPENCADC_REPO=opencadc-metadata-curation
ARG PIPE_BRANCH=master
ARG PIPE_REPO=opencadc

RUN git clone https://github.com/${OPENCADC_REPO}/caom2tools.git --branch ${OPENCADC_BRANCH} && \
    pip install ./caom2tools/caom2utils

RUN pip install git+https://github.com/${PIPE_REPO}/caom2pipe@${PIPE_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${COLLECTION_REPO}/cfhtProc2caom2@${COLLECTION_BRANCH}#egg=cfhtProc2caom2

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

