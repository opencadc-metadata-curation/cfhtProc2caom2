#!/bin/bash

# IMAGE="opencadc/ngvs2caom2"
COLLECTION="ngvs"
IMAGE="${COLLECTION}"

# echo "Get a proxy certificate"
# cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get image ${IMAGE}"
# docker pull ${IMAGE} || exit $?

echo "Run image ${IMAGE}"
docker run --rm --name ${COLLECTION}_run -v ${PWD}:/usr/src/app/ ${IMAGE} ${COLLECTION}_run || exit $?

date
exit 0
