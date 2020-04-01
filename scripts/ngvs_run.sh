#!/bin/bash

# IMAGE="bucket.canfar.net/gem2caom2"
COLLECTION="ngvs"
IMAGE="${COLLECTION}"

# echo "Get a proxy certificate"
# cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get image ${IMAGE}"
# docker pull ${IMAGE} || exit $?

echo "Run image ${IMAGE}"
docker run -m=7g --rm --name ${COLLECTION}_run -v ${PWD}:/usr/src/app/ ${IMAGE} ${COLLECTION}_run || exit $?

date
exit 0
