#!/bin/bash

IMAGE="opencadc/cfhtProc2caom2"
CMD="cfht_proc"
# echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get image ${IMAGE}"
docker pull ${IMAGE} || exit $?

echo "Run image ${IMAGE}"
docker run --rm --name ${CMD}_run -v ${PWD}:/usr/src/app/ ${IMAGE} ${CMD}_run || exit $?

date
exit 0
