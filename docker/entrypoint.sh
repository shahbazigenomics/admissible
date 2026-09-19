#!/bin/sh
# One image, two jobs.
#
#   (default)                 -> serve the project website
#   audit / identity / ...    -> run the admissible CLI
#
# So `docker run admissible` brings up the site, while
# `docker run admissible audit family/*.vcf --ped family.ped` runs the tool.
set -e

if [ "${1:-serve}" = "serve" ]; then
    exec python /app/docker/serve.py
fi

exec admissible "$@"
