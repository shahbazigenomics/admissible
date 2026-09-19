# syntax=docker/dockerfile:1
# ---------------------------------------------------------------------------
# admissible — a single image that both runs the CLI and serves the project
# website.
#
#   docker build -t admissible .
#   docker run --rm admissible audit family/*.vcf --ped family.ped   # the tool
#   docker run --rm -p 8000:8000 admissible                          # the site
#
# The package has zero required runtime dependencies by design, so there is
# nothing to compile — installing it is just copying the source and its
# metadata into place.
# ---------------------------------------------------------------------------
FROM python:3.14-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PORT=8000

WORKDIR /app

# 1 ─ the package itself (metadata first so this layer caches on source-only
#     changes to web/ or docker/).
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install .

# 2 ─ the static site and the stdlib server that hosts it.
COPY web ./web
COPY docker ./docker
RUN chmod +x docker/entrypoint.sh docker/serve.py

# 3 ─ drop privileges.
RUN useradd --system --uid 10001 admissible && chown -R admissible /app
USER admissible

EXPOSE 8000
ENTRYPOINT ["/app/docker/entrypoint.sh"]
CMD ["serve"]
