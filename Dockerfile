# Reproducible environment for the nucleoside-analogue analysis.
#   docker build -t nucleoside-analogues .
#   docker run --rm -v "$PWD:/work" nucleoside-analogues pytest -q
FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1 \
    UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

WORKDIR /work

# Dependency layer, cached independently of the source.
COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --no-install-project --group dev

COPY . .
RUN uv sync --frozen --group dev

ENV PATH="/work/.venv/bin:$PATH"
CMD ["pytest", "-q"]
