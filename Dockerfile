FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends libgomp1 \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /adapter
COPY pyproject.toml ./
COPY src ./src
RUN pip install --no-cache-dir .
ENV OMP_NUM_THREADS=1
# SimRNA and its data are supplied by a read-only local mount, never by this image.
ENTRYPOINT ["python", "-m", "simrna_adapter.cli"]
