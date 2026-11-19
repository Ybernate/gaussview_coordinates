# MolGauss — containerized Streamlit app.
#
# RDKit and Open Babel ship manylinux wheels for CPython 3.11, so a slim Python
# base image is enough (no conda required). Build and run:
#
#   docker build -t molgauss .
#   docker run --rm -p 8501:8501 molgauss
#
# then open http://localhost:8501

FROM python:3.11-slim

# libxrender1 / libxext6 / libsm6 are runtime libs some RDKit/Open Babel wheels
# link against; curl is used only by the container healthcheck.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libxrender1 libxext6 libsm6 libgomp1 curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install dependencies first so Docker caches this layer across code changes.
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# Application code.
COPY molgauss ./molgauss
COPY app.py .
COPY tests ./tests
COPY README.md .

# Streamlit serves on 8501; bind to all interfaces inside the container.
EXPOSE 8501
ENV STREAMLIT_SERVER_PORT=8501 \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD curl -fsS http://localhost:8501/_stcore/health || exit 1

CMD ["streamlit", "run", "app.py"]
