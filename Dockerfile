FROM python:3.11-slim

# libs de sistema para pysam/htslib com acesso remoto (libcurl)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-openssl-dev libssl-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .

CMD ["python", "run_pipeline.py"]
