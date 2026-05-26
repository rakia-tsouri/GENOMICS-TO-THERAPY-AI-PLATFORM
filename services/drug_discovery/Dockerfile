FROM python:3.10-slim

# Install system dependencies.
# Pocket detection uses P2Rank (a Java ML pocket predictor); we only need a JRE.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    ca-certificates \
    default-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Install P2Rank
WORKDIR /opt
RUN wget -q https://github.com/rdk/p2rank/releases/download/2.4/p2rank_2.4.tar.gz && \
    tar -xzf p2rank_2.4.tar.gz && \
    rm p2rank_2.4.tar.gz
ENV PATH="/opt/p2rank_2.4:${PATH}"

WORKDIR /app

# Install CPU-only PyTorch first (avoids ~1.5 GB of unused CUDA wheels), then
# the remaining requirements (torch-geometric resolves against the CPU torch).
RUN pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Expose port
EXPOSE 8002

# Run application
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8002"]
