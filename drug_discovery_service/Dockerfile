FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    fpocket \
    wget \
    default-jre \
    && rm -rf /var/lib/apt/lists/*

# Install P2Rank
WORKDIR /opt
RUN wget https://github.com/rdk/p2rank/releases/download/v2.4/p2rank_2.4.tar.gz && \
    tar -xzf p2rank_2.4.tar.gz && \
    rm p2rank_2.4.tar.gz
ENV PATH="/opt/p2rank_2.4:${PATH}"

WORKDIR /app

# Install Python requirements
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Expose port
EXPOSE 8002

# Run application
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8002"]
