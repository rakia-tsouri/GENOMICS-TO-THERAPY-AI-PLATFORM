FROM python:3.10-slim

WORKDIR /app

# System dependencies:
#  - build-essential: compiling any C extensions during pip install
#  - libgl1 + libglib2.0-0: required by opencv (cv2) at import time even for the
#    headless build, which still loads some shared graphics libraries.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libgl1 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Install Python requirements first to leverage Docker layer caching.
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code (build context is this service directory).
COPY . .

# Ensure the shared volume directory for Grad-CAM overlays exists.
RUN mkdir -p /data/structures/gradcam

EXPOSE 8003

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8003"]
