FROM python:3.11-slim

WORKDIR /app

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy shared models and application code
COPY shared/ /app/shared/
COPY genomics_validation/api/ /app/api/

ENV PYTHONPATH=/app
ENV REDIS_HOST=redis
ENV REDIS_PORT=6379

EXPOSE 8000

CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
