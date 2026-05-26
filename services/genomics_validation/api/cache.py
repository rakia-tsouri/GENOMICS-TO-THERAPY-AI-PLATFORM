import redis
import json
import os
from typing import Optional, Dict

# Settings from environment variables
REDIS_HOST = os.getenv("REDIS_HOST", "localhost")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
TTL_24H = 24 * 60 * 60

class CacheClient:
    def __init__(self):
        try:
            self.client = redis.Redis(host=REDIS_HOST, port=REDIS_PORT, db=0, decode_responses=True)
        except Exception:
            self.client = None

    def get_sequence(self, gene_id: str) -> Optional[str]:
        if not self.client:
            return None
        try:
            return self.client.get(f"gene:{gene_id}")
        except Exception:
            return None

    def set_sequence(self, gene_id: str, sequence: str):
        if not self.client:
            return
        try:
            self.client.setex(f"gene:{gene_id}", TTL_24H, sequence)
        except Exception:
            pass

cache = CacheClient()
