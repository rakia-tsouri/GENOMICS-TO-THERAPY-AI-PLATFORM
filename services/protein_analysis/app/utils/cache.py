import redis
import json
import os
from typing import Optional, Any
import logging

# Settings from environment variables
REDIS_HOST = os.getenv("REDIS_HOST", "redis")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
REDIS_DB = int(os.getenv("REDIS_DB", 0))

# Defined TTLs
TTL_24H = 24 * 60 * 60
TTL_7D = 7 * 24 * 60 * 60
TTL_30D = 30 * 24 * 60 * 60

logger = logging.getLogger(__name__)

class CacheClient:
    def __init__(self):
        try:
            self.client = redis.Redis(
                host=REDIS_HOST, 
                port=REDIS_PORT, 
                db=REDIS_DB, 
                decode_responses=True
            )
        except Exception as e:
            logger.error(f"Failed to connect to Redis: {e}")
            self.client = None

    def get(self, key: str) -> Optional[Any]:
        if not self.client:
            return None
        try:
            data = self.client.get(key)
            return json.loads(data) if data else None
        except Exception as e:
            logger.error(f"Redis GET error for {key}: {e}")
            return None

    def set(self, key: str, value: Any, ttl: int = TTL_24H):
        if not self.client:
            return
        try:
            self.client.setex(key, ttl, json.dumps(value))
        except Exception as e:
            logger.error(f"Redis SET error for {key}: {e}")

cache = CacheClient()
