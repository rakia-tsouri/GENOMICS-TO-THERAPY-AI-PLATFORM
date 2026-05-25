import redis
import json
import os
from datetime import timedelta

class RedisCache:
    def __init__(self):
        redis_host = os.getenv("REDIS_HOST", "redis")
        redis_port = int(os.getenv("REDIS_PORT", 6379))
        self.client = redis.Redis(host=redis_host, port=redis_port, decode_responses=True)

    def get(self, key: str):
        data = self.client.get(key)
        return json.loads(data) if data else None

    def set(self, key: str, value: dict, expire_hours: int = 24):
        self.client.setex(key, timedelta(hours=expire_hours), json.dumps(value))

cache = RedisCache()
