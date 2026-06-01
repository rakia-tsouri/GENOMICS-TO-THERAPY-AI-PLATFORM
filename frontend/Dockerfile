# node:20-slim (Debian glibc) avoids the Next.js SWC linux/arm64 musl issue
# seen on node:20-alpine.
FROM node:20-slim

WORKDIR /app

# Install dependencies
COPY package.json package-lock.json* ./
RUN npm ci

# Copy source and build.
# Next.js inlines NEXT_PUBLIC_* env vars at BUILD TIME, so the gateway URL
# must be present during `next build` — receive it from docker-compose's
# build arg and expose it to the build step.
ARG NEXT_PUBLIC_API_URL=http://localhost:8080
ENV NEXT_PUBLIC_API_URL=$NEXT_PUBLIC_API_URL
COPY . .
RUN npm run build

EXPOSE 3000

CMD ["npm", "start"]
