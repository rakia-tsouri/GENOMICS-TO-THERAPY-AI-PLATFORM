# node:20-slim (Debian glibc) avoids the Next.js SWC linux/arm64 musl issue
# seen on node:20-alpine.
FROM node:20-slim

WORKDIR /app

# Install dependencies
COPY package.json package-lock.json* ./
RUN npm ci

# Copy source and build
COPY . .
RUN npm run build

EXPOSE 3000

CMD ["npm", "start"]
