#!/usr/bin/env bash
# Capture the 6 platform screenshots used in the Results chapter of the report.
# Prereqs:
#   - docker compose stack running (frontend on :3000, gateway on :8080)
#   - public/auth_bridge.html present in the running frontend container
#     (in this repo it's at frontend/public/auth_bridge.html ; it ships into
#      the image at build time)
#   - Google Chrome installed at the standard macOS location
#
# Usage:  bash docs/report/take_screenshots.sh
set -eu

CHROME='/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
GW=http://localhost:8080/api/v1
WEB=http://localhost:3000
SHOTS="$(cd "$(dirname "$0")" && pwd)/images/screenshots"
mkdir -p "$SHOTS"
rm -f "$SHOTS"/*.png

# 0. wait for gateway + frontend to be reachable (cold-start can take 30-60s
#    while the gateway waits for postgres and runs its seed)
echo "Waiting for services..."
for i in $(seq 1 40); do
  ok_gw=$(curl -fs -m 3 http://localhost:8080/health 2>/dev/null | grep -c healthy || true)
  ok_fe=$(curl -fs -m 3 -o /dev/null -w '%{http_code}' http://localhost:3000/login)
  [ "$ok_gw" = "1" ] && [ "$ok_fe" = "200" ] && break
  sleep 2
done
echo "  gateway /health: ok=$ok_gw    frontend /login: HTTP=$ok_fe"

# 1. login as the seeded demo researcher
LOGIN=$(curl -fs -X POST "$GW/auth/login" -H 'Content-Type: application/json' \
  -d '{"email":"researcher@medconnect.dev","password":"research12345"}')
if [ -z "$LOGIN" ]; then
  echo "ERROR: gateway did not return a login response. Is 'docker compose up' done seeding?"
  echo "  try: docker logs gtt-gateway | tail -20"
  exit 1
fi
TOK=$(echo "$LOGIN" | python3 -c "import sys,json;print(json.load(sys.stdin)['access_token'])")

# 2. pick a seeded job/report (first one, which is the oldest = TP53 reference)
JID=$(curl -fs "$GW/jobs" -H "Authorization: Bearer $TOK" \
  | python3 -c "import sys,json;j=json.load(sys.stdin);print(j[-1]['id'])")
RID=$(curl -fs "$GW/reports" -H "Authorization: Bearer $TOK" \
  | python3 -c "import sys,json;r=json.load(sys.stdin);print(r[-1]['id'])")
echo "Using seeded job_id=$JID  report_id=$RID"

enc() { python3 -c "import urllib.parse,sys;print(urllib.parse.quote(sys.argv[1]))" "$1"; }

shot() {
  local name=$1 path=$2 budget=${3:-15000}
  local url="$WEB/auth_bridge.html?token=$TOK&to=$(enc "$path")"
  "$CHROME" --headless=new --hide-scrollbars --disable-gpu --no-sandbox \
    --window-size=1440,1100 --virtual-time-budget="$budget" \
    --screenshot="$SHOTS/$name.png" "$url" 2>/dev/null
  printf "  %-22s %s bytes\n" "$name.png" "$(stat -f%z "$SHOTS/$name.png" 2>/dev/null)"
}

# 3. login screenshot — no auth needed
"$CHROME" --headless=new --hide-scrollbars --disable-gpu --no-sandbox \
  --window-size=1440,900 --virtual-time-budget=6000 \
  --screenshot="$SHOTS/01_login.png" "$WEB/login" 2>/dev/null
printf "  %-22s %s bytes\n" "01_login.png" "$(stat -f%z "$SHOTS/01_login.png")"

shot 02_dashboard      /dashboard
shot 03_new_analysis   /jobs/new
# generous budget so 3Dmol.js fetches the PDB and renders
shot 04_results        "/jobs/$JID" 22000
shot 05_fusion         "/jobs/$JID" 22000
shot 06_report         "/reports/$RID" 12000

echo "Done. Screenshots saved to $SHOTS"
