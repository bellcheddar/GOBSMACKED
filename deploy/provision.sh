#!/usr/bin/env bash
# One-time droplet provisioning for GOBSMACKED. Run as root ON the droplet,
# AFTER the code is at /opt/gobsmacked (push it with deploy/deploy.sh first).
#
#   sudo SERVER_NAME=gobsmacked.mdeller.com bash /opt/gobsmacked/deploy/provision.sh
#
# Idempotent: safe to re-run. Installs system packages, a service user, the
# Python environment, the systemd units, the nginx site and a certificate.
set -euo pipefail

APP_DIR=/opt/gobsmacked
APP_USER=gobsmacked
BIND_ADDR="${BIND_ADDR:-127.0.0.1:8009}"

if [[ -f "$APP_DIR/.env" ]]; then
  set -a; # shellcheck disable=SC1091
  source "$APP_DIR/.env"; set +a
fi
SERVER_NAME="${SERVER_NAME:-gobsmacked.mdeller.com}"

echo "==> GOBSMACKED provisioning for ${SERVER_NAME}"
[[ $EUID -eq 0 ]] || { echo "Run as root (sudo)."; exit 1; }
[[ -f "$APP_DIR/wsgi.py" ]] || { echo "No code at $APP_DIR. Push it first: bash deploy/deploy.sh"; exit 1; }

echo "==> Installing system packages"
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
# No libopenbabel-dev and no swig: PLIP's openbabel dependency ships manylinux
# wheels for every Python this box would run (checked on PyPI, cp310 to cp312),
# so pip never builds it from source. Installing the headers anyway would add
# a few hundred megabytes to a disk already at 78 %.
# ffmpeg encodes the trajectory clip on the results page: libx264 is the only
# encoder asked for, and the analysis degrades to a note in the panel rather
# than an error if the binary is missing. Measured at about 1 GB installed on
# a box that was at 79 %.
apt-get install -y -qq python3-venv python3-pip python3-dev build-essential \
  nginx certbot python3-certbot-nginx rsync ffmpeg

echo "==> Creating service user '${APP_USER}'"
id -u "$APP_USER" &>/dev/null || useradd --system --home "$APP_DIR" --shell /usr/sbin/nologin "$APP_USER"
mkdir -p "$APP_DIR/data/runs" "$APP_DIR/data/structures"
chown -R "$APP_USER:$APP_USER" "$APP_DIR"

echo "==> Building the Python environment"
if [[ ! -x "$APP_DIR/.venv/bin/python" ]]; then
  sudo -u "$APP_USER" python3 -m venv "$APP_DIR/.venv"
fi
sudo -u "$APP_USER" "$APP_DIR/.venv/bin/pip" install --quiet --upgrade pip
sudo -u "$APP_USER" "$APP_DIR/.venv/bin/pip" install --quiet -r "$APP_DIR/requirements.txt"

echo "==> Initialising the database"
sudo -u "$APP_USER" "$APP_DIR/.venv/bin/python" -c \
  "import sys; sys.path.insert(0, '$APP_DIR'); from app import config, db; config.ensure_dirs(); db.init_db(); print('ready')"

echo "==> Installing systemd units"
cp "$APP_DIR/deploy/gobsmacked-web.service"   /etc/systemd/system/
cp "$APP_DIR/deploy/gobsmacked-prune.service" /etc/systemd/system/
cp "$APP_DIR/deploy/gobsmacked-prune.timer"   /etc/systemd/system/
systemctl daemon-reload
systemctl enable --now gobsmacked-web.service
systemctl enable --now gobsmacked-prune.timer

echo "==> Installing the nginx site"
sed -e "s|__SERVER_NAME__|${SERVER_NAME}|g" -e "s|__BIND_ADDR__|${BIND_ADDR}|g" \
  "$APP_DIR/deploy/nginx-gobsmacked.conf" > /etc/nginx/sites-available/gobsmacked
ln -sf /etc/nginx/sites-available/gobsmacked /etc/nginx/sites-enabled/gobsmacked
nginx -t && systemctl reload nginx

echo "==> Requesting the TLS certificate"
# Always run certbot, never skip on "a certificate already exists": the site
# file above is re-templated from source on every run with plain HTTP only, so
# a skip would leave a freshly written vhost with no SSL block at all. nginx
# then serves some other vhost's certificate on 443 for this hostname, which is
# a real incident a sibling app on this droplet already had. certbot --nginx
# reuses a still-valid certificate rather than re-requesting one, so running it
# unconditionally is both safe and correct.
#
# Retried because the droplet-wide certbot renewal timer can hold certbot's
# lock at the same moment ("Another instance of Certbot is already running").
certbot_ok=0
for attempt in 1 2 3; do
  if certbot --nginx -d "$SERVER_NAME" --non-interactive --agree-tos \
       -m "${CERTBOT_EMAIL:-marc@marcdeller.com}" --redirect; then
    certbot_ok=1; break
  fi
  echo "    certbot attempt ${attempt}/3 failed, retrying in 10s..."
  sleep 10
done
if [[ "$certbot_ok" -ne 1 ]]; then
  echo "    certbot failed after 3 attempts. Is DNS for ${SERVER_NAME} pointed here yet?"
  echo "    ${SERVER_NAME} has NO SSL config right now. Re-run: certbot --nginx -d ${SERVER_NAME}"
  exit 1
fi

# certbot's `listen 443 ssl;` does not enable HTTP/2 on nginx 1.24, so add it
# here, idempotently, which also fixes an already-provisioned site on a re-run.
if grep -q "listen.*443 ssl" /etc/nginx/sites-available/gobsmacked && \
   ! grep -q "listen.*443 ssl http2" /etc/nginx/sites-available/gobsmacked; then
  echo "==> Enabling HTTP/2"
  python3 - <<'PYEOF'
import re
path = "/etc/nginx/sites-available/gobsmacked"
text = open(path).read()
text = re.sub(r'listen ((?:\[::\]:)?443) ssl( ipv6only=on)?;',
              lambda m: f'listen {m.group(1)} ssl http2{m.group(2) or ""};', text)
open(path, "w").write(text)
PYEOF
  nginx -t && systemctl reload nginx
fi

echo "==> Done."
systemctl --no-pager --lines=3 status gobsmacked-web || true
echo "    Site:  https://${SERVER_NAME}/"
echo "    Timer: $(systemctl list-timers gobsmacked-prune.timer --no-pager | sed -n 2p)"
