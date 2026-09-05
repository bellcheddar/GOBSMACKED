#!/usr/bin/env bash
# Push GOBSMACKED from the Mac to the droplet and restart the service.
#
#   bash deploy/deploy.sh
#
# Reads DROPLET_SSH / DROPLET_PATH from .env. Excludes the venv, the runtime
# data, the database and secrets, so the server's own state is never clobbered.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if [[ -f .env ]]; then set -a; source .env; set +a; fi
DROPLET_SSH="${DROPLET_SSH:-}"
DROPLET_PATH="${DROPLET_PATH:-/opt/gobsmacked}"
SSH_KEY="${SSH_KEY:-}"

if [[ -z "$DROPLET_SSH" ]]; then
  echo "DROPLET_SSH is not set. Copy .env.example to .env and fill it in."; exit 1
fi

SSH_OPTS=()
[[ -n "$SSH_KEY" ]] && SSH_OPTS=(-e "ssh -i ${SSH_KEY/#\~/$HOME}")

echo "==> Syncing to ${DROPLET_SSH}:${DROPLET_PATH}"
# ${arr[@]+"${arr[@]}"} expands to nothing when empty without tripping `set -u`,
# which macOS's bash 3.2 needs.
rsync -az --delete ${SSH_OPTS[@]+"${SSH_OPTS[@]}"} \
  --exclude '.venv/' --exclude 'data/' --exclude '__pycache__/' \
  --exclude '*.pyc' --exclude '.git/' --exclude '.env' \
  --exclude 'gobsmacked.db' --exclude 'gobsmacked.db-*' \
  --exclude 'tests/fixtures/_structures/' \
  ./ "${DROPLET_SSH}:${DROPLET_PATH}/"

echo "==> Installing dependencies and restarting"
SSH_CMD=(ssh)
[[ -n "$SSH_KEY" ]] && SSH_CMD=(ssh -i "${SSH_KEY/#\~/$HOME}")
"${SSH_CMD[@]}" "$DROPLET_SSH" bash -s <<REMOTE
set -euo pipefail
cd "${DROPLET_PATH}"
if [[ ! -x .venv/bin/python ]]; then
  echo "No venv yet: run deploy/provision.sh as root first."; exit 0
fi
sudo -u gobsmacked env PIP_NO_CACHE_DIR=1 ./.venv/bin/pip install --quiet -r requirements.txt
# rsync runs as root and leaves new files root-owned. Chown them, but prune the
# venv, the runtime data and the live database: a blanket chown over an open
# SQLite file has taken a sibling app down mid-run.
sudo find "${DROPLET_PATH}" \
  -path "${DROPLET_PATH}/data" -prune -o \
  -path "${DROPLET_PATH}/.venv" -prune -o \
  -name 'gobsmacked.db*' -prune -o \
  -exec chown gobsmacked:gobsmacked {} +
sudo systemctl restart gobsmacked-web.service
sleep 2
sudo systemctl is-active gobsmacked-web.service
REMOTE

echo "==> Verifying the live site"
# The deploy script exiting 0 says the commands ran, not that the site serves
# the new build. Fetch a page back and check it.
SERVER_NAME="${SERVER_NAME:-gobsmacked.mdeller.com}"
if curl -sf --max-time 20 "https://${SERVER_NAME}/healthz" | grep -q '"ok"'; then
  echo "    https://${SERVER_NAME}/ is up: $(curl -s --max-time 20 https://${SERVER_NAME}/healthz)"
else
  echo "    WARNING: https://${SERVER_NAME}/healthz did not answer as expected."
  exit 1
fi
