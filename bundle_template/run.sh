#!/usr/bin/env bash
# =============================================================================
#  GOBSMACKED run bundle
#
#      ./run.sh              run every stage that has not finished
#      ./run.sh --list       show what would run, and what is already done
#      ./run.sh --stage dock rerun from a stage onward
#
#  One step, and no prerequisites beyond curl and bash: this installs pixi if
#  the machine has not got it, builds the environment from the lock file
#  shipped beside this script, and runs the pipeline.
#
#  Re-running resumes: finished stages are skipped, so a failure costs only the
#  stage that failed.
# =============================================================================
set -euo pipefail

PHOS=$'\033[38;2;93;225;230m'
GREEN=$'\033[38;2;126;226;168m'
AMBER=$'\033[38;2;255;180;84m'
RED=$'\033[38;2;255;92;92m'
BOLD=$'\033[1m'
RESET=$'\033[0m'

step() { printf "%s\n" "${PHOS}→${RESET} $1"; }
ok()   { printf "%s\n" "${GREEN}✓${RESET} $1"; }
warn() { printf "%s\n" "${AMBER}⚠${RESET} $1"; }
die()  { printf "%s\n" "${RED}✗${RESET} $1" >&2; exit 1; }

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

JOB=$(awk '/^job_id:/ { print $2; exit }' campaign.yaml 2>/dev/null || echo "unknown")
printf "\n%s\n" "${BOLD}${PHOS}GOBSMACKED${RESET} ${BOLD}${JOB}${RESET}"
printf "%s\n\n" "fold, dock, relax, summarise: on this machine"

# ---------------------------------------------------------------------------
#  Checks worth making before a job that takes an hour
# ---------------------------------------------------------------------------
case "$(uname -s)" in
  Darwin|Linux) ;;
  *) warn "unrecognised platform '$(uname -s)'; the environment is locked for linux-64 and osx-arm64." ;;
esac

# The environment is about 2.5 GB and the results a few hundred megabytes more.
if command -v df >/dev/null 2>&1; then
  FREE_GB=$(df -Pk "$HERE" | awk 'NR==2 { printf "%d", $4 / 1048576 }')
  if [ -n "${FREE_GB:-}" ] && [ "$FREE_GB" -lt 5 ]; then
    warn "only ${FREE_GB} GB free here. The environment alone is about 2.5 GB."
  fi
fi

# Writing a trajectory into an iCloud-synced folder is a real problem: the sync
# daemon competes for I/O and "Optimise Mac Storage" can evict files mid-run.
case "$HERE" in
  "$HOME"/Documents/*|"$HOME"/Desktop/*|*"/Library/Mobile Documents/"*)
    if [ -d "$HOME/Library/Mobile Documents/com~apple~CloudDocs" ]; then
      warn "this bundle is in a folder macOS may sync to iCloud."
      warn "iCloud can evict files mid-run and its daemon competes for disk I/O."
    fi
    ;;
esac

# ---------------------------------------------------------------------------
#  pixi
# ---------------------------------------------------------------------------
if ! command -v pixi >/dev/null 2>&1; then
  if [ -x "$HOME/.pixi/bin/pixi" ]; then
    export PATH="$HOME/.pixi/bin:$PATH"
  else
    step "installing pixi (pixi.sh), once, into ~/.pixi"
    curl -fsSL https://pixi.sh/install.sh | bash >/dev/null \
      || die "could not install pixi. Install it yourself (https://pixi.sh) and run this again."
    export PATH="$HOME/.pixi/bin:$PATH"
    command -v pixi >/dev/null 2>&1 \
      || die "pixi installed but is not on PATH. Open a new terminal and run ./run.sh again."
  fi
fi
ok "pixi $(pixi --version 2>/dev/null | awk '{print $2}')"

# ---------------------------------------------------------------------------
#  Environment, from the lock file shipped with this bundle
# ---------------------------------------------------------------------------
if [ ! -d ".pixi/envs/default" ]; then
  step "fetching the environment: about 2.5 GB, once, then it is reused"
  printf "  %s\n" "no solve is needed, the versions are pinned in pixi.lock"
fi
# --locked refuses to run if pixi.lock and pixi.toml disagree, which is the
# point of shipping the lock: the environment is the one this was tested with,
# not whatever solves today.
pixi install --locked \
  || die "the environment could not be installed. Check the network and try again."
ok "environment ready"

# The affinity stage shells into its own environment, which pixi would otherwise
# install the moment that stage runs: a multi-gigabyte download landing an hour
# into a job, with a progress bar for molecular dynamics on screen. Installed
# here instead, where a long download is what the reader is already waiting for.
# Read from inside the affinity block rather than by grepping the whole file:
# `include:` is unique in campaign.yaml today, and a bare grep would silently
# start matching the day it is not.
WANTS_AFFINITY=$(awk '
  /^affinity:/ { inblock = 1; next }
  /^[a-z_]+:/  { inblock = 0 }
  inblock && /include:[[:space:]]*true/ { print "yes"; exit }
' campaign.yaml 2>/dev/null || true)
if [ "${WANTS_AFFINITY:-}" = "yes" ]; then
  if [ ! -d ".pixi/envs/affinity" ]; then
    step "fetching the affinity environment (Boltz-2 and torch), about 3 GB, once"
  fi
  if ! pixi install --locked -e affinity; then
    warn "the affinity environment could not be installed."
    warn "the run continues; the affinity stage will record why it was skipped."
  else
    ok "affinity environment ready"
  fi
fi

# ---------------------------------------------------------------------------
#  Run
# ---------------------------------------------------------------------------
printf "\n"
pixi run --locked gobsmacked "$@"
