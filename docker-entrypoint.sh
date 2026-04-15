#!/usr/bin/env bash
#
# Re-exec as a user matching the host-side mount owner so output files
# don't end up root-owned. Detects UID/GID from /data automatically;
# HOST_UID/HOST_GID env vars override the auto-detection.
#
set -euo pipefail

HOST_UID=${HOST_UID:-$(stat -c '%u' /data 2>/dev/null || echo 1000)}
HOST_GID=${HOST_GID:-$(stat -c '%g' /data 2>/dev/null || echo 1000)}

if [ "$(id -u)" = "0" ]; then
    groupadd -o -g "$HOST_GID" dnmb 2>/dev/null || true
    useradd  -o -u "$HOST_UID" -g "$HOST_GID" -M -d /data -s /bin/bash dnmb 2>/dev/null || true
    exec gosu dnmb:dnmb "$@"
fi

exec "$@"
