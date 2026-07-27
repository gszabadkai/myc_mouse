#!/usr/bin/env bash
#
# Delete the macOS custom-folder-icon stubs -- files named exactly "Icon" plus a
# carriage return -- from the paper project tree, and keep deleting them for as
# long as a render is in flight.
#
# WHY THIS EXISTS
#
# This repo lives on a Google Drive for Desktop mount, and Drive stamps a
# zero-byte "Icon\r" into every folder it syncs. Quarto walks the supporting
# directory it creates for the render (paper/myc_mito_files/) and calls readdir
# on each entry; the Icon stub carries a resource fork, so readdir fails and the
# whole render dies at the very last step, after every chunk has run:
#
#   ERROR: NotADirectory (os error 20): readdir '.../paper/myc_mito_files/Icon'
#
# A pre-render clean is NOT sufficient, and that is the point of the sweeper.
# Measured 2026-07-27: Quarto creates myc_mito_files/ mid-render and Drive
# stamps it within seconds, so the file that kills the render is written inside
# the render window -- after pre-render has finished and before post-render
# begins. One render succeeded and the next failed on identical inputs; it is a
# race, not a state.
#
# HOW
#
#   pre-render   ./clean_icon_files.sh watch   sweep, then leave a sweeper
#                                              running for the render window
#   post-render  ./clean_icon_files.sh stop    retire the sweeper, sweep again
#   (no args)    ./clean_icon_files.sh         one sweep, for use by hand
#
# The sweeper is controlled by a lock file holding a token, not by a PID: it
# exits within one interval of the token changing or the lock disappearing, and
# in any case at MAX_SECONDS. Nothing can be left behind to kill, and no stale
# PID can ever be signalled.
#
# SCOPE. The match is the exact four-character name Icon + CR, so it cannot
# touch a real file called Icons or Icon.png. Every match is a Finder icon stub;
# none is project content. This is the only thing the script deletes.

set -euo pipefail

readonly ICON_NAME=$'Icon\r'
readonly INTERVAL=0.5
readonly MAX_SECONDS=1800

# QUARTO_PROJECT_ROOT is set for post-render but NOT for pre-render (verified
# 2026-07-27); in both cases the cwd is the project directory, so derive it from
# the script's own location and use that as the fallback.
ROOT="${QUARTO_PROJECT_ROOT:-$(cd "$(dirname "$0")" && pwd)}"
readonly ROOT
readonly LOCK="${TMPDIR:-/tmp}/myc_mito_icon_sweeper.lock"

sweep() {
  find "$ROOT" -type f -name "$ICON_NAME" -delete 2>/dev/null || true
}

start_sweeper() {
  local token
  token="$$-$(date +%s)"
  printf '%s' "$token" > "$LOCK"

  (
    local deadline
    deadline=$(( $(date +%s) + MAX_SECONDS ))
    while [ "$(cat "$LOCK" 2>/dev/null || true)" = "$token" ] &&
          [ "$(date +%s)" -lt "$deadline" ]; do
      sweep
      sleep "$INTERVAL"
    done
  ) >/dev/null 2>&1 &
}

case "${1:-once}" in
  watch)
    rm -f "$LOCK"          # retire any sweeper left over from a crashed render
    sweep
    start_sweeper
    ;;
  stop)
    rm -f "$LOCK"          # the sweeper notices within one interval and exits
    sweep
    ;;
  once)
    sweep
    ;;
  *)
    echo "usage: $(basename "$0") [watch|stop]" >&2
    exit 2
    ;;
esac

exit 0
