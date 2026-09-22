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
#   (no args)    ./clean_icon_files.sh         one sweep of paper/, by hand
#   after reboot ./clean_icon_files.sh repo    one sweep of the WHOLE repo,
#                                              including .git
#
# The sweeper is controlled by a lock file holding a token, not by a PID: it
# exits within one interval of the token changing or the lock disappearing, and
# in any case at MAX_SECONDS. Nothing can be left behind to kill, and no stale
# PID can ever be signalled.
#
# WHY `repo` EXISTS. The stubs are harmless in the working tree (gitignored) but
# not inside .git: on 2026-07-27 a Drive re-sync left 272 of them there, one in
# refs/, and git read it as a ref -- `git fetch` died with "fatal: bad object
# refs/Icon?". A reboot re-syncs the whole mount, which is exactly when that
# happens, so sweep the repo before the first git command after a restart.
#
# SCOPE. Two conditions, both required: the exact four-character name Icon + CR
# (so it cannot touch a real file called Icons or Icon.png) AND zero bytes. No
# legitimate git file can be both -- a ref holds a 40-character SHA and objects
# are zlib-compressed -- so the `repo` sweep cannot damage the repository. Every
# match is a Finder icon stub. This is the only thing the script deletes.

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

# -size 0 is a safety condition, not an optimisation: it is what makes the
# `repo` sweep safe to point at .git. See SCOPE in the header.
sweep() {
  find "${1:-$ROOT}" -type f -name "$ICON_NAME" -size 0 -delete 2>/dev/null || true
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
  repo)
    top="$(git -C "$ROOT" rev-parse --show-toplevel 2>/dev/null || true)"
    if [ -z "$top" ]; then
      echo "$(basename "$0"): not inside a git repository ($ROOT)" >&2
      exit 2
    fi
    n_tree=$(find "$top" -path "$top/.git" -prune -o \
                  -type f -name "$ICON_NAME" -size 0 -print 2>/dev/null | wc -l | tr -d ' ')
    n_git=$(find "$top/.git" -type f -name "$ICON_NAME" -size 0 2>/dev/null | wc -l | tr -d ' ')
    sweep "$top"
    echo "swept $n_tree stub(s) from the working tree and $n_git from .git"
    ;;
  *)
    echo "usage: $(basename "$0") [watch|stop|once|repo]" >&2
    exit 2
    ;;
esac

exit 0
