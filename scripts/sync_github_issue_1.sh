#!/usr/bin/env bash
# Sync GitHub Issue #1 body from reports/github_issue_1_body.md
# Requires: gh CLI — https://cli.github.com/  then: gh auth login
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BODY="$ROOT/reports/github_issue_1_body.md"

if ! command -v gh >/dev/null 2>&1; then
  echo "Error: gh not found. Install: brew install gh && gh auth login"
  exit 1
fi

if [[ ! -f "$BODY" ]]; then
  echo "Error: missing $BODY"
  exit 1
fi

gh issue edit 1 --repo DmitriiZakharenko/tcp-modeling-gbm --body-file "$BODY"
echo "Updated https://github.com/DmitriiZakharenko/tcp-modeling-gbm/issues/1"
