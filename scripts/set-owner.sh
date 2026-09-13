#!/usr/bin/env bash
# Fill in the GitHub owner across the repo, once, after you know the account name.
#   ./scripts/set-owner.sh your-github-username
set -euo pipefail
[ $# -eq 1 ] || { echo "usage: $0 <github-username>" >&2; exit 2; }
OWNER="$1"
FILES=(README.md pyproject.toml CITATION.cff CONTRIBUTING.md)
for f in "${FILES[@]}"; do
  [ -f "$f" ] && perl -pi -e "s/\bOWNER\b/$OWNER/g" "$f"
done
echo "Set owner to '$OWNER' in: ${FILES[*]}"
grep -rn "github.com/$OWNER/admissible" "${FILES[@]}" | head
