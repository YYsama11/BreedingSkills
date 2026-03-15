#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
skills_dir="${repo_root}/skills"
skills_home="${LOCAL_SKILLS_HOME:-${HOME}/.local/skills}"
target_root="${skills_home}"
mode="copy"

usage() {
  cat <<'EOF'
Usage:
  bash scripts/install_to_local_skills.sh [--symlink] [skill-name ...]

Behavior:
  - With no skill-name arguments, installs all skills from this repository
  - By default, copies each skill into $LOCAL_SKILLS_HOME
  - With --symlink, creates symbolic links instead of copying
EOF
}

skills=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --symlink)
      mode="symlink"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      skills+=("$1")
      shift
      ;;
  esac
done

mkdir -p "${target_root}"

if [[ ${#skills[@]} -eq 0 ]]; then
  while IFS= read -r -d '' skill_dir_path; do
    skills+=("$(basename "${skill_dir_path}")")
  done < <(find "${skills_dir}" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)
fi

for skill_name in "${skills[@]}"; do
  src="${skills_dir}/${skill_name}"
  dst="${target_root}/${skill_name}"
  if [[ ! -f "${src}/SKILL.md" ]]; then
    echo "Skipping ${skill_name}: SKILL.md not found" >&2
    continue
  fi
  rm -rf "${dst}"
  if [[ "${mode}" == "symlink" ]]; then
    ln -s "${src}" "${dst}"
    echo "Linked ${skill_name} -> ${dst}"
  else
    cp -a "${src}" "${dst}"
    echo "Copied ${skill_name} -> ${dst}"
  fi
done
