#!/usr/bin/env bash

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

INSTALL_MODE="global"
PROJECT_PATH=""
VALIDATE_ONLY=false
UPDATE_MODE=false
UNINSTALL_MODE=false
VERBOSE=false
DRY_RUN=false
FORCE=false
LIST_ONLY=false
SKILL_FILTER=""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_SKILLS_DIR="$SCRIPT_DIR/skills"
MARKER_FILE=".breedingskills-origin"

print_common_options() {
  echo "  --global              Install to the default global location (default)"
  echo "  --project [PATH]      Install to a project/workspace directory"
  echo "  --skills LIST         Install only specified skills (comma-separated)"
  echo "  --categories LIST     Alias of --skills for bioSkills-like usage"
  echo "  --list                List available skills"
  echo "  --validate            Validate skills before installing"
  echo "  --update              Only update installed skills whose source changed"
  echo "  --uninstall           Remove BreedingSkills installed by this installer"
  echo "  --dry-run             Preview installation without copying files"
  echo "  --verbose             Show detailed output"
  echo "  --force               Overwrite existing non-BreedingSkills directories"
  echo "  --help                Show this help message"
}

available_skill_dirs() {
  find "$REPO_SKILLS_DIR" -mindepth 1 -maxdepth 1 -type d | sort
}

normalize_csv() {
  echo "$1" | tr -d ' ' | tr ',' '\n' | sed '/^$/d'
}

skill_is_selected() {
  local skill_name="$1"
  if [[ -z "$SKILL_FILTER" ]]; then
    return 0
  fi
  while IFS= read -r selected; do
    [[ "$selected" == "$skill_name" ]] && return 0
  done < <(normalize_csv "$SKILL_FILTER")
  return 1
}

list_skills() {
  echo "Available BreedingSkills:"
  echo ""
  while IFS= read -r skill_dir; do
    local skill_name
    skill_name="$(basename "$skill_dir")"
    local description
    description="$(grep '^description:' "$skill_dir/SKILL.md" | head -1 | sed 's/^description:[[:space:]]*//')"
    echo "  $skill_name"
    [[ -n "$description" ]] && echo "    ${description:0:100}"
    echo ""
  done < <(available_skill_dirs)
}

validate_skill() {
  local skill_dir="$1"
  local skill_name
  skill_name="$(basename "$skill_dir")"
  local errors=()

  [[ -f "$skill_dir/SKILL.md" ]] || errors+=("Missing SKILL.md")
  if [[ -f "$skill_dir/SKILL.md" ]]; then
    grep -q '^---$' "$skill_dir/SKILL.md" || errors+=("Missing YAML frontmatter")
    grep -q '^name:' "$skill_dir/SKILL.md" || errors+=("Missing name field")
    grep -q '^description:' "$skill_dir/SKILL.md" || errors+=("Missing description field")
  fi

  for optional_dir in scripts references assets; do
    if [[ -d "$skill_dir/$optional_dir" ]] && [[ -z "$(find "$skill_dir/$optional_dir" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
      errors+=("Empty optional directory: $optional_dir")
    fi
  done

  if (( ${#errors[@]} > 0 )); then
    echo -e "  ${RED}FAIL${NC} $skill_name"
    for error in "${errors[@]}"; do
      echo "       - $error"
    done
    return 1
  fi

  [[ "$VERBOSE" == true ]] && echo -e "  ${GREEN}PASS${NC} $skill_name"
  return 0
}

validate_all_skills() {
  echo "Validating BreedingSkills..."
  echo ""
  local passed=0
  local failed=0

  while IFS= read -r skill_dir; do
    local skill_name
    skill_name="$(basename "$skill_dir")"
    skill_is_selected "$skill_name" || continue
    if validate_skill "$skill_dir"; then
      passed=$((passed + 1))
    else
      failed=$((failed + 1))
    fi
  done < <(available_skill_dirs)

  echo ""
  echo "Validation complete: $passed passed, $failed failed"
  [[ $failed -eq 0 ]]
}

resolve_target_dir() {
  if [[ "$INSTALL_MODE" == "global" ]]; then
    echo "$DEFAULT_TARGET_DIR"
    return
  fi

  local project_root="${PROJECT_PATH:-$PWD}"
  echo "$project_root/$PROJECT_SUBDIR"
}

prepare_target_dir() {
  local target_dir="$1"
  if [[ "$DRY_RUN" == true ]]; then
    return 0
  fi
  mkdir -p "$target_dir"
  if [[ ! -w "$target_dir" ]]; then
    echo -e "${RED}Error: Cannot write to target directory: $target_dir${NC}"
    exit 1
  fi
}

source_mtime() {
  local dir="$1"
  find "$dir" -type f ! -name '*.pyc' -printf '%T@\n' | sort -n | tail -1
}

target_mtime() {
  local dir="$1"
  if [[ ! -d "$dir" ]]; then
    echo "0"
    return
  fi
  find "$dir" -type f -printf '%T@\n' | sort -n | tail -1
}

copy_skill_tree() {
  local src_dir="$1"
  local target_dir="$2"

  rm -rf "$target_dir"
  mkdir -p "$target_dir"

  rsync -a \
    --exclude='__pycache__' \
    --exclude='*.pyc' \
    --exclude='.DS_Store' \
    "$src_dir/" "$target_dir/"

  cat > "$target_dir/$MARKER_FILE" <<EOF
repo=BreedingSkills
source=$src_dir
installed_at=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
tool=${TOOL_NAME}
EOF
}

install_skills() {
  local target_dir="$1"
  local installed=0
  local skipped=0

  [[ "$DRY_RUN" == true ]] && echo "Dry run - would install to: $target_dir" || echo "Installing BreedingSkills to: $target_dir"
  echo ""

  prepare_target_dir "$target_dir"

  while IFS= read -r skill_dir; do
    local skill_name dest_dir
    skill_name="$(basename "$skill_dir")"
    skill_is_selected "$skill_name" || continue
    dest_dir="$target_dir/$skill_name"

    if [[ -d "$dest_dir" ]] && [[ ! -f "$dest_dir/$MARKER_FILE" ]] && [[ "$FORCE" != true ]]; then
      echo -e "  ${YELLOW}Skip existing non-BreedingSkills directory:${NC} $skill_name"
      echo "    Use --force to overwrite."
      skipped=$((skipped + 1))
      continue
    fi

    if [[ "$UPDATE_MODE" == true ]] && [[ -d "$dest_dir" ]]; then
      local src_ts dst_ts
      src_ts="$(source_mtime "$skill_dir")"
      dst_ts="$(target_mtime "$dest_dir")"
      if awk "BEGIN {exit !($src_ts <= $dst_ts)}"; then
        [[ "$VERBOSE" == true ]] && echo "  Skipped (unchanged): $skill_name"
        skipped=$((skipped + 1))
        continue
      fi
    fi

    if [[ "$DRY_RUN" == true ]]; then
      echo "  Would install: $skill_name"
      installed=$((installed + 1))
      continue
    fi

    copy_skill_tree "$skill_dir" "$dest_dir"
    echo "  Installed: $skill_name"
    installed=$((installed + 1))
  done < <(available_skill_dirs)

  echo ""
  if [[ "$DRY_RUN" == true ]]; then
    echo "Would install: $installed skills"
  else
    echo "Installation complete."
    echo "  Installed: $installed"
    if [[ "$UPDATE_MODE" == true ]]; then
      echo "  Skipped (unchanged): $skipped"
    fi
  fi

  return 0
}

uninstall_skills() {
  local target_dir="$1"
  if [[ ! -d "$target_dir" ]]; then
    echo "No skills directory found at: $target_dir"
    return 0
  fi

  echo "Removing BreedingSkills from: $target_dir"
  echo ""

  local removed=0
  while IFS= read -r dir; do
    [[ -d "$dir" ]] || continue
    rm -rf "$dir"
    echo "  Removed: $(basename "$dir")"
    removed=$((removed + 1))
  done < <(find "$target_dir" -mindepth 1 -maxdepth 1 -type d -exec test -f "{}/$MARKER_FILE" \; -print | sort)

  echo ""
  echo "Uninstall complete. Removed $removed skills."
}

run_installer() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --global)
        INSTALL_MODE="global"
        shift
        ;;
      --project)
        INSTALL_MODE="project"
        if [[ $# -ge 2 && ! "$2" =~ ^-- ]]; then
          PROJECT_PATH="$2"
          shift
        fi
        shift
        ;;
      --skills|--categories)
        if [[ $# -lt 2 || "$2" =~ ^-- ]]; then
          echo "Error: $1 requires a comma-separated list"
          exit 1
        fi
        SKILL_FILTER="$2"
        shift 2
        ;;
      --list)
        LIST_ONLY=true
        shift
        ;;
      --validate)
        VALIDATE_ONLY=true
        shift
        ;;
      --update)
        UPDATE_MODE=true
        shift
        ;;
      --uninstall)
        UNINSTALL_MODE=true
        shift
        ;;
      --dry-run)
        DRY_RUN=true
        shift
        ;;
      --verbose)
        VERBOSE=true
        shift
        ;;
      --force)
        FORCE=true
        shift
        ;;
      --help|-h)
        print_usage
        exit 0
        ;;
      *)
        if type parse_extra_arg &>/dev/null && parse_extra_arg "$1"; then
          shift "${EXTRA_SHIFT:-1}"
        else
          echo "Unknown option: $1"
          echo ""
          print_usage
          exit 1
        fi
        ;;
    esac
  done

  if [[ "$LIST_ONLY" == true ]]; then
    list_skills
    exit 0
  fi

  if [[ "$VALIDATE_ONLY" == true ]]; then
    validate_all_skills
    exit $?
  fi

  local target_dir
  target_dir="$(resolve_target_dir)"

  if [[ "$UNINSTALL_MODE" == true ]]; then
    uninstall_skills "$target_dir"
    exit 0
  fi

  validate_all_skills
  install_skills "$target_dir"

  if type post_install_message &>/dev/null; then
    echo ""
    post_install_message
  fi
}
