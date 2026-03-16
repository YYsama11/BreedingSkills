#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/install-common.sh"

TOOL_NAME="Codex CLI"
DEFAULT_TARGET_DIR="$HOME/.codex/skills"
PROJECT_SUBDIR=".codex/skills"

print_usage() {
  echo "Usage: $0 [OPTIONS]"
  echo ""
  echo "Install BreedingSkills to Codex CLI"
  echo ""
  echo "Options:"
  print_common_options
}

post_install_message() {
  echo "Skills are now available to Codex CLI."
  echo "Start a new Codex session if they do not appear immediately."
}

run_installer "$@"
