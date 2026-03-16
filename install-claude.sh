#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/install-common.sh"

TOOL_NAME="Claude Code"
DEFAULT_TARGET_DIR="$HOME/.claude/skills"
PROJECT_SUBDIR=".claude/skills"

print_usage() {
  echo "Usage: $0 [OPTIONS]"
  echo ""
  echo "Install BreedingSkills to Claude Code"
  echo ""
  echo "Options:"
  print_common_options
}

post_install_message() {
  echo "Skills are now available to Claude Code."
  echo "Start a new Claude Code session if they do not appear immediately."
}

run_installer "$@"
