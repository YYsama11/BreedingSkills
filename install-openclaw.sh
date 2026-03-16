#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/install-common.sh"

TOOL_NAME="OpenClaw"
DEFAULT_TARGET_DIR="$HOME/.openclaw/skills"
PROJECT_SUBDIR="skills"

print_usage() {
  echo "Usage: $0 [OPTIONS]"
  echo ""
  echo "Install BreedingSkills to OpenClaw"
  echo ""
  echo "Options:"
  print_common_options
}

post_install_message() {
  echo "Skills are now available to OpenClaw."
  echo "Start a new OpenClaw session if they do not appear immediately."
}

run_installer "$@"
