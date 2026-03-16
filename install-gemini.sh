#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/install-common.sh"

TOOL_NAME="Gemini CLI"
DEFAULT_TARGET_DIR="$HOME/.gemini/skills"
PROJECT_SUBDIR=".gemini/skills"

print_usage() {
  echo "Usage: $0 [OPTIONS]"
  echo ""
  echo "Install BreedingSkills to Gemini CLI"
  echo ""
  echo "Options:"
  print_common_options
}

post_install_message() {
  echo "Skills are now available to Gemini CLI."
  echo "Start a new Gemini session if they do not appear immediately."
}

run_installer "$@"
