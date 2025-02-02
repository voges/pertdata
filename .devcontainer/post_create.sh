#!/usr/bin/env bash

# Remove all editor settings for Python from the VS Code remote settings.
sed --in-place '/"\[python\]": {/,/}/d' /home/vscode/.vscode-server/data/Machine/settings.json

# Install additional packages.
sudo apt-get update && sudo apt-get install --yes shellcheck

# Install Hatch.
pip --disable-pip-version-check install hatch

# Install pertdata in editable mode.
pip --disable-pip-version-check install --editable .
