#!/usr/bin/env bash

# Remove all editor settings for Python from the VS Code remote settings.
sed --in-place '/"\[python\]": {/,/}/d' /home/vscode/.vscode-server/data/Machine/settings.json

# Install additional packages.
sudo apt-get update && sudo apt-get install --yes shellcheck

# Install the pertdata package.
python3 -m venv .venv
# shellcheck disable=SC1091
source .venv/bin/activate
pip3 --disable-pip-version-check install --requirement requirements.txt
pip3 install --editable .
