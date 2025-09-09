#!/bin/bash

# Get the path of the evodesign repository
SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

# Install the anaconda environment
conda env create -f $SCRIPT_DIR/environment.yml

# Path to the evodesign source code
NEW_PATH=$SCRIPT_DIR/evodesign

# Line to add in ~/.bashrc
LINE="export PYTHONPATH=\$PYTHONPATH:$NEW_PATH"

# Check if ~/.bashrc already contains this exact path
if grep -Fxq "$LINE" ~/.bashrc; then
    echo "PYTHONPATH already includes $NEW_PATH"
else
    echo "$LINE" >> ~/.bashrc
    echo "Added $NEW_PATH to PYTHONPATH in ~/.bashrc"
fi

# Reload ~/.bashrc for current session
source ~/.bashrc
