#!/bin/bash
#Author: M. El Aabaribaoune

KSHRC="/home/${USER}/.kshrc"
KSHRC_OR="/home/${USER}/.kshrc.or"

# Check if the .kshrc file exists
if [ -f "$KSHRC" ]; then
    # If the .kshrc file exists, move it to .kshrc.or
    mv "$KSHRC" "$KSHRC_OR"
    echo "The file .kshrc has been moved to .kshrc.or."
fi

cat <<EOL > "$KSHRC"
# .kshrc

# Source global definitions
if [ -f /etc/kshrc ]; then
        . /etc/kshrc
fi

# Use emacs editing mode by default
set -o emacs

# User specific aliases and functions
module purge
module list
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1-GCC-7.3.0-2.30
module load netCDF/4.6.1-foss-2018b
module load netCDF-Fortran/4.4.4-foss-2018b
module list
EOL

echo "The file .kshrc has been created successfully."
