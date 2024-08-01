#!/bin/sh
# Linux wrapper for all EnzyDock Python scripts. The following must be set correctly:
# $PYTHON4ENZYDOCK must point to a Python interpreter with the following installed features:
# rdkit numpy openbabel pandas pymol sklearn kneed 
# Written by Dan T Major
# Date: 12/09/2022 
# Copyright © 2022 Dan T. Major

PYTHON4ENZYDOCK=/home/qnt/majort/anaconda3/envs/my-rdkit-env/bin/python3.9
SCRIPT=$1
shift
exec "$PYTHON4ENZYDOCK" "$SCRIPT" "$@"
