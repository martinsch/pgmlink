#!/bin/bash

# we assume that this script resides in BUILDEM_DIR
export BUILDEM_DIR=$(cd `dirname $0` && pwd)
bash $BUILDEM_DIR/bin/ilastik_gui "$@"
