# Script to create a report from mzTab file

#!/bin/sh

# function for absolute path of a file (replacement for readlink in Linux)
# see https://stackoverflow.com/questions/3572030/bash-script-absolute-path-with-osx
realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

# script directory
SCRIPT_DIRECTORY=$(dirname $(realpath $0))
${SCRIPT_DIRECTORY}/MQ2mzTab "$@"
