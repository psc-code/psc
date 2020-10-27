#!/usr/bin/env bash

set -e
set -x
shopt -s dotglob

readonly name="gtensor"
readonly ownership="gtensor upstream <robot@gene>"
readonly subtree="external/gtensor/gtensor"
readonly repo="https://github.com/wdmapp/gtensor.git"
readonly tag="master"
readonly shortlog="true"
readonly paths="
"

extract_source () {
    git_archive
}

. "${BASH_SOURCE%/*}/../update-common.sh"
