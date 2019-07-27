#!/usr/bin/env bash

export SOURCE_DIR=${TRAVIS_BUILD_DIR}
export BUILD_DIR=$(readlink -f ${SOURCE_DIR}/..)
export COMMIT_RANGE="${TRAVIS_COMMIT_RANGE/.../ }"
if [ -n "${TRAVIS_PULL_REQUEST_BRANCH}" ]
then
  export BUILD_LABEL="${TRAVIS_PULL_REQUEST_BRANCH}_${TRAVIS_BUILD_NUMBER}"
else
  export BUILD_LABEL="${TRAVIS_BRANCH}_${TRAVIS_BUILD_NUMBER}"
fi

cd ${SOURCE_DIR}

# Check C and C++ code with clang-format
echo "Checking formatting for commit range: ${COMMIT_RANGE}"
DIFF="$(./scripts/developer/git/git-clang-format --diff ${COMMIT_RANGE})"
if [ -n "${DIFF}" ] && [ "${DIFF}" != "no modified files to format" ]
then
  echo "clang-format:"
  echo "  Code format checks failed."
  echo "  Please run clang-format (or git clang-format) on your changes"
  echo "  before committing."
  echo "  The following changes are suggested:"
  echo "${DIFF}"
  exit 1
fi

exit 0
