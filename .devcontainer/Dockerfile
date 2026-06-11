#FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04
ARG VARIANT="noble"
FROM nvidia/cuda:13.2.1-devel-ubuntu24.04
#FROM buildpack-deps:${VARIANT}-curl

LABEL dev.containers.features="common"

ARG VARIANT
RUN if [ "$VARIANT" = "noble" ] || [ "$VARIANT" = "resolute" ]; then \
    if id "ubuntu" &>/dev/null; then \
        echo "Deleting user 'ubuntu'  for $VARIANT" && userdel -f -r ubuntu || echo "Failed to delete ubuntu user for $VARIANT"; \
    else \
         echo "User 'ubuntu' does not exist for $VARIANT"; \
    fi; \
    fi

# [Optional] Uncomment this section to install additional OS packages.
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
        git \
        libopenmpi-dev \
        cmake \
		libhdf5-dev \
		hdf5-tools \
		less \
		python3-pip \
		libpython3-dev \
		clang-format \
		ninja-build

# Local fix for issue with timezone data at 'buildpack-deps'. Ref# https://github.com/devcontainers/images/issues/1709
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y reinstall --no-install-recommends tzdata \
    && echo "Etc/UTC" > /etc/timezone \
    && rm -rf /var/lib/apt/lists/*

# Remove NVIDIA CUDA apt sources to prevent devcontainer feature installation failures.
# CUDA is already provided by the base image; the apt repo is not needed and causes
# intermittent "mirror sync in progress" errors when devcontainer features run apt-get update.
RUN rm -f /etc/apt/sources.list.d/cuda*.list /etc/apt/sources.list.d/nvidia*.list \
    && apt-get clean


