#!/usr/bin/env bash
# Generic Apptainer builder: turns any conda-style spec file (same schema
# as environment.yml -- name/channels/dependencies, whatever it's called
# or extensioned) into a .sif image, using the shared apptainer.def
# template next to this script.
#
# Usage:
#   ./build_apptainer.sh <name_of_sif> <name_of_specification_file> [env_name] [base_image]
#
# Examples:
#   ./build_apptainer.sh bioenv.sif environment.yml
#   ./build_apptainer.sh bioenv_sruge.sif sruge_scripts.spec bioenv
#
# Run this with bash, not `sh` -- it uses bash-only features (arrays,
# BASH_SOURCE). Either `./build_apptainer.sh ...` (after chmod +x) or
# `bash build_apptainer.sh ...` both work; plain `sh build_apptainer.sh`
# may invoke dash and fail.
set -euo pipefail
 
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEF_FILE="$SCRIPT_DIR/apptainer.def"
 
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <name_of_sif> <name_of_specification_file> [env_name] [base_image]" >&2
    exit 1
fi
 
SIF_NAME="$1"
SPEC_FILE="$2"
ENV_NAME="${3:-env}"
BASE_IMAGE="${4:-condaforge/miniforge3:latest}"
 
# Append .sif if the caller didn't include it.
case "$SIF_NAME" in
    *.sif) ;;
    *) SIF_NAME="${SIF_NAME}.sif" ;;
esac
 
if [ ! -f "$DEF_FILE" ]; then
    echo "Error: $DEF_FILE not found next to this script." >&2
    exit 1
fi
 
if [ ! -f "$SPEC_FILE" ]; then
    echo "Error: spec file '$SPEC_FILE' not found." >&2
    exit 1
fi
SPEC_FILE="$(realpath "$SPEC_FILE")"
 
if ! command -v apptainer &> /dev/null; then
    echo "Error: 'apptainer' not found in PATH." >&2
    echo "Install it first: https://apptainer.org/docs/admin/main/installation.html" >&2
    exit 1
fi
# Resolve the full path now, while we still have the invoking user's PATH
# (e.g. a conda env's bin/) -- sudo resets PATH to its own secure_path and
# won't see a conda-installed binary otherwise.
APPTAINER_BIN="$(command -v apptainer)"
 
BUILD_ARGS=(--build-arg "SPEC_FILE=$SPEC_FILE" --build-arg "ENV_NAME=$ENV_NAME" --build-arg "BASE_IMAGE=$BASE_IMAGE")
 
echo "Building $SIF_NAME"
echo "  spec file : $SPEC_FILE"
echo "  env name  : $ENV_NAME"
echo "  base image: $BASE_IMAGE"
echo "  apptainer : $APPTAINER_BIN"
echo ""
 
if [ "$(id -u)" -eq 0 ]; then
    # Already root (e.g. invoked via `sudo -E ./build_apptainer.sh ...`).
    "$APPTAINER_BIN" build "${BUILD_ARGS[@]}" "$SIF_NAME" "$DEF_FILE"
elif sudo -n true 2>/dev/null || sudo -v 2>/dev/null; then
    # Real root via sudo -- no --fakeroot needed. Preserve PATH so sudo
    # can find the resolved binary above instead of "command not found".
    echo "Using sudo (real root) -- no --fakeroot needed."
    sudo env "PATH=$PATH" "$APPTAINER_BIN" build "${BUILD_ARGS[@]}" "$SIF_NAME" "$DEF_FILE"
else
    # No sudo available (the common case on shared HPC login nodes) --
    # fall back to unprivileged fakeroot. Requires fakeroot to be enabled
    # for your account (`apptainer config fakeroot --add <user>`, run by
    # an admin if you don't have sudo yourself).
    echo "No sudo access -- building with --fakeroot."
    "$APPTAINER_BIN" build --fakeroot "${BUILD_ARGS[@]}" "$SIF_NAME" "$DEF_FILE"
fi
 
echo ""
echo "Build complete: $SIF_NAME"
echo "Verifying environment '$ENV_NAME' inside the container..."
# Deliberately avoid `conda activate` here: under `bash -c` inside a
# container it can silently block on stdin (e.g. a first-run channel
# Terms-of-Service prompt some conda/mamba versions show, with no visible
# cue that it's waiting). `conda list -p <prefix>` reads the env directly
# with no activation needed. `< /dev/null` guarantees nothing can wait on
# input, and `timeout` guarantees this step can't hang the script forever
# even so -- if it fails or times out, that's reported as a warning, not
# a fatal error, since the .sif was already built successfully above.
ENV_PREFIX="/opt/conda/envs/$ENV_NAME"
if ! timeout 120 "$APPTAINER_BIN" exec "$SIF_NAME" "$ENV_PREFIX/bin/python3" -V < /dev/null; then
    echo "Warning: could not run python3 in '$ENV_NAME' (failed or timed out after 120s)." >&2
fi
if ! timeout 120 "$APPTAINER_BIN" exec "$SIF_NAME" /opt/conda/bin/conda list -p "$ENV_PREFIX" < /dev/null; then
    echo "" >&2
    echo "Warning: 'conda list' verification failed or timed out after 120s." >&2
    echo "The .sif file was still built successfully: $SIF_NAME" >&2
    echo "Debug manually with, e.g.:" >&2
    echo "  $APPTAINER_BIN exec $SIF_NAME $ENV_PREFIX/bin/python3 -V" >&2
fi
