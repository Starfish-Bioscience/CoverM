#!/bin/bash
set -e

# Load user environment (provides ORAS_TOKEN etc.)
# shellcheck source=/dev/null
[ -f "$HOME/.bashrc" ] && source "$HOME/.bashrc"

# Build script for CoverM Apptainer container
# Requires: cargo build --release to be run first
#
# Usage:
#   ./build-container.sh                          # Build → container/
#   ./build-container.sh --install                # Build → /HDD/singularity_apptainer/
#   ./build-container.sh --install --push         # Build, install, and push to GHCR
#   ./build-container.sh --dry-run [--install] [--push]  # Show what would happen, no build/push

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"

# Derive version from Cargo.toml + git branch suffix (feature/spatial-metrics → spatial-metrics)
CARGO_VERSION=$(grep '^version' "$REPO_DIR/Cargo.toml" | head -1 | sed 's/version = "\(.*\)"/\1/')
GIT_BRANCH=$(git -C "$REPO_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null || echo "")
if [[ "$GIT_BRANCH" =~ ^(feature|bugfix|fix|release|hotfix)/ ]]; then
    BRANCH_SUFFIX="${GIT_BRANCH#*/}"
    VERSION="${CARGO_VERSION}-${BRANCH_SUFFIX}"
else
    VERSION="${CARGO_VERSION}"
fi

# GHCR configuration
GHCR_OWNER="starfish-bioscience"
GHCR_IMAGE="coverm"
GHCR_REGISTRY="ghcr.io/${GHCR_OWNER}/${GHCR_IMAGE}"

# Parse arguments
PUSH_TO_GHCR=false
INSTALL=false
DRY_RUN=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --push)
            PUSH_TO_GHCR=true
            shift
            ;;
        --install)
            INSTALL=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--install] [--push] [--dry-run]"
            exit 1
            ;;
    esac
done

# SIF destination: container/ by default, /HDD/singularity_apptainer/ with --install
if [ "$INSTALL" = true ]; then
    SIF_DEST="/HDD/singularity_apptainer"
else
    SIF_DEST="$SCRIPT_DIR"
fi
SIF_FILE="${SIF_DEST}/coverm_${VERSION}.sif"

if [ "$DRY_RUN" = true ]; then
    echo "=== Dry run — nothing will be built or pushed ==="
    echo "  Version   : $VERSION"
    echo "  SIF file  : $SIF_FILE"
    if [ "$PUSH_TO_GHCR" = true ]; then
        echo "  GHCR push : oras://${GHCR_REGISTRY}:${VERSION}"
        echo "              oras://${GHCR_REGISTRY}:latest"
    else
        echo "  GHCR push : (pass --push to enable)"
    fi
    exit 0
fi

echo "=== Building CoverM $VERSION Apptainer container ==="

# Check if binary exists
if [ ! -f "$REPO_DIR/target/release/coverm" ]; then
    echo "ERROR: Binary not found at $REPO_DIR/target/release/coverm"
    echo "Please run 'cargo build --release' first"
    exit 1
fi

# Create temporary directory in /tmp (accessible in fakeroot mode)
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Copy necessary files to temp directory
cp "$REPO_DIR/target/release/coverm" "$TMPDIR/"
cp "$SCRIPT_DIR/coverm.def" "$TMPDIR/"
cp "$SCRIPT_DIR/conda.yml" "$TMPDIR/"

cd "$TMPDIR"
# Unset APPTAINER_BINDPATH to prevent mounting non-existent paths
unset APPTAINER_BINDPATH
# apptainer build may exit with 1 due to a getcwd() failure after tmpdir cleanup
# (known false positive). Check the SIF was actually created before continuing.
apptainer build --fakeroot "coverm_${VERSION}.sif" coverm.def || true
if [ ! -f "coverm_${VERSION}.sif" ]; then
    echo "ERROR: apptainer build failed — SIF file not created."
    exit 1
fi

# Move result to permanent location
mv "coverm_${VERSION}.sif" "$SIF_FILE"

echo ""
echo "=== Build complete ==="
echo "Apptainer image: $SIF_FILE"
echo ""
echo "Test the container with:"
echo "  apptainer exec $SIF_FILE coverm --version"

# Push to GHCR if requested
if [ "$PUSH_TO_GHCR" = true ]; then
    echo ""
    echo "=== Pushing to GitHub Container Registry ==="
    echo "Registry: $GHCR_REGISTRY:$VERSION"
    echo ""

    # ORAS_TOKEN is a dedicated PAT in ~/.bashrc for container registry pushes,
    # separate from the gh CLI token. APPTAINER_LOGIN_PASSWORD / _USERNAME are
    # the native Apptainer env vars read automatically by `apptainer remote login`.
    if [ -z "$ORAS_TOKEN" ]; then
        echo "ERROR: ORAS_TOKEN is not set. Define it in ~/.bashrc."
        exit 1
    fi
    GHCR_USERNAME=$(gh api user --jq '.login' 2>/dev/null || echo "")
    if [ -z "$GHCR_USERNAME" ]; then
        echo "ERROR: Could not determine GitHub username. Make sure 'gh' is installed and authenticated."
        exit 1
    fi
    echo "$ORAS_TOKEN" | apptainer registry login --username "$GHCR_USERNAME" --password-stdin oras://ghcr.io

    apptainer push "$SIF_FILE" "oras://${GHCR_REGISTRY}:${VERSION}"
    apptainer push "$SIF_FILE" "oras://${GHCR_REGISTRY}:latest"

    echo ""
    echo "=== Push complete ==="
    echo "Image available at:"
    echo "  oras://${GHCR_REGISTRY}:${VERSION}"
    echo "  oras://${GHCR_REGISTRY}:latest"
    echo ""
    echo "Pull with:"
    echo "  apptainer pull oras://${GHCR_REGISTRY}:${VERSION}"
fi
