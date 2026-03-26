#!/bin/bash
set -e

# Build script for CoverM containers (Docker and Apptainer)
# Requires: cargo build --release to be run first

VERSION="0.10.0-spatial"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Building CoverM $VERSION containers ==="

# Check if binary exists
if [ ! -f "$SCRIPT_DIR/target/release/coverm" ]; then
    echo "ERROR: Binary not found at $SCRIPT_DIR/target/release/coverm"
    echo "Please run 'cargo build --release' first"
    exit 1
fi

echo ""
echo "=== Building Docker image ==="
cd "$SCRIPT_DIR"
docker build -t coverm:$VERSION .
echo "Docker image built: coverm:$VERSION"

echo ""
echo "=== Building Apptainer image ==="
# Create temporary directory in /tmp (accessible in fakeroot mode)
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Copy necessary files to temp directory
cp "$SCRIPT_DIR/target/release/coverm" "$TMPDIR/"
cp "$SCRIPT_DIR/coverm.def" "$TMPDIR/"
cp "$SCRIPT_DIR/conda.yml" "$TMPDIR/"

cd "$TMPDIR"
# Unset APPTAINER_BINDPATH to prevent mounting non-existent paths
unset APPTAINER_BINDPATH
apptainer build --fakeroot "coverm_${VERSION}.sif" coverm.def

# Move result back
mv "coverm_${VERSION}.sif" "$SCRIPT_DIR/"

echo ""
echo "=== Build complete ==="
echo "Docker image: coverm:$VERSION"
echo "Apptainer image: $SCRIPT_DIR/coverm_${VERSION}.sif"
