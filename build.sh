#!/bin/bash
set -euo pipefail

# Quick local cross-build helper. Tagged releases are built by GoReleaser
# (see .goreleaser.yaml / .github/workflows/release.yml).

mkdir -p dist/linux dist/windows

LDFLAGS="-s -w"

echo "Building for Linux (amd64)..."
GOOS=linux GOARCH=amd64 go build -ldflags "$LDFLAGS" -o dist/linux/admix ./cmd/admix
GOOS=linux GOARCH=amd64 go build -ldflags "$LDFLAGS" -o dist/linux/server ./cmd/server

echo "Building for Windows (amd64)..."
GOOS=windows GOARCH=amd64 go build -ldflags "$LDFLAGS" -o dist/windows/admix.exe ./cmd/admix
GOOS=windows GOARCH=amd64 go build -ldflags "$LDFLAGS" -o dist/windows/server.exe ./cmd/server

echo "Build complete! Check the 'dist' directory."
