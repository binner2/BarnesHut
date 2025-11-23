#!/bin/bash
# Script to download and setup GLAD for OpenGL 4.5 Core Profile

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
GLAD_DIR="$SCRIPT_DIR/glad"

echo "========================================="
echo "  GLAD OpenGL Loader Setup"
echo "========================================="

# Check if glad directory already exists
if [ -d "$GLAD_DIR" ] && [ -f "$GLAD_DIR/glad.h" ]; then
    echo "GLAD already exists at $GLAD_DIR"
    read -p "Do you want to re-download? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Using existing GLAD installation."
        exit 0
    fi
    rm -rf "$GLAD_DIR"
fi

# Create glad directory
mkdir -p "$GLAD_DIR"
cd "$GLAD_DIR"

echo ""
echo "Checking for Python glad-generator..."

# Check if Python glad is installed
if command -v glad &> /dev/null; then
    echo "Found glad-generator, generating GLAD files..."
    glad --api="gl:core=4.5" --out-path=. --generator=c

    # Check if generation was successful
    if [ -f "glad/glad.h" ] && [ -f "src/glad.c" ]; then
        echo "Moving files to correct location..."
        mv glad/glad.h .
        mv glad/KHR .
        mv src/glad.c .
        rm -rf glad src
        echo "✓ GLAD generated successfully!"
        exit 0
    else
        echo "✗ GLAD generation failed."
    fi
else
    echo "glad-generator not found."
    echo "Install with: pip install glad"
    echo ""
fi

echo ""
echo "========================================="
echo "  Manual GLAD Setup Instructions"
echo "========================================="
echo ""
echo "Please download GLAD manually:"
echo ""
echo "1. Visit: https://glad.dav1d.de/"
echo ""
echo "2. Configure:"
echo "   - Language: C/C++"
echo "   - Specification: OpenGL"
echo "   - gl: Version 4.5 (or higher)"
echo "   - Profile: Core"
echo "   - Extensions: (none needed)"
echo "   - Options: ☑ Generate a loader"
echo ""
echo "3. Click 'Generate' and download the ZIP file"
echo ""
echo "4. Extract and copy files:"
echo "   - glad.h → $GLAD_DIR/glad.h"
echo "   - glad.c → $GLAD_DIR/glad.c"
echo "   - KHR/   → $GLAD_DIR/KHR/"
echo ""
echo "5. Alternatively, use the provided minimal GLAD (OpenGL 4.5)"
echo ""
echo "========================================="

# Create minimal GLAD header as fallback
cat > "$GLAD_DIR/glad_minimal_notice.txt" << 'EOF'
For production use, download the full GLAD loader from:
https://glad.dav1d.de/

Settings:
- API: OpenGL Core 4.5
- Generate a loader: YES

Or use:
pip install glad
glad --api="gl:core=4.5" --out-path=. --generator=c
EOF

echo ""
echo "Notice file created at: $GLAD_DIR/glad_minimal_notice.txt"
echo ""
echo "Alternative: Install glad-generator with:"
echo "  pip install glad"
echo "  $SCRIPT_DIR/download_glad.sh"
echo ""

exit 1
