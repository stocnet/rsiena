#!/bin/bash
set -e # Stop execution immediately if any command fails

echo "Step 1: Running update script from within the src/ directory..."
(cd src && ./update_sources_list.sh)

echo "Step 2: Running autoconf from the root directory..."
# Generates the fresh Mac/Linux Makevars
autoconf
./configure

echo "Step 3: Synchronizing the SOURCES line to Makevars.win..."
# 1. Extract the clean SOURCES line from the freshly generated Mac Makevars
NEW_SOURCES_LINE=$(grep "^SOURCES =" src/Makevars)

# 2. Use sed to replace the old SOURCES line inside your tracked Makevars.win
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS-specific sed configuration
    sed -i "" "s|^SOURCES =.*|$NEW_SOURCES_LINE|" src/Makevars.win
else
    # Linux-specific sed configuration
    sed -i "s|^SOURCES =.*|$NEW_SOURCES_LINE|" src/Makevars.win
fi

echo "Success! All build files are updated and synchronized perfectly."
