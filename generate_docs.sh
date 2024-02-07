#!/bin/bash

# Package name
PACKAGE_NAME="autostreamtree"

# Generate documentation in HTML
pdoc --html --output-dir ./docs $PACKAGE_NAME

# Assuming pdoc generates a folder with the package name, move it to the root of /docs
mv ./docs/$PACKAGE_NAME/* ./docs/

echo "HTML documentation generated in ./docs/"


