#!/bin/bash

# Package name
PACKAGE_NAME="autostreamtree"

# Generate documentation in HTML
pdoc3 --html --output-dir ./docs ./$PACKAGE_NAME
mv ./docs/autostreamtree/* ./docs
rm -rf ./docs/autostreamtree

echo "HTML documentation generated in ./docs/"


