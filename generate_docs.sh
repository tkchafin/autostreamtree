#!/bin/bash

# Package name
PACKAGE_NAME="autostreamtree"

# Version number 
VERSION_NUMBER="v1.1.0"

PDF_OUTPUT_NAME="${PACKAGE_NAME}_${VERSION_NUMBER}_$(date +%Y-%m-%d).pdf"

# Generate documentation and convert it to PDF
pdoc --pdf $PACKAGE_NAME > pdf.md 
pandoc \
    --metadata=title:"${PACKAGE_NAME} ${VERSION_NUMBER} $(date +%Y-%m-%d)" \
    --from=markdown+abbreviations+tex_math_single_backslash \
    --pdf-engine=xelatex \
    --variable=mainfont:"Helvetica" \
    --toc \
    --toc-depth=4 \
    --output=./docs/$PDF_OUTPUT_NAME pdf.md
rm pdf.md

echo "PDF documentation generated at ./docs/$PDF_OUTPUT_NAME"


