#!/bin/zsh
# Shell script to convert source files into a markdown
# document. Most of the heavy lifting done by pandoc:
#   https://pandoc.org

pandoc \
    --from markdown \
    --to gfm-raw_html \
    --wrap=preserve \
    --citeproc \
    --bibliography ../Source/References.bib \
    --lua-filter=relocate_figs.lua \
    --output ../chapter.md \
    ../Source/text.md
