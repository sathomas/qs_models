#!/bin/zsh
# Shell script to convert source files into a MS-Word
# document. Most of the heavy lifting done by pandoc:
#   https://pandoc.org

pandoc \
    --from markdown \
    --to docx \
    --reference-doc reference.docx \
    --citeproc \
    --csl springer-basic-brackets.csl \
    --bibliography ../Source/References.bib \
    --lua-filter=remove_figs.lua \
    --output ../Submission/chapter.docx \
    ../Source/text.md
