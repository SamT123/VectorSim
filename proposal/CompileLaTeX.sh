#!/bin/bash

# Author:  Sam Turner sat19@ic.ac.uk
# Script:  CompileLaTeX.sh
# Desc:    Compiles miniproject write up, including insertion of word count and bibliography.
# Arguments: none
# Date: March 2020


# Compile LaTeX document
xelatex "proposal.tex" > junk.txt
xelatex "proposal.tex" > junk.txt
bibtex "proposal"
xelatex "proposal.tex" > junk.txt
xelatex "proposal.tex" > junk.txt


# Move document to specified location
#mv $1.pdf $2$1.pdf 

# Open document
# open $2$1.pdf

# Remove extra files created by pdflatex
rm -f  *~ 
rm -f *.aux 
rm -f *.dvi 
rm -f *.log 
rm -f *.nav 
rm -f *.out 
rm -f *.snm 
rm -f *.toc 
rm -f *.bbl 
rm -f *.blg 
rm -f *.synctex*
