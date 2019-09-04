# Copyright (C) 2009 International Business Machines. All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# Author: Andreas Waechter, 2010-05-26


DOCS = sipopt_manual.pdf

all: $(DOCS)

.SUFFIXES: .tex .pdf

.tex.pdf:
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*

clean:
	rm -f *.aux *.log *.bbl *.blg *.toc
