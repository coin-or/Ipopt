# Copyright (C) 2009 International Business Machines. All Rights Reserved.
# This code is published under the Common Public License.
#
# Author: Andreas Waechter, 2010-05-26


DOCS = overview.pdf usage.pdf

all: $(DOCS)

.SUFFIXES: .tex .pdf

.tex.pdf:
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*

clean:
	rm -f *.aux *.log *.bbl *.blg
