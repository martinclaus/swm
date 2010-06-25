sw.pdf: sw.tex bibliography.bib command.tex clean
	pdflatex sw
	pdflatex sw
	bibtex sw
	pdflatex sw
	pdflatex sw

clean:
	rm -fv *.aux *.bbl *.blg *.log *.pdf *.ps *.dvi *.tgz

archive: clean
	tar -czf sw.tgz *


