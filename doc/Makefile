sw.pdf: sw.tex bibliography.bib command.tex clean
	pdflatex -shell-escape sw
	pdflatex -shell-escape sw
	bibtex sw
	pdflatex -shell-escape sw
	pdflatex -shell-escape sw

clean:
	rm -fv *.aux *.bbl *.blg *.log *.pdf *.ps *.dvi *.tgz

archive: clean
	tar -czf sw.tgz *


