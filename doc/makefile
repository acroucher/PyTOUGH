FILENAME=PyTOUGH-guide

$(FILENAME).pdf: $(FILENAME).dvi
	dvipdf $(FILENAME).dvi

$(FILENAME).dvi: *.tex *.bib
	latex $(FILENAME)
	makeindex $(FILENAME)
	bibtex $(FILENAME)
	latex $(FILENAME)
	latex $(FILENAME)