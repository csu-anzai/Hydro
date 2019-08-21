OKULAREXE=C:/Users/tah09e/AppData/KDE/bin/okular.EXE
MAINFILE=polyhedraVolume_writeup

all: clean
	
	pdflatex  ${MAINFILE}
#	bibtex	  ${MAINFILE}
	pdflatex  ${MAINFILE}
#	pdflatex  ${MAINFILE}
	okular ${MAINFILE}.pdf

pdf: all
	dvips -Ppdf main.dvi
	ps2pdf main.ps
#	${OKULAREXE} main.pdf
	okular main.pdf

dvi: all
	okular main.dvi
	
clean:
	
	rm -rf *aux *spl *pdf *log *bbl *blg *dvi *ps
	
