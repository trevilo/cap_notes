all:	all_texs

clean: clean_texs clean_exec

clean_texs:
	rm -f *.aux
	rm -f *.log
	rm -f *.dvi
	rm -f *.toc
	rm -f *.blg
	rm -f *.bbl
	rm -f resamplig_slides*.ps
	rm -f *.out
	rm -f *.nav
	rm -f *.snm
	rm -f *.synctex.gz

clean_exec:
	rm -f *~

#OBJS_convert = \
#	rawfigs/xyz.eps \
#	rawfigs/abc.eps

# $(OBJS_convert)

all_texs:
	pdflatex talk.tex
	pdflatex talk.tex



.SUFFIXES: .eps .png .jpg

%.eps: %.png
	convert $< $(basename $<).eps

%.eps: %.jpg
	convert $< $(basename $<).eps
