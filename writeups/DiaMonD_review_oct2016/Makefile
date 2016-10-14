all:	all_texs

clean: clean_texs clean_exec

clean_texs:
	rm -f *.aux
	rm -f *.log
	rm -f *.dvi
	rm -f *.toc
	rm -f *.blg
	rm -f *.bbl
	rm -f DiaMonD_supercap*.ps
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
	pdflatex DiaMonD_supercap.tex
	pdflatex DiaMonD_supercap.tex
	pdflatex DiaMonD_supercap.tex
	pdflatex DiaMonD_supercap.tex



.SUFFIXES: .eps .png .jpg

%.eps: %.png
	convert $< $(basename $<).eps

%.eps: %.jpg
	convert $< $(basename $<).eps
