bbr.pdf: bbr.tex \
R/slide.tex overview.tex algebra.tex descript/slide.tex htest/slide.tex \
prop/slide.tex nonpar/slide.tex corr/slide.tex rmsintro/slide.tex \
reg/slide.tex multgroup/slide.tex infreview.tex ancova/slide.tex \
change/slide.tex serial/slide.tex obsvar/slide.tex propensity.tex \
info/slide.tex dx/slide.tex hdata/slide.tex repro/slide.tex
	pdflatex bbr

R/slide.tex: R/slide.Rnw
	knitrnf $(<D) $(<F)

descript/slide.tex: descript/slide.Rnw
	knitrnf $(<D) $(<F)

htest/slide.tex: htest/slide.Rnw
	knitrnf $(<D) $(<F)

prop/slide.tex: prop/slide.Rnw
	knitrnf $(<D) $(<F)

nonpar/slide.tex: nonpar/slide.Rnw
	knitrnf $(<D) $(<F)

corr/slide.tex: corr/slide.Rnw
	knitrnf $(<D) $(<F)

rmsintro/slide.tex: rmsintro/slide.Rnw
	knitrnf $(<D) $(<F)

reg/slide.tex: reg/slide.Rnw
	knitrnf $(<D) $(<F)

multgroup/slide.tex: multgroup/slide.Rnw
	knitrnf $(<D) $(<F)

ancova/slide.tex: ancova/slide.Rnw
	knitrnf $(<D) $(<F)

change/slide.tex: change/slide.Rnw
	knitrnf $(<D) $(<F)

serial/slide.tex: serial/slide.Rnw
	knitrnf $(<D) $(<F)

obsvar/slide.tex: obsvar/slide.Rnw
	knitrnf $(<D) $(<F)

info/slide.tex: info/slide.Rnw
	knitrnf $(<D) $(<F)

dx/slide.tex: dx/slide.Rnw
	knitrnf $(<D) $(<F)

hdata/slide.tex: hdata/slide.Rnw
	knitrnf $(<D) $(<F)

repro/slide.tex: repro/slide.Rnw
	knitrnf $(<D) $(<F)
