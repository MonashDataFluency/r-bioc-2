
.PHONY : all

all : docs/workshop.R docs/slides.html docs/index.html docs/slides.html docs/workshop.html

docs/%.html : markdown/%.Rmd
	Rscript -e 'rmarkdown::render("$<", output_dir="docs")'

docs/workshop.R : markdown/workshop.Rmd scripts/purify.py
	python3 scripts/purify.py <markdown/workshop.Rmd >docs/workshop.R

