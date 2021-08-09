
.PHONY : all

all : docs/index.html

docs/%.html : markdown/%.Rmd
	Rscript -e 'rmarkdown::render("$<", "all", output_dir="docs")'
