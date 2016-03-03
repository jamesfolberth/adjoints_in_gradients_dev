
SVGS := $(wildcard *.svg)
PNGS := $(SVGS:%.svg=%.png)

.PHONY: clean convert all

all: convert
	./trim.sh

convert: $(PNGS)

%.png: %.svg
	convert $< $@

clean:
	rm -f *.png;
