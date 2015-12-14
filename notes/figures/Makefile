
SVGS := $(wildcard *.svg)
PNGS := $(SVGS:%.svg=%.png)

.PHONY: clean convert all

all: convert

convert: $(PNGS)

%.png: %.svg
	convert $< $@

clean:
	rm -f *.png;
