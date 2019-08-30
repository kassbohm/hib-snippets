DIR=/home/kai/A_Sphinx/Teach/wb-hibbeler/source/hib-snippets/

PY_FILES = $(shell find $(DIR) -name '*.py')
IPYNB_FILES = $(patsubst $(DIR)%.py,$(DIR)%.ipynb,$(PY_FILES))

.PHONY: all

all: $(IPYNB_FILES)

./%.ipynb: ./%.py
	/usr/local/bin/ipynb-py-convert $< $@

clean:
	find  -type f -name "*.py" -exec touch {} +
