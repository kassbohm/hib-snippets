PY_DIR=/home/kai/A_Sphinx/Teach/wb-hibbeler/source/hib-snippets/de/py/
IPYNB_DIR=/home/kai/A_Sphinx/Teach/wb-hibbeler/source/hib-snippets/de/ipynb/

PY_FILES = $(shell find $(PY_DIR) -name '*.py')
IPYNB_FILES = $(patsubst $(PY_DIR)%.py, $(IPYNB_DIR)%.ipynb, $(PY_FILES))

.PHONY: all

all: $(IPYNB_FILES)

$(IPYNB_DIR)%.ipynb: $(PY_DIR)%.py
	/usr/local/bin/ipynb-py-convert $< $@

clean:
	find  -type f -name "*.py" -exec touch {} +
