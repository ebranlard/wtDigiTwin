all: test

test:
	python -m unittest discover -v

install:
	python -m pip install -e .

