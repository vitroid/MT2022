all: Figure1.pdf Figure2.pdf Figure3.pdf Figure4.pdf
Figure%.pdf: Figure%.py LJparam.py
	python $<

pep8:
	autopep8 -r -a -a -i .
