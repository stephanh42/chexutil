# Makefile so that we can do :make form within Vim
.PHONY: target
target: chexutil.cpp
	python3 setup.py build_ext --inplace

