all:
	python setup.py build_ext -i

clean:
	rm -r build; rm *.c; rm *.cpp; rm *.so; rm *.html; rm *.pyc

dist_clean:
	cd apf/base && rm *.c *.cpp *.so *.html *.pyc ; rm -r __pycache__  &

	cd apf/models && rm *.c *.cpp *.so *.html *.pyc ; rm -r __pycache__  &

	cd apf && rm -r __pycache__
	
