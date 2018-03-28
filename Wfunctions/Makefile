objects = WM.so WMP2.so WP1.so WP2.so WS1.so WS2.so WS1D.so WD.so


all: $(objects)

clean:
	rm $(objects)

%.so: %.f
#	f2py -c -L/anaconda/lib/ -lgfortran $^ -m $(basename $@) 
	python2 -m numpy.f2py -c $^ -m $(basename $@) 

