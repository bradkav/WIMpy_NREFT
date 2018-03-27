objects = WM.so WMP2.so WP1.so WP2.so WS1.so WS2.so WS1D.so WD.so


all: $(objects)

clean:
	rm $(objects)

%.so: %.f
<<<<<<< HEAD
	f2py -c -L/anaconda/lib/ -lgfortran $^ -m $(basename $@) 
=======
	python2 -m numpy.f2py -c $^ -m $(basename $@) 
>>>>>>> b5c746db130b0d48e6621d1f3b8374e28db66f32
