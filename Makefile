all dep clean tests::
	cd brute && $(MAKE) $@ && cd .. 
	cd rosenbrock && $(MAKE) $@ && cd ..
	cd advcoordesc && $(MAKE) $@ && cd ..
	cd gridlip && $(MAKE) $@ && cd ..

doc: indent doxy

doxy:
	mkdir -p doc/html &&\
	doxygen doxy.conf

clean::
	rm -rf *~ PI* core bin/* obj/* tmp *.log
