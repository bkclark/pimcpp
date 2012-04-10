all:	libpspline.a

libpspline.a: evtricub.o  herm3ev.o  mktricubw.o  tcspline.o bcspeval.o \
	bcspline.o ibc_ck.o splinck.o tcspeval.o v_spline.o zonfind.o 
	ar r libpspline.a *.o
	ranlib libpspline.a

all:	evtricub.o  herm3ev.o  mktricubw.o  tcspline.o ibc_ck.o splinck.o\
        zonfind.o tcspeval.o v_spline.o bcspline.o bcspeval.o

clean:
	rm -f *.o

.f.o:
	$(F77) -g -O3 -c -o $*.o $<
