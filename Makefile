CXXFLAGS = -O3 -lpthread

all:
	@cd src && ${MAKE} CXXFLAGS="${CXXFLAGS}"

prefix = /usr/local
exec_prefix = ${prefix}
bindir = ${exec_prefix}/bin

install: all
	mkdir -p ${bindir}
	cp src/ ${bindir}

profile:
	@cd src && ${MAKE} profile CXXFLAGS="${CXXFLAGS}"

clean:
	@cd src && ${MAKE} clean