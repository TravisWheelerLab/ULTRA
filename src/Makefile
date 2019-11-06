all: ultra

ultra: *.cpp *.hpp
	@$(CXX) -std=c++11 $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o ../$@ *.cpp

clean:
	@rm -f ../ultra
	@rm -f ../profile