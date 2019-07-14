all:
	g++ -std=c++11 -O2  -I ~/parser Q2a.cpp ~/parser/lexer.cpp -o b.out
	g++ -std=c++11 -O2  -I ~/parser Q2b.cpp ~/parser/lexer.cpp -o m.out
	
