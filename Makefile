indices_test: indices_test.c indices.c
	gcc indices_test.c indices.c -o indices_test.exe
	./indices_test.exe