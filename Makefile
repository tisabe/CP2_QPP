indices_test: indices_test.c indices.c
	gcc indices_test.c indices.c -o indices_test.exe -02 -std=c99
	./indices_test.exe