IPM_OBJ = IPM.o Main.o
AES_IPM_OBJ = common.o aes.o IPM.o AES_IPM.o Main.o

MAIN = Main.out
BLOCK = Block_IPM.out
PORTABLE = portable.out

C_FLAGS = -Wall -pedantic -g -O3
C_FLAGS_PORTABLE = -Wall -pedantic -g -O3 -static -static-libgcc -lpthread
SODIUM_FLAG = -lsodium
PROF_FILE = profile.txt
PROF_FLAG = -pg

all : $(MAIN)

$(MAIN) : $(AES_IPM_OBJ)
	gcc -o $@ $(AES_IPM_OBJ) $(SODIUM_FLAG) $(PROF_FLAG) $(C_FLAGS)

portable :$(PORTABLE)

$(PORTABLE): $(AES_IPM_OBJ)
	gcc -o $@ $(AES_IPM_OBJ) $(SODIUM_FLAG) $(PROF_FLAG) $(C_FLAGS_PORTABLE)

%.o: %.c 
	gcc -o $@ -c $< $(SODIUM_FLAG) $(PROF_FLAG) $(C_FLAGS)


run: $(MAIN)
	./$(MAIN)

prof: gmon.out #Generate profile file with GNU PROFILING
	gprof $(MAIN) gmon.out > $(PROF_FILE)

gmon.out:
	echo "run '$(MAIN)' first to get '$@'."

clean:
	rm -f *.o *.out
