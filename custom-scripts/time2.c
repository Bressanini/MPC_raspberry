#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>


int main (int argc, char **argv) {
	unsigned long x = 0xFFFFFFFF;
	struct timespec start, finish;
	long delta_usecs;

	clock_gettime(CLOCK_MONOTONIC, &start);

	while (x--);

	clock_gettime(CLOCK_MONOTONIC, &finish);
	
	delta_usecs = (finish.tv_sec - start.tv_sec) * 1000000 + (finish.tv_nsec - start.tv_nsec) / 1000;
		
	printf("delta_usecs %dus\n", delta_usecs);
	
	return 0;
}
