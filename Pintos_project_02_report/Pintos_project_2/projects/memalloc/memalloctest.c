#include <stdio.h>
#include <string.h>

#include "threads/thread.h"
#include "threads/malloc.h"
#include "threads/palloc.h"

#include "projects/memalloc/memalloctest.h"

void run_memalloc_test(char **argv UNUSED)
{
	size_t i;
	char* dynamicmem[4];

	for (i=0; i<3; i++) {
		dynamicmem[i] = (char *) malloc (145000);
		memset (dynamicmem[i], 0x00, 145000);
		printf ("dump page status : \n");
		palloc_get_status (0);
	}

	free(dynamicmem[1]);
	printf ("dump page status : \n");
	palloc_get_status (0);

	dynamicmem[3] = (char *) malloc (16000);
	memset (dynamicmem[3], 0x00, 16000);
	printf ("dump page status : \n");
	palloc_get_status (0);

	dynamicmem[1] = (char *) malloc (145000);
	memset (dynamicmem[1], 0x00, 145000);
	printf ("dump page status : \n");
	palloc_get_status (0);

	thread_sleep (100);

	for (i=0; i<4; i++) {
		free(dynamicmem[i]);
		printf ("dump page status : \n");
		palloc_get_status (0);
	}
}
