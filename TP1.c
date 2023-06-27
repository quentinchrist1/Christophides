#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "TP1Functions.h"

int main()
{
	//Structure
	graph G;

	// Nom du fichier en lecture
	char* instancename = (char*)malloc(256*sizeof(char));
	strcpy(instancename,"K50.csv");

	readTP1Instance(&G,instancename);
	printf("Starting graph\n");
	display_graph(&G);
	Christophides(&G);
	free_graph(&G);
	free(instancename);

	return 0;
}

