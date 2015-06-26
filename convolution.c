#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

void convolution_blury(int *redDest, int *greenDest, int *blueDest,
		const int *redSource, const int *greenSource, const int *blueSource,
		const int n, const int m){
	int i, j, index;
	for (i = 0; i < n; i++)
    	for (j = 0; j < m; j++){
            index = j + i * m;
			redDest[index] = (redSource[index]*0.9+redSource[index-1]*0.9+redSource[index+1]*0.9)/2.7;
			greenDest[index] = (greenSource[index]*0.9+greenSource[index-1]*0.9+greenSource[index+1]*0.9)/2.7;
			blueDest[index] = (blueSource[index]*0.9+blueSource[index-1]*0.9+blueSource[index+1]*0.9)/2.7;
    	}
}

void convolution_sharp(int *redDest, int *greenDest, int *blueDest,
		const int *redSource, const int *greenSource, const int *blueSource,
		const int n, const int m){
	int i, j, index;
	for (i = 0; i < n; i++)
    	for (j = 0; j < m; j++){
            index = j + i * m;
			redDest[index] = (redSource[index]*0.9+redSource[index-1]*0.9+redSource[index+1]*0.9)/2.7;
			greenDest[index] = (greenSource[index]*0.9+greenSource[index-1]*0.9+greenSource[index+1]*0.9)/2.7;
			blueDest[index] = (blueSource[index]*0.9+blueSource[index-1]*0.9+blueSource[index+1]*0.9)/2.7;
    	}
}

void securite(int *redDest, int *greenDest, int *blueDest, const int n, const int m)
	{
		int i, j, index;
		for (i = 0; i < n; i++)
			for (j = 0; j < m; j++)
				{
					index = j + i * m;
					if (redDest[index]<0)
						redDest[index]=0;
					if (greenDest[index]<0)
						greenDest[index]=0;
					if (blueDest[index]<0)
						blueDest[index]=0;
				}
				
	}

int main(int argc, char *argv[] )
{
	int n, m, colorDepth, i, j, index;
	double temps_fetch, temps_convolution, temps_ecriture;
	double stop1, stop2, stop3, start1, start2, start3;
	int *redSource, *greenSource, *blueSource;
	int *redDest, *greenDest, *blueDest;
	FILE *inFile, *outFile;
	if (argc < 2 ){
		printf("Usage: convolution input_image output_image\n");
		exit(1);
	}

	start1 = omp_get_wtime();

	// Opening input file
	inFile = fopen(argv[1], "r");
	if (inFile == NULL) {
	  printf("Can't open input file %s\n", argv[1]);
	  exit(1);
	}

	// Reading magic number
	char magic_number[5];
    fscanf(inFile, "%s", magic_number);
    if (strcmp(magic_number, "P3")){
  	  printf("Error while reading file %s\n", argv[1]);
  	  exit(1);
    }
    // Reading image size
    fscanf(inFile, "%d", &m);
    fscanf(inFile, "%d", &n);
    // Reading color depth
    fscanf(inFile, "%d", &colorDepth);

    // Allocating memory
    redSource = (int*)malloc(n * m * sizeof(int));
    greenSource = (int*)malloc(n * m * sizeof(int));
    blueSource = (int*)malloc(n * m * sizeof(int));
    redDest = (int*)malloc(n * m * sizeof(int));
    greenDest = (int*)malloc(n * m * sizeof(int));
    blueDest = (int*)malloc(n * m * sizeof(int));

    // Reading pixel data from file
    for (i = 0; i < n; i++)
    	for (j = 0; j < m; j++){
            index = j + i * m;
    		fscanf(inFile, "%d", &redSource[index]);
    		fscanf(inFile, "%d", &greenSource[index]);
    		fscanf(inFile, "%d", &blueSource[index]);
    	}

	stop1 = omp_get_wtime();
	start2 = omp_get_wtime();

    convolution_blury(redDest, greenDest, blueDest,
    		redSource, greenSource, blueSource, n, m);
    securite(redDest, greenDest, blueDest, n, m);
    		
    stop2 = omp_get_wtime();
	start3 = omp_get_wtime();
    // Opening file
	outFile = fopen(argv[2], "w");
	if (outFile == NULL) {
	  printf("Can't open output file %s\n", argv[2]);
	  exit(1);
	}
    // Writing metadata to file
	fprintf(outFile, "P3\n");
    fprintf(outFile, "%d %d\n", m, n);
    fprintf(outFile, "%d\n", colorDepth);

    // Writing pixel data to file
    for (i = 0; i < n; i++){
    	for (j = 0; j < m; j++){
            index = j + i * m;
            fprintf(outFile, "%d ", redDest[index]);
            fprintf(outFile, "%d ", greenDest[index]);
            fprintf(outFile, "%d ", blueDest[index]);
    	}
    	fprintf(outFile,"\n");
    }

	// Closing files
	fclose(inFile);
	fclose(outFile);

	// Deallocating memory
	free(redSource);
	free(greenSource);
	free(blueSource);
	free(redDest);
	free(greenDest);
	free(blueDest);
	stop3 = omp_get_wtime();

	temps_fetch=stop1-start1;
	temps_convolution =stop2-start2;
	temps_ecriture=stop3-start3;
	printf("\ntemps fetch=%f\ntemps convolution=%f\ntemps ecriture=%f\n\n",temps_fetch,temps_convolution,temps_ecriture);
	
	return 0;
}

