#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>       


/*void convolution_blur(int *redDest, int *greenDest, int *blueDest,
		const int *redSource, const int *greenSource, const int *blueSource,
		const int n, const int m){
	int i, j, index;
	for (i = 1; i < n-1; i++)
    	for (j = 1; j < m-1; j++){
            index = j + i * m;
           	    redDest[index]   =  (1./9.)*redSource[index-m-1]+(1./9.)*redSource[index-m]+(1./9.)*redSource[index-m+1]+
									(1./9.)*redSource[index-1  ]+(1./9.)*redSource[index  ]+(1./9.)*redSource[index+1  ]+
									(1./9.)*redSource[index+m-1]+(1./9.)*redSource[index+m]+(1./9.)*redSource[index+m+1];
							    
			    greenDest[index] =  (1./9.)*greenSource[index-m-1]+(1./9.)*greenSource[index-m]+(1./9.)*greenSource[index-m+1]+
									(1./9.)*greenSource[index-1  ]+(1./9.)*greenSource[index  ]+(1./9.)*greenSource[index+1  ]+
									(1./9.)*greenSource[index+m-1]+(1./9.)*greenSource[index+m]+(1./9.)*greenSource[index+m+1];
								
			    blueDest[index]  =  (1./9.)*blueSource[index-m-1]+(1./9.)*blueSource[index-m]+(1./9.)*blueSource[index-m+1]+
									(1./9.)*blueSource[index-1  ]+(1./9.)*blueSource[index  ]+(1./9.)*blueSource[index+1  ]+
									(1./9.)*blueSource[index+m-1]+(1./9.)*blueSource[index+m]+(1./9.)*blueSource[index+m+1];
    	}
}*/

void convolution_base(int *redDest, int *greenDest, int *blueDest,
		const int *redSource, const int *greenSource, const int *blueSource,
		const int n, const int m, float k[9]){
	int i, j, index;
	#pragma omp parallel for schedule(runtime)
	for (i = 1; i < n-1; i++)
    	for (j = 1; j < m-1; j++){
            index = j + i * m;
           	    redDest[index]   =  k[0]*redSource[index-m-1]+k[1]*redSource[index-m]+k[2]*redSource[index-m+1]+
									k[3]*redSource[index-1  ]+k[4]*redSource[index  ]+k[5]*redSource[index+1  ]+
									k[6]*redSource[index+m-1]+k[7]*redSource[index+m]+k[8]*redSource[index+m+1];
							    
			    greenDest[index] =  k[0]*greenSource[index-m-1]+k[1]*greenSource[index-m]+k[2]*greenSource[index-m+1]+
									k[3]*greenSource[index-1  ]+k[4]*greenSource[index  ]+k[5]*greenSource[index+1  ]+
									k[6]*greenSource[index+m-1]+k[7]*greenSource[index+m]+k[8]*greenSource[index+m+1];
								
			    blueDest[index]  =  k[0]*blueSource[index-m-1]+k[1]*blueSource[index-m]+k[2]*blueSource[index-m+1]+
									k[3]*blueSource[index-1  ]+k[4]*blueSource[index  ]+k[5]*blueSource[index+1  ]+
									k[6]*blueSource[index+m-1]+k[7]*blueSource[index+m]+k[8]*blueSource[index+m+1];
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
					if (redDest[index]>255)
						redDest[index]=255;
					if (greenDest[index]>255)
						greenDest[index]=255;
					if (blueDest[index]>255)
						blueDest[index]=255;
				}
				
	}
	
//Kernels for image processing
float blur[9]={1./9., 1./9., 1./9., 
	           1./9., 1./9., 1./9., 
	           1./9., 1./9., 1./9.};
	           
float copy[9]={0, 0, 0, 
			   0, 1, 0, 
			   0, 0, 0};	    
	           
float edge[9]={0, -1, 0, 
			   -1, 4, -1,
			   0, -1, 0};
			  
float sharp[9]={0, -1, 0,
			    -1, 5, -1,	
			    0, -1, 0};		
				
float bizarre[9]={-2, 0, -2, 
				  7, -5, 7,
				  -2, 0, -2};

int main(int argc, char *argv[] )
{	  
	int n, m, colorDepth, i, j, index;
	double temps_fetch, temps_convolution, temps_ecriture, temps_reflexion;
	double stop1, stop2, stop3, start1, start2, start3, stop4, start4;
	int *redSource, *greenSource, *blueSource;
	int *redDest, *greenDest, *blueDest, *greyDest;
	FILE *inFile, *outFile;
	if (argc < 2 ){
		printf("Usage: convolution input_image output_image\n");
		exit(1);
	}
//----------------------------------------------------------------------
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
    greyDest=(int*)malloc(n*m*sizeof(int));

    // Reading pixel data from file
    for (i = 0; i < n; i++)
    	for (j = 0; j < m; j++){
            index = j + i * m;
    		fscanf(inFile, "%d", &redSource[index]);
    		fscanf(inFile, "%d", &greenSource[index]);
    		fscanf(inFile, "%d", &blueSource[index]);
    	}

	stop1 = omp_get_wtime();
//----------------------------------------------------------------------
	start4 = omp_get_wtime();
	int mode=1000;
	while(mode!=1 && mode!=2 && mode!=3 && mode!=4 && mode!=5 && mode!=6)
	{
	printf("\nMode Selection\n \n1 copy	\n2 blur \n3 edge \n4 sharp \n5 noise\n6 graymap conversion\n");
	scanf("%d",&mode);
	printf("\n");
	if(mode!=1 && mode!=2 && mode!=3 && mode!=4 && mode!=5 && mode!=6)
		{
			printf("\n\n\n\n\n\n");
			printf("-------------------------------------");
			printf("\nI'd rather you put in a proper number\n");
			printf("-------------------------------------\n");
		}}
	stop4 = omp_get_wtime();
	switch(mode)
	{
		case 1 :
			{start2 = omp_get_wtime();
			convolution_base(redDest, greenDest, blueDest,
							redSource, greenSource, blueSource, n, m, copy);
			securite(redDest, greenDest, blueDest, n, m);    		
			stop2 = omp_get_wtime();}
		break;
		case 2 :
			{start2 = omp_get_wtime();
			convolution_base(redDest, greenDest, blueDest,
							redSource, greenSource, blueSource, n, m, blur);
			securite(redDest, greenDest, blueDest, n, m);    		
			stop2 = omp_get_wtime();}
		break;
		case 3 :
			{start2 = omp_get_wtime();
			convolution_base(redDest, greenDest, blueDest,
							redSource, greenSource, blueSource, n, m, edge);
			securite(redDest, greenDest, blueDest, n, m);    		
			stop2 = omp_get_wtime();}
		break;
		case 4 :
			{start2 = omp_get_wtime();
			convolution_base(redDest, greenDest, blueDest,
							redSource, greenSource, blueSource, n, m, sharp);
			securite(redDest, greenDest, blueDest, n, m);    		
			stop2 = omp_get_wtime();}
		break;
		case 5 :
			{start2 = omp_get_wtime();
			convolution_base(redDest, greenDest, blueDest,
							redSource, greenSource, blueSource, n, m, bizarre);
			securite(redDest, greenDest, blueDest, n, m);    		
			stop2 = omp_get_wtime();}
		break;
		case 6 :
			{
			#pragma omp parallel for
			for (i = 0; i < n; i++)
	    		for (j = 0; j < m; j++){
				index = j + i * m;
	    			greyDest[index] = 0.21*redSource[index]+0.72*greenSource[index]+0.07*blueSource[index];
			}
			}

		break;
}
	
//----------------------------------------------------------------------
	start3 = omp_get_wtime();
    // Opening file
	outFile = fopen(argv[2], "w");
	if (outFile == NULL) {
	  printf("Can't open output file %s\n", argv[2]);
	  exit(1);
	}
	if(mode==6){
		// Writing metadata to file for graymap
		 fprintf(outFile, "P2\n");
    		 fprintf(outFile, "%d %d\n", m, n);
   		 fprintf(outFile, "%d\n", colorDepth);

		// Writing pixel data to file
		 for (i = 0; i < n; i++){
    			for (j = 0; j < m; j++){
            			index = j + i * m;
            			fprintf(outFile, "%d ", greyDest[index]);
    			}
    		 fprintf(outFile,"\n");
    		 }
	}
	else
	{
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
	free(greyDest);

	stop3 = omp_get_wtime();
//----------------------------------------------------------------------
	temps_fetch=stop1-start1;
	temps_convolution =stop2-start2;
	temps_ecriture=stop3-start3;
	temps_reflexion=stop4-start4;
	printf("\ntemps fetch=%f s\ntemps convolution=%f s\ntemps ecriture=%f s\n\n",temps_fetch,temps_convolution,temps_ecriture);
	printf("almost forgot, it took you %f seconds to know what you wanted me to do, a bit long no ? \n\n",temps_reflexion);
	
	return 0;
}

