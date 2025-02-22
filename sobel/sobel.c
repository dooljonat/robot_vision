#include <stdio.h>                          /* Sobel.c */
#include <math.h>
#include <stdlib.h>

// Note: use HI = 100 and LO = 40 for assignment

int pic[256][256];
int outpicx[256][256];
int outpicy[256][256];
int maskx[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
int masky[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};
double ival[256][256],maxival;

int main(int argc, char** argv)
{
    // Initialize variables
    int i,j,p,q,mr,sum1,sum2;
    double hi_threshold, lo_threshold;
    FILE *magnitude, *hi, *lo, *fp1, *fopen();
    char *foobar;

    // Open input file in read (binary) mode
    argc--; argv++;
    foobar = *argv;
    fp1=fopen(foobar,"rb");

    // Open output file(s) in write (binary) mode
	argc--; argv++;
	foobar = *argv;
	magnitude=fopen(foobar,"wb");

	argc--; argv++;
	foobar = *argv;
	hi=fopen(foobar,"wb");

	argc--; argv++;
	foobar = *argv;
	lo=fopen(foobar,"wb");
    
    // Get Hi & Lo threshold from args
    argc--; argv++;
	foobar = *argv;
	hi_threshold = atof(foobar);

    printf("%lf\n", hi_threshold);

    argc--; argv++;
	foobar = *argv;
	lo_threshold = atof(foobar);

    printf("%lf\n", lo_threshold);

    // Ignore the first 15 characters of the .pgm input,
    //  as these are for the header of the .pgm files
    for (i = 0; i < 15; i++)
    {
        getc(fp1);
    }

    // Read picture into array
    for (i=0;i<256;i++)
    { for (j=0;j<256;j++)
            {
                pic[i][j]  =  getc (fp1);
            }
    }

    mr = 1; // Mask radius
    // Compute convolution
    for (i=mr;i<256-mr;i++)
    { for (j=mr;j<256-mr;j++)
        {
            sum1 = 0;
            sum2 = 0;
            for (p=-mr;p<=mr;p++)
            {
                for (q=-mr;q<=mr;q++)
                {
                    sum1 += pic[i+p][j+q] * maskx[p+mr][q+mr];
                    sum2 += pic[i+p][j+q] * masky[p+mr][q+mr];
                }
            }

            // printf("%d %d \n", sum1, sum2);

            outpicx[i][j] = sum1;
            outpicy[i][j] = sum2;
        }
    }

    // Compute the magnitude of the gradient
    maxival = 0;
    for (i=mr;i<256-mr;i++)
    { 
        for (j=mr;j<256-mr;j++)
        {
            ival[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                    (outpicy[i][j]*outpicy[i][j])));
            if (ival[i][j] > maxival)
            {
                maxival = ival[i][j];
            }
        }
    }

    // Output the PGM headers
    fprintf(magnitude, "P5\n");              // val for PGM files
    fprintf(magnitude, "%d %d\n", 256, 256); // width and height
    fprintf(magnitude, "255\n");             // maximum grey value for image

    fprintf(hi, "P5\n");              // val for PGM files
    fprintf(hi, "%d %d\n", 256, 256); // width and height
    fprintf(hi, "255\n");             // maximum grey value for image

    fprintf(lo, "P5\n");              // val for PGM files
    fprintf(lo, "%d %d\n", 256, 256); // width and height
    fprintf(lo, "255\n");             // maximum grey value for image

    // Write magnitude, hi, and lo outputs
    for (i = 0; i < 256; i++) 
    {
        for (j = 0; j < 256; j++)
        {
            // Scale the image to the maximum grey value
            ival[i][j] = (ival[i][j] / maxival) * 255; 

            // Output magnitude
            fprintf(magnitude,"%c", (char)((int)(ival[i][j])));

            // If val is above hi threshold, add it to hi.pgm output
            if (ival[i][j] > hi_threshold)
            {
                fprintf(hi, "%c", (char)255);
            }
            else
            {
                fprintf(hi, "%c", (char)0);
            }

            // if val is above lo threshold, add it to lo.pgm output
            if(ival[i][j] > lo_threshold)
            {
                fprintf(lo, "%c", (char)255);
            }
            else
            {
                fprintf(lo, "%c", (char)0);
            }
        }
    }

    fclose(fp1);
    fclose(magnitude);
    fclose(hi);
    fclose(lo);
}