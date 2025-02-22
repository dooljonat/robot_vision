/* Jonathan Dooley - January 20th 2024 */

#include <stdio.h>                  /*  Canny algorithm (based on marrh.c) */
#include <math.h>
#include <stdlib.h>
#define  PICSIZE 256
#define  MAXMASK 100

#define  LOW_PERCENT 0.35

int    pic[PICSIZE][PICSIZE];
double outpic1x[PICSIZE][PICSIZE];
double outpic1y[PICSIZE][PICSIZE];
double outpic2[PICSIZE][PICSIZE];
int    peaks[PICSIZE][PICSIZE];
double maskx[MAXMASK][MAXMASK];
double masky[MAXMASK][MAXMASK];
double convx[PICSIZE][PICSIZE];
double convy[PICSIZE][PICSIZE];
double magnitude[PICSIZE][PICSIZE];
int    histogram[PICSIZE];
double final[PICSIZE][PICSIZE];

int main(int argc, char **argv)
{
    int     i, j, p, q, s, t, mr, centx, centy, area, boolean_moretodo;
    double  xmaskval, ymaskval, xsum, ysum, sigma, maxival, minival, maxval, slope, cutOff, HI, LO;
    FILE    *fo1, *fo2, *fo3, *fp1, *fopen();
    char    *foobar;

    // Check for correct number of arguments
    if (argc != 3)
    {
        printf("Error! Incorrect number of arguments passed to the program.\n");
        printf("Usage: ./canny <input>.pgm <sigma_value>\n");
        return 1;
    }

    // Scan in input file
    argc--; argv++;
    foobar = *argv;
    fp1 = fopen(foobar,"rb");

    // Scan in sigma value
    argc--; argv++;
    foobar = *argv;
    sigma = atof(foobar);

    // Open two output files
    fo1 = fopen("mag_out.pgm", "wb");
    fo2 = fopen("peaks_out.pgm", "wb");
    fo3 = fopen("final_out.pgm", "wb");

    // Create headers for output files
    fprintf(fo1, "P5\n");
    fprintf(fo1, "%d %d\n", PICSIZE, PICSIZE);
    fprintf(fo1, "255\n");

    fprintf(fo2, "P5\n");
    fprintf(fo2, "%d %d\n", PICSIZE, PICSIZE);
    fprintf(fo2, "255\n");

    fprintf(fo3, "P5\n");
    fprintf(fo3, "%d %d\n", PICSIZE, PICSIZE);
    fprintf(fo3, "255\n");

    // Flexible size masking
    mr = (int)(sigma * 3);
    centx = (MAXMASK / 2);
    centy = (MAXMASK / 2);

    // Ignore the first 15 characters of the .pgm input,
    //  as these are for the header of the .pgm files
    for (i = 0; i < 15; i++)
    {
        getc(fp1);
    }

    // Scan pixel values of input into pic array
    for (i = 0; i < PICSIZE; i++)
    {
        for (j = 0; j < PICSIZE; j++)
        {
            pic[i][j] = getc(fp1);
        }
    }

    // ** Part 1 ** //

    // Create the flexible mask table with the first 
    // x and y derivatives of the gaussian function
    for (p = -mr; p <= mr; p++)
    {
        for (q = -mr; q <= mr; q++)
        {
            // Compute masks using first x and y derivative
            // of Gaussian
            xmaskval = q * (exp(-1*( ( (p*p) + (q*q) ) / (2*(sigma*sigma)))));
            ymaskval = p * (exp(-1*( ( (p*p) + (q*q) ) / (2*(sigma*sigma)))));

            maskx[p+centy][q+centx] = xmaskval;
            masky[p+centy][q+centx] = ymaskval;
        }
    }

    // Convolution of the file with the masks created above
    for (i = mr; i <= 255 - mr; i++)
    {
        for (j = mr; j <= 255 - mr; j++)
        {
            xsum = 0;
            ysum = 0;

            for (p = -mr; p <= mr; p++)
            {
                for (q = -mr; q <= mr; q++)
                {
                    xsum += pic[i+p][j+q] * maskx[p+centy][q+centx];
                    ysum += pic[i+p][j+q] * masky[p+centy][q+centx];
                }
            }

            outpic1x[i][j] = xsum;
            outpic1y[i][j] = ysum;

            convx[i][j] = xsum;
            convy[i][j] = ysum;
        }
    }

    // Calculating the magnitude based on the x and y
    // components of convolutions
    // (maxival is the largest value in the magnitude array,
    //  is used later for scaling)

    maxival = 0;
    for (i = mr; i < PICSIZE-mr; i++)
    {
        for (j = mr; j < PICSIZE-mr; j++)
        {
            magnitude[i][j] = sqrt((double) ((outpic1x[i][j]*outpic1x[i][j]) + (outpic1y[i][j]*outpic1y[i][j])));
            
            // find max value iteratively
            if (magnitude[i][j] > maxival)
            {
                maxival = magnitude[i][j];
            }
        }
    }

    // Output magnitude file to PGM image
    for (i = 0; i < PICSIZE; i++)
    {
        for (j = 0; j < PICSIZE; j++)
        {
            magnitude[i][j] = (magnitude[i][j] / maxival) * 255;
            fprintf(fo1, "%c", (char)((int)(magnitude[i][j])));
        }
    }

    // ** Part 2: ** //

    // Find and mark the peaks in the magnitude array
    // NOTE: we have taken the tangent of these values,
    //       to avoid the arctan calculation typically required
    for (i = mr; i < PICSIZE - mr; i++)
    {
        for (j = mr; j < PICSIZE - mr; j++)
        {
            if ( (convx[i][j]) == 0.0)
            {
                convx[i][j] = .00001;
            }

            slope = convy[i][j] / convx[i][j];

            if ( (slope <= .4142) && (slope > -.4142) )
            {
                if ( (magnitude[i][j] > magnitude[i][j-1]) && (magnitude[i][j] > magnitude[i][j+1]))
                {
                    peaks[i][j] = 255;
                }
            }

            else if ( (slope <= 2.4142) && (slope > .4142))
            {
                if ( (magnitude[i][j] > magnitude[i-1][j-1]) && (magnitude[i][j] > magnitude[i+1][j+1]))
                {
                    peaks[i][j] = 255;
                }
            }

            else if ( (slope <= -.4142) && (slope > -2.4142))
            {
                if ( (magnitude[i][j] > magnitude[i+1][j-1]) && (magnitude[i][j] > magnitude[i-1][j+1]))
                {
                    peaks[i][j] = 255;
                }
            }

            else 
            {
                if ( (magnitude[i][j] > magnitude[i-1][j]) && (magnitude[i][j] > magnitude[i+1][j]))
                {
                    peaks[i][j] = 255;
                }
            }
        }
    }

    // Output peaks file
    for (i = 0; i < PICSIZE; i++)
    {
        for (j = 0; j < PICSIZE; j++)
        {
            fprintf(fo2, "%c", (char)((int)peaks[i][j]));
        }
    }

    // ** Part 3 & 4: ** //

    // Init histogram
    for (i = 0; i < PICSIZE; i++)
    {
        histogram[i] = 0;
    }

    // Populate histogram to automatically
    //  calculate hi and lo threshold that
    //  will be used in double threshold 
    for (i = 0; i < PICSIZE; i++)
    {
        for (j = 0; j < PICSIZE; j++)
        {
            (histogram[(int)magnitude[i][j]])++;
        }
    }

    cutOff = 0.1 * PICSIZE * PICSIZE;
    area = 0;
    for (i = PICSIZE - 1; i >= 0; i--)
    {
        area += histogram[i];
        if (area > cutOff)
        {
            break;
        }
    }

    HI = i;
    LO = LOW_PERCENT * HI;

    // Apply double threshold and create final output
    for (i = 0; i < PICSIZE; i++)
    {
        for (j = 0; j < PICSIZE; j++)
        {
            if (peaks[i][j] == 255)
            {
                if (magnitude[i][j] > HI)
                {
                    peaks[i][j] = 0;
                    final[i][j] = 255;
                }

                else if (magnitude[i][j] < LO)
                {
                    peaks[i][j] = 0;
                    final[i][j] = 0;
                }
            }
        }
    }

    // Pick which values will stay in the final output array
    // based on the values of their neighbors
    boolean_moretodo = 1;

    while (boolean_moretodo == 1)
    {
        boolean_moretodo = 0;

        for (i = 0; i < PICSIZE; i++)
        {
            for (j = 0; j < PICSIZE; j++)
            {
                if (peaks[i][j] == 255)
                {
                    for (p = -1; p <= 1; p++)
                    {
                        for (q = -1; q <= 1; q++)
                        {
                            if (final[i+p][j+q] == 255)
                            {
                                peaks[i][j] = 0;
                                final[i][j] = 255;
                                boolean_moretodo = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    // Outputting final .pgm
    for (i = 0; i < PICSIZE; i++)
    {
        for (j = 0; j < PICSIZE; j++)
        {
            fprintf(fo3, "%c", (char)((int)final[i][j]));
        }
    }

    fclose(fp1);
    fclose(fo1);
    fclose(fo2);
    fclose(fo3);

    return 0;
}

