#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

int *randomWalk(int N, int step);
int **histogram(int N, int *positions, int *diff);
double *distribution(int N, double *x, int points);

double pi;

int main(int argc, char const *argv[])
{
    srand(time(NULL));

    pi = acos(-1);

    /*
        primer punto
    */

    int i, diff, *rw, rep = 10000, N = 10000;
    int *last = calloc(rep, sizeof(int));

    for(i = 0; i < rep; i++)
    {
        rw = randomWalk(N, 1.0);
        last[i] = rw[N-1];

        free(rw);
    }

    int **hist = histogram(rep, last, &diff);

    double *x = malloc(diff*sizeof(double)), dx = (hist[0][diff - 1] - hist[0][0])/(double) (diff - 1.0);
    for(i = 0; i < diff; i++)
    {
        x[i] = hist[0][0] + i*dx;
    }

    double *dist = distribution(N, x, diff);

    FILE *file = fopen("hist.dat", "w");
    for(i = 0; i < diff; i++)
    {
        fprintf(file, "%d %d %f %f\n", hist[0][i], hist[1][i], x[i], dist[i]);
    }

    fclose(file);

    free(x);
    free(dist);
    free(last);
    free(hist[0]);
    free(hist[1]);
    free(hist);

    /*
        segundo punto
    */
    int j, k, Nvalues = 100, Nrepetitions = 10000;
    double dn = 5.0/Nvalues, mean, mean2;

    file = fopen("N.dat", "w");

    for(k = 0; k < Nvalues; k++)
    {
        N = pow(10, dn*k);

        mean = 0;
        mean2 = 0;

        for(j = 0; j < Nrepetitions; j++)
        {   
            rw = randomWalk(N, 1.0);

            for(i = 0; i < N; i++)
            {
                mean += rw[i];
                mean2 += rw[i]*rw[i];
            }
            free(rw);
        }
        
        mean *= 1.0/(Nrepetitions * N);
        mean2 *= 1.0/(Nrepetitions * N);

        fprintf(file, "%d %.0f %.0f\n", N, mean, mean2);
    }
    fclose(file);

    return 0;
}

double *distribution(int N, double *x, int points)
{
    int i;
    double factor;
    factor = 1/sqrt(2.0*pi*N*N);


    double *p = malloc(points*sizeof(double));
    for(i = 0; i < points; i++)
    {
        p[i] = factor*exp(-pow(x[i], 2.0)/(2.0*N));
    }

    return p;
}

int *randomWalk(int N, int step)
{
    int i, direction, *x = calloc(N, sizeof(int));

    x[0] = 0;

    for(i = 0; i < N-1; i++)
    {
        if(rand()%2 == 1){direction = step;}
        else{direction = -step;}

        x[i + 1] = x[i] + direction;
    }
    return x;
}

int **histogram(int N, int *positions, int *diff)
{
    int i, j, index, temp;

    for(i = 0; i < N - 1; i++)
    {
        index = i;
        for(j = i + 1; j < N; j++)
        {
            if(positions[j] < positions[index])
            {
                index = j;
            }
        }

        temp = positions[i];
        positions[i] = positions[index];
        positions[index] = temp;
    }

    int rep;
    int **hist = malloc(2*sizeof(int *));

    hist[0] = malloc(sizeof(int));
    hist[1] = malloc(sizeof(int));

    *diff = 0;

    for(i = 0; i < N; i++)
    {
        rep = 1;
        for(j = i + 1; j < N; j++)
        {
            if(positions[i] == positions[j]){rep += 1;}
        }

        if(i == 0)
        {

            hist[0][*diff] = positions[i];
            hist[1][*diff] = rep;
        }

        else if(hist[0][*diff] != positions[i])
        {
            *diff += 1;

            hist[0] = realloc(hist[0], (*diff + 1)*sizeof(int));
            hist[1] = realloc(hist[1], (*diff + 1)*sizeof(int));

            hist[0][*diff] = positions[i];
            hist[1][*diff] = rep;
        }
    }

    return hist;
}
