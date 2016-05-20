#include <stdio.h>
#include <stdlib.h>

#define _FILE_OFFSET_BITS 64
#include "graph.h"


int main(int argc, char *argv[])
{

        char* filename = argv[1];
        Graph g = graph_create(1);

	unsigned int siz = (sizeof filename) / (sizeof *filename);

        if (filename[5] == 'd')
            g = graph_loader_directed(filename);
        else if (filename[5] == 'u')
            g = graph_loader_undirected(filename);
        else {
	    printf("%d",siz);
	    printf(" --------------- ");
            printf("File name must start with 'u' for undirected graph or 'd' for directed graph");
            return 0;
        }

	

        //Graph g = graph_loader_undirected("u_Email-Enron.txt");


        printf("Average path distance: %lf\n", average_path_distance(g));
        printf("Network diameter: %d\n", network_diameter_optimized(g));
        printf("Globale efficiency: %lf\n", globale_efficiency(g));
        printf("Clustering coefficient at vertex 3: %lf\n", clustering_coefficient(g, 3));
        printf("Transitivity: %lf\n", transitivity(g));
        printf("Betweenness of vertex 3: %lf\n", betweenness_centrality(g, 3));
        printf("Closeness of vertex 3: %lf\n", closeness_centrality(g, 3));
        printf("Degree Distribution of vertex 3: %lf\n", degree_distribution(g, 3));
        printf("Pearson correlation coefficient: %lf\n", pearson_correlation_coefficient(g));

}
