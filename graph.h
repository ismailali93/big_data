#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED

/* basic directed graph type */
typedef struct graph *Graph;

/* create a new graph with n vertices labeled 0..n-1 and no edges */
Graph graph_create(int);

/* free all space used by graph */
void graph_destroy(Graph);

/* add an edge to an existing graph */
/* doing this more than once may have unpredictable results */
void graph_add_edge(Graph, int, int);

Graph graph_loader_directed(char* adresse);

Graph graph_loader_undirected(char* adresse);

int check_edge(Graph,int,int);

int get_vernum(Graph);

int get_edgenum(Graph);

int get_degree(Graph, int);

void change_edgenum(Graph, int);

void print_graph(Graph);

int shortest_path(Graph, int, int);

int sum_shortests_path(Graph, int);

int triangle_num(Graph, int);

double average_path_distance(Graph);

int network_diameter(Graph);

int network_diameter_optimized(Graph);

double clustering_coefficient(Graph, int);

double globale_efficiency(Graph);

double transitivity(Graph);

double betweenness_centrality(Graph, int);

double closeness_centrality(Graph, int);

double PageRank(Graph);

double degree_distribution(Graph, int);

double mean_degree_neighbors(Graph, int);

double pearson_correlation_coefficient(Graph);

#endif // GRAPH_H_INCLUDED
