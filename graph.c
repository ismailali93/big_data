#include <stdlib.h>
#include <stdio.h>
#include "graph.h"

#define _FILE_OFFSET_BITS 64
#define MAX 2147483647

/* structure of a graph
*/
struct graph{
    int ver_num;
    int edge_num;
    struct successor{
        int d;      /* degree */
        int len;    /* number of slots in array */
        int list[1];/* actual list of successors */
    } *alist[1];
};

/* This function is to allocate memory to a new graph
*/
Graph graph_create(int n){
    Graph g;
    int i;

    //allocates memory for the graph structure
    g = malloc(sizeof(struct graph) + sizeof(struct successor *) * (n-1));

    //initializes values
    g->ver_num = n;
    g->edge_num = 0;

    //allocates and intializes the vector of successors
    for(i = 0; i < n; i++){
        g->alist[i] = malloc(sizeof(struct successor));
        g->alist[i]->d = 0;
        g->alist[i]->len = 1;
    }

    return g;
}

/* This function is to free the memory for a graph
*/
void graph_destroy(Graph g){
    int i;
    for(i = 0; i < g->ver_num; i++)
        free(g->alist[i]);
    free(g);
}

/* This function is to add a new edge to a graph
*/
void graph_add_edge(Graph g, int u, int v){
    //reallocates memory if vector of successors is full
	printf("%d",g->alist[u]->d);
	int a = g->alist[u]->d;
	int b = g->alist[u]->len;
    while(g->alist[u]->d >= g->alist[u]->len){
	
        g->alist[u]->len *= 2;
        g->alist[u] = realloc(g->alist[u],sizeof(struct successor) + sizeof(int) * (g->alist[u]->len - 1));
    }
    /* now add edge */
    g->alist[u]->list[g->alist[u]->d++] = v;
    g->edge_num++;
}

/* Load a directed graph from a file
 * @author Antoine
 */
Graph graph_loader_directed(char* adresse)
{
    
    FILE* fichier = NULL;
    Graph g;
    fichier = fopen(adresse, "r");
    
    int i = 0;
    int imax = 0;
    int j = 0;
    int max = 0;
    int temp = 0;

    if (fichier != NULL)
    {
        char actuel;

        printf("\ncreation du graphe g ...");

        do
        {
            fscanf(fichier, "%d", &temp);
            actuel = fgetc(fichier);
            if(temp>max) { max = temp; }
            i++;
            if(i>imax) { imax = i; }
            if(actuel=='\n' || actuel=='\r') {
                i = 0;
                j++;
            }
        } while(actuel != EOF);

        g = graph_create(max+1);

        rewind(fichier);

        int t1 = 0;
        int t2 = 0;

        for(i = 0; i < j; i++) {
            fscanf(fichier, "%d %d\n", &t1, &t2);
	    if(t1 < 0 || t2 < 0) { printf("Neg values : %d --- %d",t1,t2); return; }
	    else graph_add_edge(g,t1,t2);
        }

        printf("\nGraphe cree avec succes\nnombre de sommets : %d\nnombre d'arretes : %d\n\n",g->ver_num,g->edge_num);

        fclose(fichier);
    } else {
        printf("Erreur : fichier introuvale");
    }

    return g;
}

/* Load an undirected graph from a file
 * @author Antoine
 */
Graph graph_loader_undirected(char* adresse)
{
    FILE* fichier = NULL;
    Graph g;

    fichier = fopen(adresse, "r");

    int i = 0;
    int imax = 0;
    int j = 0;
    int max = 0;
    int temp = 0;

    if (fichier != NULL)
    {
        char actuel;

        printf("\ncreation du graphe g ...");

        do
        {
            fscanf(fichier, "%d", &temp);
            actuel = fgetc(fichier);
            if(temp>max) { max = temp; }
            i++;
            if(i>imax) { imax = i; }
            if(actuel=='\n'|| actuel=='\r') {
                i = 0;
                j++;
            }
        } while(actuel != EOF);

        g = graph_create(max+1);

        rewind(fichier);

        int t1 = 0;
        int t2 = 0;

        for(i = 0; i < j; i++) {
            fscanf(fichier, "%d %d\n", &t1, &t2);
            graph_add_edge(g,t1,t2);
            graph_add_edge(g,t2,t1);
        }

        printf("\nGraphe cree avec succes\nnombre de sommets : %d\nnombre d'arretes : %d\n\n",g->ver_num,g->edge_num/2);

        fclose(fichier);
    } else {
        printf("Erreur : fichier introuvale");
    }

    return g;
}

/* this function takes in a graph and
 * returns the number of vertexes
 * @author: Claudionor
 */
int get_vernum(Graph g){
    return g->ver_num;
}

/* this function takes in a graph and
 * returns the number of edges
 * @author: Phuoc
 */
int get_edgenum(Graph g){
    return g->edge_num;
}

/* this function takes in a graph and a vertex and
 * return the degree of that vertex
 * @author: Phuoc
 */
int get_degree(Graph g, int i){
    return g->alist[i]->d;
}

/* this function adds to the total number of edges
 * @author: Phuoc
 */
void change_edgenum(Graph g, int i){
    g->edge_num += i;
}

/* This function prints out the graph in adjacency list
 * @author Ta Duy
  */
void print_graph(Graph g){
    int n = g->ver_num;
    int i,j;
    for (i = 0; i < n; i++){
        printf("%d: ", i);
        for (j = 0; j < g->alist[i]->d; j++)
            printf("%d ", g->alist[i]->list[j]);
        printf("\n");
    }
    printf("\n");
    printf("total: %d edges\n", g->edge_num);
}

/* this function takes in 2 nodes and
 * returns the length of the shortest path between the 2 nodes
 * @author Ta Duy
 * modified by Claudionor
 */
int shortest_path(Graph g, int i, int j) {
    //declaration of a queue
    int *queue;
    queue = malloc(sizeof(int) * g->ver_num);
    int queue_size;
    int index = 0;

    //answer
    int sp;

    //vector of the distances from vertex i initialized with max int value
    int *dist = malloc(sizeof(int) * g->ver_num);
    int h;
    for (h = 0; h < g->ver_num; h++)
        dist[h] = MAX;

    //adds vertex i to the queue and updates his distance from himself to 0
    queue[0] = i;
    dist[i] = 0;
    queue_size=1;

    //updates the distances while there is a vertex in the queue using dijkstra's algorithm
    //stops when find the target vertex
    while (index < queue_size){
        if (queue[index] == j)
            break;
        else{
            int s = queue[index];
            int m;
            for (m = 0; m < g->alist[s]->d; m++)
                if (dist[g->alist[s]->list[m]] == MAX){
                    queue[queue_size] = g->alist[s]->list[m];
                    queue_size++;
                    dist[g->alist[s]->list[m]] = dist[s] + 1;
                }
            index++;
        }
    }

    //saves answer and free memory
    sp=dist[j];
    free(queue);
    free(dist);

    return sp;
}

/* this function takes in node and runs Dijkstra's algorithm
 * to return the shortest paths from it to the other nodes
 * @author Ta Duy
 * modified by Claudionor
 */
int* shortests_paths(Graph g, int i) {
    int n=g->ver_num;
    //declaration of a queue
    int *queue = malloc(sizeof(int) * n);
    int queue_size;
    int index = 0;

    //vector of the distances from vertex i initialized with max int value
    int *dist = malloc(sizeof(int) * n);
    int h;
    for (h = 0; h < n; h++)
        dist[h] = MAX;

    //adds vertex i to the queue and updates his distance from himself to 0
    queue[0] = i;
    dist[i] = 0;
    queue_size=1;

    //runs dijkstra's algorithm
    while (index < queue_size){
        int s = queue[index];
        for (h = 0; h < g->alist[s]->d; h++){
            int aux=g->alist[s]->list[h];
            if (dist[aux] == MAX){
                queue[queue_size] = aux;
                queue_size++;
                dist[aux] = dist[s] + 1;
            }
        }
        index++;
    }

    //free memory
    free(queue);

    return dist;
}

/* This function  is to check whether an edge is present in the graph
 * @author Ta Duy
 */
int check_edge(Graph g, int u, int v) {
    int i;
    for (i = 0; i < g->alist[u]->d; i++)
        if (v == g->alist[u]->list[i])
            return 1; //linear search in the vector os successors of u
    return 0;
}

/* This function  is to calculate the number of triangles
 * with a vertex at node i
 * @author Ta Duy
 * modified by Claudionor
 */
int triangle_num (Graph g, int i) {
    int u, v;
    int s = 0;
    //uses the check_edge function to count how many neighbors of u are also neighbors of one another
    for (u = 0; u < g->alist[i]->d; u++)
        if(check_edge(g,g->alist[i]->list[u],i)==1)
            for (v = u + 1; v < g->alist[i]->d; v++)
                if(check_edge(g,g->alist[i]->list[v],i)==1)
                        s += check_edge(g, g->alist[i]->list[u],g->alist[i]->list[v])*
                            check_edge(g, g->alist[i]->list[v],g->alist[i]->list[u]);
    return s;

}

/* This functions returns the average path distance in the graph
 * @author Ta Duy and Claudionor
 */
double average_path_distance(Graph g) {
    int n = g -> ver_num;
    int d = 0;
    int i, j;
    int * dist;
    //sums the distance from each vertex to all other vertexes
    for (i = 0; i < n; i++){
        dist=shortests_paths(g,i);
        for(j=0;j<n;j++){
            if(dist[j]==MAX)
                return MAX;
            d+=dist[j];
        }
        //free memory
        free(dist);
    }
    //divides the sum of the distances by the number of pair of vertexes
    return (double)d/(n*(n-1));
}

/* This functions returns the global efficiency in the graph
 * @author Ta Duy and Claudionor
 */
double globale_efficiency(Graph g) {
    int n = g -> ver_num;
    double d = 0;
    int i, j;
    int* dist;
    //sums the inverse of the distance from each vertex to all others
    for (i = 0; i < n; i++){
        dist=shortests_paths(g,i);
        for (j = 0; j < n; j++)
            if(j!=i)
                d+=1.0/dist[j];
        //free memory
        free(dist);
    }
    //divides the sum of the inverses by the number of pair of vertexes
    return d/(n * (n - 1));
}

/* This functions returns the network diameter in the graph
 * @author Claudionor
 */
int network_diameter(Graph g) {
    int max=-1;
    int n=g -> ver_num;
    int i, j;
    int* dist;
    //takes the maximum of all the distances from one vertex to another one
    for (i = 0; i < n; i++){
        dist=shortests_paths(g,i);
        for (j=0; j < n; j++)
            if(dist[j]>max)
                max=dist[j];
        //free memory
        free(dist);
    }

    return max;
}
/* This functions returns the network diameter in the graph
 * optimized version
 * @author Antoine
 */
int network_diameter_optimized(Graph g) {

    // Mise en place d'une structure de file dont l'\E9l\E9ment
    // de base est un triplet (sommet 1, sommet 2, distance
    //entre les sommets 1 et 2)

    typedef struct Element Element;
    struct Element
    {
        int sommet_1;
        int sommet_2;
        int dist;
        Element *suivant;
    };

    typedef struct File File;
    struct File
    {
        Element *premier;
        Element *dernier;
    };

    File *initialisation(i)
    {
        File *file = malloc(sizeof(*file));
        Element *element = malloc(sizeof(*element));

        if (file == NULL || element == NULL)
        {
            exit(EXIT_FAILURE);
        }

        element->sommet_1 = 0;
        element->sommet_2 = 0;
        element->dist = 0;
        element->suivant = NULL;
        file->premier = element;
        file->dernier = element;
        return file;
    }

    void enfiler(File *file, int s_1, int s_2, int d)
    {
        Element *nouveau = malloc(sizeof(*nouveau));
        if (file == NULL || nouveau == NULL)
        {
            exit(EXIT_FAILURE);
        }

        nouveau->sommet_1 = s_1;
        nouveau->sommet_2 = s_2;
        nouveau->dist = d;
        nouveau->suivant = NULL;
        file->dernier->suivant = nouveau;
        file->dernier = nouveau;
    }

    void defiler(File *file)
    {
        if (file == NULL)
        {
            exit(EXIT_FAILURE);
        }

        if (file != NULL && file->premier != NULL)
        {
            Element *elementDefile = file->premier;
            file->premier = file->premier->suivant;
            free(elementDefile);
        }

        return 0;
    }

    int n = g->ver_num; // nombre de sommets du graphe
    int i,j; // variables de boucle
    int max = 0; // Variable portant la r\E9ponse

    File *file = initialisation(); // Initialisation de la file

    int** distances = (int**) malloc(n * n * sizeof(int*));
    for (i = 0; i < n; i++) { distances[i] = (int*) malloc(n * sizeof(int)); }

    int** traite = (int**) malloc(n * n * sizeof(int*));
    for (i = 0; i < n; i++) { traite[i] = (int*) malloc(n * sizeof(int)); }

    int** ajoute = (int**) malloc(n * n * sizeof(int*));
    for (i = 0; i < n; i++) { ajoute[i] = (int*) malloc(n * sizeof(int)); }

    //int distances[n][n]; // Matrice des distances.
    //int traite[n][n];   // Matrice de bool\E9en pr\E9sentant si le couple (i,j) a \E9t\E9 trait\E9
    //int ajoute[n][n];   // Matrice de bool\E9en pr\E9sentant si le couple (i,j) a \E9t\E9 d\E9j\E0 \E9t\E9 ajout\E9 \E0 la pile

    // Initialisation des variables

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            traite[i][j] = 0; // Le couple (i,j) n'a pas encore \E9t\E9 trait\E9
            ajoute[i][j] = 0; // Le couple (i,j) n'a pas encore \E9t\E9 ajout\E9 \E0 la pile
            distances[i][j] = -1; // La distance entre i et j est initialis\E9e \E0 -1

        }
        distances[i][i] = 0; // La distance entre i et i est initialis\E9e \E0 0
        traite[i][i] = 1;    // Le couple (i,i) a \E9t\E9 trait\E9
        ajoute[i][i] = 1;    // Le couple (i,i) a \E9t\E9 ajout\E9 \E0 la pile
        enfiler(file,i,i,0); // Le couple (i,i) est ajout\E9 \E0 la pile

    }

    int temp_i, temp_j, temp_d;
    int sommet_etudie;
    int n1;

    while(file->premier != NULL) { // Tant que la file d'attente n'est pas vide
        temp_i = file->premier->sommet_1;
        temp_j = file->premier->sommet_2;
        temp_d = file->premier->dist;
        defiler(file);  // On extrait le triplet (sommet i, sommet j, distance d de i vers j)

        if(traite[temp_i][temp_j] == 0) {
            distances[temp_i][temp_j] = temp_d;
            traite[temp_i][temp_j] = 1;
        }

        // Pour chaque voisin k du sommet temp_j, on consid\E8re le chemin de temp_i vers k

        n1 = g->alist[temp_j]->d;
        for(i = 0; i < n1; i++) {

            sommet_etudie = g->alist[temp_j]->list[i];

            if (ajoute[temp_i][sommet_etudie] == 0) {
                ajoute[temp_i][sommet_etudie] = 1;
                enfiler(file,temp_i,sommet_etudie,temp_d+1);
            }
        }
    }

    int convexe = 0;

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            if (distances[i][j] > max) { max = distances[i][j]; }
            if (distances[i][j] == -1 && distances[j][i] == -1) { convexe = 1; }
        }
    }

    if (convexe == 1) {max = -1; }

    return max;
}

/* This function returns the clustering coefficient
 * @author Ta Duy
 */
double clustering_coefficient(Graph g, int i) {
    int E = triangle_num(g, i);
    int a = g->alist[i]->d;
    //if there are sufficient edges, return the formula described in the paper
    if (a > 1)
        return (double)E*2.0/(a * (a - 1));
    return 0;
}

/* This function is to calculate the transitivity of the graph
 * @author Ta Duy
 */
double transitivity(Graph g) {
    int i, s = 0, n = g->ver_num;
    //sums the total of triangles for the graph
    for (i = 0; i < n; i++)
        s += triangle_num(g, i);
    //returns the formula given in paper for graph transitivity
    double d=s;
    d/=n;
    d/=n-1;
    d/=n-2;
    return d;
}

/* This function is to find the number of shortest paths
 * between 2 nodes
 * @author Ta Duy
 * modified: Claudionor
 */

int number_of_shortest_path(Graph g, int i, int j) {
    int n=g->ver_num;

    //declares a queue
    int *queue;
    queue = malloc(sizeof(int) * n);
    int queue_size = 0;
    int index = 0;

    //declares a vector of the amount of shortest paths to each vertex
    int *nsp;
    nsp = malloc(sizeof(int) * n);
    int h;
    for (h = 0; h < n; h++)
        nsp[h] = 0;
    nsp[i]=1;

    //declares a vector of distances from i and initializes all distances to INT_MAX
    int *dist = malloc(sizeof(int) * n);
    for (h = 0; h < n; h++)
        dist[h] = MAX;

    //adds i to the queue and makes the distance to himself equal to 0
    queue[0] = i;
    queue_size++;
    dist[i] = 0;

    //while running dijkstra's algorithm, counts the amount of shortests paths to each vertex
    //and uses it to calculate the number of shortests paths to the next one
    while (index < queue_size && dist[queue[index]] != dist[j]) {
        int s = queue[index];
        int m;
        for (m = 0; m < g->alist[s]->d; m++){
            int aux=g->alist[s]->list[m];
            if (dist[aux] == MAX){
                queue[queue_size] = aux;
                queue_size++;
                dist[aux] = dist[s] + 1;
                nsp[aux] = nsp[s];
            }
            else if(dist[aux]==dist[s]+1)
                nsp[aux]+=nsp[s];
        }
        index++;
    }

    //saves answer and frees memory
    int ans=nsp[j];
    free(queue);
    free(nsp);
    free(dist);

    return ans;

}

/* This function is to calculate the betweenness centrality of a node
 * @author Ta Duy and Claudionor
 */
double betweenness_centrality(Graph g, int x) {
    int i, j, n = g->ver_num;
    double s = 0;
    int* distx=shortests_paths(g,x);
    int* disti;
    //calculates the number os shortests paths from i to j that passes by x
    //and divides by the number of all shortests paths from i to j, to all pairs i,j
    for (i = 0; i < n; i++){
        disti=shortests_paths(g,i);
        for (j = i + 1; j < n; j++)
            if(disti[j] < MAX && disti[j] == disti[x] + distx[j])
                s += number_of_shortest_path(g, i, x) *
                number_of_shortest_path(g, x, j) *
                1.0/number_of_shortest_path(g, i, j);
        free(disti);
    }
    free(distx);
    return s;
}

/* This function is to calculate the closeness centrality
 * @author Ta Duy and Claudionor
 */
double closeness_centrality(Graph g, int x) {
    int i, s = 0;
    int n=g->ver_num;
    int* distx=shortests_paths(g,x);
    //calculates the sum of the distances from x to all other vertexes
    for (i = 0; i < n ; i++)
        s += distx[i];
    free(distx);
    //returns the inverse of the average path length from x
    return (double)n/s;
}

/* This function calculates the fraction of each degree
 * @author Claudionor
 */
double degree_distribution(Graph g, int degree){
    int n=g->ver_num;
    int cont=0;
    int i;
    //counts how many vertexes have degree degree
    for(i=0;i<n;i++)
        if(g->alist[i]->d==degree)
            cont++;
    //returns the proportion of the vertexes of that degree
    return (double)cont/n;
}

/* This function is to calculate the Pearson correlation coefficient
 * @author Claudionor
 */
double pearson_correlation_coefficient(Graph g){
    int m=g->edge_num;
    int n=g->ver_num;
    int sum1=0;
    int sum2=0;
    int sum3=0;
    double d1,d2,d3;
    int i,j;
    //calculates intermidiate sums used in the formula described in the paper
    for(i=0;i<n;i++)
        for(j=i+1;j<n;j++)
            if(check_edge(g,i,j)==1){
                sum1+=g->alist[i]->d*g->alist[j]->d;
                sum2+=g->alist[i]->d+g->alist[j]->d;
                sum3+=(g->alist[i]->d)*(g->alist[i]->d)+(g->alist[j]->d)*(g->alist[j]->d);
            }
    //multiplies the sums by their constants
    d1=(double)sum1/m;
    d2=(double)sum2/(2.0*m);
    d3=(double)sum3/m;
    //calculates the final formula
    return (d1-d2*d2)/(d3-d2*d2);
}
