#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <assert.h>
#define LARGE_NUMBER 1e9



typedef struct edge 
{
	int id;
	int i;
	int j;
	double cost;

	//Indicates if the edge belongs to the spanning tree
	int is_in_tree;


} edge;

typedef struct node 
{
	int id;

	//Indicates if the node already belongs to the spanning tree
	int is_in_tree;

	//number of incident edges
	int degree;
	//Array of pointers to edges: each one points to an incident edge
	edge** edges;

	//temp index: indicates the current number of neighbors stored
	int tmp_degree;



	//DFS 
	//father indicates the pointer on the edge that goes to the father in the DFS tree
	edge* father;
	int visited;
	int beg;
	int end;


	//Dijkstra
	double d;
	int in_S;//0: untouched, 1: in S, 2: permanently labeled



} node;



typedef struct graph
{
	int n;
	int m;

	edge* edges;
	node* nodes;

	//maximum edge cost
	double max_edge_cost;


	//DFS parameter: current dfs "period"
	int t;


} graph;

typedef struct tree 
{
	//The size of the tree is n-1 edges
	int n;
	int * parity; //Christophides : 0 if even 1 if odd 0 is default, indice is the id of the node

	//array of pointers to the edges in the tree
	edge** edges;

	//Tree total cost
	double cost;
} tree;


int readTP1Instance(graph* Gptr, char* instanceFileName);
int initialize_tree(tree* Tptr,int n);
int prim(graph* Gptr,int origin_node, tree*Tptr);
int display_graph(graph*Gptr);
int display_tree(tree* Tptr);
void display_cycle(int cycle[], int n);
void display_edges(edge** e,int m);
int free_graph(graph*Gptr);
int free_tree(tree*Tptr);
int search_id_node(int i, graph* Gptr,int n);
int DFS_explore(graph*Gptr, int s);
int DFS(graph*Gptr,int s,tree*Tptr);
int select_min_d_node(graph* Gptr, node**nptr);
int check_dijkstra_sanity(graph*Gptr);
int dijkstra(graph*Gptr,int s,tree*Tptr);
int edge_search(int i, int j, graph * Gptr);
void matching(graph* Gptr, edge** M);
void parite(tree * Tptr);
void Complete_construction(tree * Tptr, graph * Gptr,graph * Cptr);
void eulerian_multigraph(edge** M, int len, tree * Tptr, graph * Eptr);
int inList(int x, int list[],int n);
int listNodes(edge** edges,int m, int nodes[]);
void createGraph(edge** edges,int m,graph* Gptr);
void copyEdges(edge* destE,edge* sourceE);
int isABridge(edge** edges, int m, int n, int indE,int nextNode);
void supprEdges(edge** e,int m, int k);
void eulerianCycle (graph* Gptr, int cycle[]);
void supprRep(int list[], int n, int sortie[]);
void saveFile(graph* Gptr,int sortie[]);
void Christophides(graph * Gptr);


