#include "TP1Functions.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include<stdio.h>

int readTP1Instance(graph* Gptr, char* instanceFileName){
	int rval = 0;

	//We open the file
	FILE* fin = fopen(instanceFileName,"r");

	//We read the number of nodes and vertices
	int n,m;
	int e,i,j;

	rval = fscanf(fin,"%d,%d\n", &n,&m);
	Gptr->n = n;
	Gptr->m = m;

	//We malloc and fill the nodes in their structure
	Gptr->nodes= (node*)malloc(sizeof(node)*n);
	node* iptr;
	for( i = 0 ;  i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		iptr->id = i;
		//The node is initially not in the tree
		iptr->is_in_tree = 0;
		//The node has initially degree zero
		iptr->degree= 0;
		iptr->tmp_degree= 0;
	}

	//We malloc, read and fill the edges in their structure
	Gptr->edges = (edge*)malloc(sizeof(edge)*m);
	double cost;
	edge* eptr;
	double max_edge_cost = -LARGE_NUMBER;
	for( e = 0 ;  e < m ; e++)
	{
		eptr = &(Gptr->edges[e]);
		rval=fscanf(fin,"%d,%d,%lf\n", &i,&j,&cost);
		eptr->id = e;
		eptr->i = i;
		eptr->j = j;
		eptr->cost = cost;
		//We update the largest edge cost
		if(max_edge_cost < cost)
			max_edge_cost = cost;
		//The edge is initially not in the tree
		eptr->is_in_tree = 0;

		//We update the degree of the edge's endpoints
		//tail node
		iptr = &(Gptr->nodes[i]);
		iptr->degree++;
		//head node
		iptr = &(Gptr->nodes[j]);
		iptr->degree++;

		//fprintf(stderr,"%d,%d,%d,%lf\n",eptr->ideptr->i,eptr->j,eptr->cost);
	}
	//We store the largest edge cost in the graph structure
	Gptr->max_edge_cost = max_edge_cost;


	//We close the file
	fclose(fin);

	//We fill in the neighbors for each node:
	//Note that we are working with undirected graphs
	//We malloc their arrays
	int degree;
	for( i = 0 ; i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		degree = iptr->degree;
		iptr->edges = (edge**)malloc(sizeof(edge*)*degree);
	}
	//We loop over the edges to fill the neighbors arrays
	int tmp;
	for( e = 0 ;  e < m ; e++)
	{
		eptr = &(Gptr->edges[e]);
		i = eptr->i;
		j = eptr->j;

		//tail node
		iptr = &(Gptr->nodes[i]);
		tmp = iptr->tmp_degree;
		iptr->edges[tmp] = eptr;
		iptr->tmp_degree++;

		//head node
		iptr = &(Gptr->nodes[j]);
		tmp = iptr->tmp_degree;
		iptr->edges[tmp] = eptr;
		iptr->tmp_degree++;
	}
	return rval;
}

int initialize_tree(tree* Tptr,int n){
	int rval = 0;

	Tptr->n = n;
	Tptr->parity = (int *)calloc(n, sizeof(int));
	Tptr->edges = (edge**)malloc(sizeof(edge*) * (n-1));
	Tptr->cost = 0;

	return rval;
}

int prim(graph* Gptr, int origin_node, tree*Tptr){
	int rval = 0;

	int m = Gptr->m;
	int n = Gptr->n;
	int e,i,j;

	//We reinitialize the edges and nodes to not be in the tree
	for( e = 0 ;  e < m ; e++)
		Gptr->edges[e].is_in_tree = 0;
	for( i = 0 ;  i < n ; i++)
		Gptr->nodes[i].is_in_tree = 0;
	//We mark the origin node as in the tree
	Gptr->nodes[origin_node].is_in_tree = 1;

	//Prim's main body
	double mincost;
	double maxcost = Gptr->max_edge_cost;
	edge*eptr;
	node*iptr,*jptr;
	int candidate_e;
	int tree_size = 0;
	for(tree_size = 0 ; tree_size < n-1 ; tree_size++)
	{
		mincost = maxcost; 
		candidate_e = -1;

		for( e = 0 ;  e < m ; e++)
		{
			//We fetch the pointer to the edge
			eptr = &(Gptr->edges[e]);

			//If it is already in the tree, skip it
			if(eptr->is_in_tree)
				continue;

			//We fetch the pointers to the endpoint nodes
			i = eptr->i;
			iptr = &(Gptr->nodes[i]);
			j = eptr->j;
			jptr = &(Gptr->nodes[j]);

			//If i and j are in the tree, skip it
			if(	(iptr->is_in_tree) && (jptr->is_in_tree))
				continue;

			//If i and j are both outside the tree skip it
			if(	!(iptr->is_in_tree) && !(jptr->is_in_tree))
				continue;

			//Otherwise, we consider it and check its cost
			if(eptr->cost < mincost)
			{
				mincost = eptr->cost;	
				candidate_e = eptr->id;
			}
		}

		//If we did not find a candidate, the graph is not connected 
		if(candidate_e == -1)
		{
			fprintf(stderr,"The graph is disconnected. Aborting Prim\n");
			rval = 1;
		}
		//Otherwise, we add the least cost edge we found to the tree
		//We fetch the pointer to the edge
		eptr = &(Gptr->edges[candidate_e]);
		//We mark it as in the tree
		eptr->is_in_tree = 1;

		//We fetch the pointers to the endpoint nodes
		i = eptr->i;
		iptr = &(Gptr->nodes[i]);
		j = eptr->j;
		jptr = &(Gptr->nodes[j]);
		//And we mark both nodes as marked (one of them already is)
		iptr->is_in_tree = 1;
		jptr->is_in_tree = 1;
	}	

	//We fill the tree 
	Tptr->cost = 0;
	int count = 0;
	//We fill the tree structure
	for( e = 0 ; e < m ; e++)
	{
		eptr = &(Gptr->edges[e]);
		//If the edge is not marked, we skip it
		if(!eptr->is_in_tree)
			continue;
		//Otherwise, we put it in the tree
		Tptr->edges[count++] = eptr;
		Tptr->cost += eptr->cost;
	}
	return rval;
}

int display_graph(graph*Gptr){
	int rval = 0; 
	int i,e,n,m;

	n = Gptr->n;
	m = Gptr->m;
	fprintf(stderr,"The graph has %d nodes and %d edges:\n",n,m);
	fprintf(stderr,"Displaying nodes:\n");

	node*iptr;
	edge*eptr;

	// printf("0 : %d nb arêtes%d\n",Gptr->nodes[0].id,Gptr->nodes[0].degree);
	// printf("1 : %d nb arêtes%d\n",Gptr->nodes[1].id,Gptr->nodes[1].degree);
	// printf("2 : %d\n",Gptr->nodes[2].id);

	for( i = 0 ; i < n ; i++){
		iptr = &(Gptr->nodes[i]);
		for( e = 0 ; e < iptr->degree ; e++){
			eptr = iptr->edges[e];
			fprintf(stderr,"\t\t(%d,%d)\t[w=%lf]\n",eptr->i,eptr->j,eptr->cost);
		}
	}	
	// printf("End Display\n");
	return rval;
}

int display_tree(tree* Tptr){
	int rval = 0;
	int n = Tptr->n;
	int e;
	edge*eptr;

	fprintf(stderr,"Tree\n");
	fprintf(stderr,"WEIGHT:\t%lf\n",Tptr->cost);
	for( e = 0 ; e < n-1 ; e++)
	{
		eptr = Tptr->edges[e];
		fprintf(stderr,	"Edge #%d (%d,%d)\tcost:\t%lf\n",
				eptr->id,eptr->i, eptr->j, eptr->cost);	
	}
	fprintf(stderr,"\n\n\n");
	return rval;
}

void display_cycle(int cycle[], int n){
	for(int i=0; i<n;i++){
		printf("%d ",cycle[i]);
	}
	printf("\n");
}

void display_edges(edge** e,int m){
	edge* eptr;
	printf("Nombre d'arêtes : %d\n",m);
	for(int i=0;i<m;i++){
		eptr=e[i];
		printf("(%d,%d)\n",eptr->i,eptr->j);
	}
}

int free_graph(graph*Gptr)
{
	int rval=0;
	int n = Gptr->n;
	int i;
	node* iptr;
	for( i = 0 ; i < n ; i++){
		iptr = &(Gptr->nodes[i]);
		free(iptr->edges);
	}
	free(Gptr->edges);
	free(Gptr->nodes);

	return rval;
}

int free_tree(tree*Tptr){
	int rval=0;
	free(Tptr->edges);
	return rval;
}

int search_id_node(int i, graph* Gptr,int n){
	// Search the new index of node of id i
	int k=0;
	while((k<n)&&(Gptr->nodes[k].id!=i)){
		k++;
	}
	return k;
}

int DFS_explore(graph*Gptr, int s){
	// printf("In DFS_explore\n");
	int rval = 0;
	
	node*iptr = &(Gptr->nodes[search_id_node(s,Gptr,Gptr->n)]);
	iptr->visited = 1;
	Gptr->t++;
	iptr->beg=Gptr->t;
	int degree = iptr->degree;
	int e;
	edge*eptr;
	int i,j;
	node*jptr;

	for( e  = 0 ;  e  < degree ; e++){
		eptr = iptr->edges[e];
		i = eptr->i;
		j = eptr->j;
		//We select the neighboring node
		if(i != s)
			j=i;
		jptr = &(Gptr->nodes[search_id_node(j,Gptr,Gptr->n)]);

		//If neighbor already visited, skip it
		if(jptr->visited)
			continue;
		//Otherwise, we explore it
		jptr->visited = 1;
		jptr->father = eptr;
		// printf("Val de j %d\n",j);
		DFS_explore(Gptr,j);
		
	}
	iptr->visited = 2;
	Gptr->t++;
	iptr->end = Gptr->t;

	// printf("DFS_explore end\n");

	return rval;
}



int DFS(graph*Gptr,int s, tree*Tptr){
	// printf("Help\n");
	int rval = 0;

	//We reinitialize the DFS parameters
	Gptr->t = 0;

	int n = Gptr->n;
	int i;
	node*iptr;
	// Here i is an index
	for( i = 0 ; i < n ; i++){
		iptr = &(Gptr->nodes[i]);
		iptr->father = NULL;
		iptr->visited= 0;
		iptr->beg= -1;
		iptr->end= -1;
	}

	// printf("Hello2\n");
	// Here i is an id (the name of the node)
	i = s;
	int j;
	while(i != -1){
		//We explore the node (of label i but different index than i)
		iptr = &(Gptr->nodes[search_id_node(i,Gptr,Gptr->n)]);
		if(iptr->visited == 0){
			DFS_explore(Gptr,i);
			// printf("DFS_explore success\n");
		}

		//We pick the next node to explore
		i=-1;
		for( j = 0 ; j < n ; j++)
		{
			iptr = &(Gptr->nodes[j]);
			if(iptr->visited)
				continue;
			i=Gptr->nodes[j].id; // i take the label of the node
			break;
		}
	}
	// printf("Hello3\n");

	//We now find the set of edges composing the tree
	int cnt = 0;
	Tptr->cost = 0;
	for( i = 0 ; i < n ; i++){
		iptr = &(Gptr->nodes[i]);
		//If we are at the root
		if(iptr->father == NULL)
			continue;

		Tptr->edges[cnt++]=iptr->father;
		Tptr->cost += iptr->father->cost;
	}

	// printf("Hello4\n");
	// Number of nodes visited in Gptr
	int visited=0;
	for(int i=0;i<n;i++){
		// printf("%d\n",i);
		iptr = &(Gptr->nodes[i]);
		visited = visited+ iptr->visited/2;
	}
	rval=visited;

	return rval;
}


int select_min_d_node(graph* Gptr, node**nptr){
	int rval = 0;

	int i;
	int n = Gptr->n;

	double min_d = LARGE_NUMBER;; 
	node*iptr;

	for( i = 0 ; i < n ; i++){
		iptr = &(Gptr->nodes[i]);
		if(iptr->in_S != 1)
			continue;

		if(iptr->d > min_d)
			continue;

		min_d = iptr->d;
		*nptr = iptr;		
	}
	return rval;
}

int check_dijkstra_sanity(graph*Gptr){
	int rval = 0;

	int m = Gptr->m;
	int e,i,j;
	edge*eptr;
	node*iptr,*jptr;

	int sanity = 1;

	for( e = 0 ; e < m ; e++){
		eptr = &(Gptr->edges[e]);
		i = eptr->i;
		iptr = &(Gptr->nodes[i]);
		j = eptr->j;
		jptr = &(Gptr->nodes[j]);

		if(iptr->d + eptr->cost < jptr->d){
			fprintf(stderr,
					"Dijkstra failed for edge %d between nodes %d and %d: %lf+%lf = %lf < %lf\n",
					e,i,j,iptr->d,eptr->cost,iptr->d + eptr->cost,jptr->d);
			sanity = 0;
		}
		if(jptr->d + eptr->cost < iptr->d){
			fprintf(stderr,
					"Dijkstra failed for edge %d between nodes %d and %d: %lf+%lf = %lf < %lf\n",
					e,j,i,jptr->d,eptr->cost,jptr->d + eptr->cost,iptr->d);
			sanity = 0;
		}
	}
	if(sanity)
		fprintf(stderr,"Potentials are consistent. Disjktra ended successfully.\n");
	else
		fprintf(stderr,"Potentials are NOT consistent. Disjktra failed.\n");

	return rval;
}


int dijkstra(graph*Gptr,int s,tree*Tptr){
	int rval = 0;

	int n = Gptr->n;
	int i,j;
	node*iptr,*jptr;
	for( i = 0 ; i < n ; i++){
		iptr = &(Gptr->nodes[i]);
		iptr->father = NULL;
		iptr->d = LARGE_NUMBER;
		//All nodes are outside S at first
		iptr->in_S = 0;
	}

	i=s;
	iptr = &(Gptr->nodes[i]);
	iptr->father = NULL;
	iptr->d = 0; 
	iptr->in_S = 1;

	edge*eptr;
	int e;
	int m = Gptr->m;
	for( e = 0 ; e < m ; e++){
		eptr = &(Gptr->edges[e]);
		eptr->is_in_tree = 0;
	}

	//Dijkstra loop
	while(1){
		//We select the node with lowest label in S
		iptr = NULL;
		select_min_d_node(Gptr, &iptr);
		if(iptr == NULL)
			break;

		i = iptr->id;

		//We delete it from S
		iptr->in_S = 2;

		//We check all incident edges
		for( e = 0 ; e < iptr->degree ; e++)
		{
			eptr = iptr->edges[e];
			//We catch the index of the neighboring node
			if(eptr->i == i)
				j = eptr->j;
			else
				j = eptr->i;
			//and its pointer
			jptr = &(Gptr->nodes[j]);

			//if already permanently labeled, skip it
			if(jptr->in_S == 2)
				continue;

			//If untouched yet, add it to S
			if(jptr->in_S == 0)
				jptr->in_S =1;

			//If improvable, do it
			if(jptr->d > iptr->d + eptr->cost){
				if(jptr->father != NULL)
					jptr->father->is_in_tree = 0;
				jptr->father = eptr;
				jptr->father->is_in_tree = 1;
				jptr->d = iptr->d + eptr->cost;
			}
		}
	}

	check_dijkstra_sanity(Gptr);

	//We fill the tree 
	Tptr->cost = 0;
	int count = 0;
	//We fill the tree structure
	for( e = 0 ; e < m ; e++){
		eptr = &(Gptr->edges[e]);
		//If the edge is not marked, we skip it
		if(!eptr->is_in_tree)
			continue;
		//Otherwise, we put it in the tree
		Tptr->edges[count++] = eptr;
		Tptr->cost += eptr->cost;
	}

	return rval;
}

int edge_search(int i, int j, graph * Gptr){
	// Search the edge (i,j) in graph Gptr
	int k = 0;
	while((k<Gptr->m)&&(Gptr->edges[k].i != i || Gptr->edges[k].j != j)&&(Gptr->edges[k].i != j || Gptr->edges[k].j != i)){
		k++;
	}
	//printf("Valeur de k %d\n",k);
	return k;
}

void matching(graph* Gptr, edge** M){
	// Calculation of a simple matching in Gptr
	// We add every alternative edge to the matching
	int i=0; edge* eptr; int ind; int incr = 0;
	while(i<Gptr->n){
		ind = edge_search(Gptr->nodes[i].id,Gptr->nodes[i+1].id,Gptr);
		// printf("Indince %d\n",ind);
		eptr = &(Gptr->edges[ind]);
		M[incr]=eptr;
		incr++;
		i +=2;
	}
}

void parite(tree * Tptr){	
	// Check the parity for each node of a tree
	int n = Tptr->n;
	int m = n - 1; // Number of edges
	int node_id; int k;
	for(k = 0; k < m; k++){
		node_id = Tptr->edges[k]->j;
		Tptr->parity[node_id] += 1;
		node_id = Tptr->edges[k]->i;
		Tptr->parity[node_id] += 1;
	}
	for(k = 0; k < n; k++){
		Tptr->parity[k] = Tptr->parity[k] % 2;
	}
}




void Complete_construction(tree * Tptr, graph * Gptr,graph * Cptr){ 
	// Construct a complete graph based on the odd nodes of a tree
	int n = Tptr->n; int i; int num_nodes = 0; //number of nodes

	// Number of uneven nodes
	int length=0;
	for(i = 0; i < n; i++){
		length += Tptr->parity[i];
	}
	int nodes_id[length];

	// List of the uneven nodes
	for(i = 0; i < n; i++){
		if(Tptr->parity[i]){
			nodes_id[num_nodes] = i;
			num_nodes += 1;
		}
	}

	// Nodes construction
	Cptr->n = num_nodes;
	Cptr->nodes = (node *) malloc(sizeof(node) * num_nodes);
	node* iptr;
	for(i = 0; i < num_nodes; i++){
		iptr = &(Cptr->nodes[i]);
		iptr->id = nodes_id[i];
		iptr->degree = num_nodes-1;
		//printf("Noeuds dans le n graphe %d\n",Cptr->nodes[i].id);
	}

	// Edges construction
	Cptr->m = num_nodes * (num_nodes - 1) / 2; // Number of edges for a complete graph
	Cptr->edges = (edge *) malloc(sizeof(edge)*Cptr->m);
	int j; int inc = 0; int id;
	double max_edge_cost = -LARGE_NUMBER;

	for(i = 0; i < num_nodes; i++){
		for(j = i+1; j < num_nodes; j++){
			// Search of the edge (i,j)
			id = edge_search(Cptr->nodes[i].id,Cptr->nodes[j].id, Gptr);
			Cptr->edges[inc].id = inc;
			// printf("Valeur id %d\n",Cptr->edges[inc].id);
			Cptr->edges[inc].i = Cptr->nodes[i].id;
			Cptr->edges[inc].j = Cptr->nodes[j].id;
			//printf("Ajout de l'arête (%d,%d)\n",Cptr->edges[inc].i,Cptr->edges[inc].j);
			Cptr->edges[inc].cost = Gptr->edges[id].cost;
			//printf("Recup cout\n");

			if(max_edge_cost < Cptr->edges[inc].cost){
				max_edge_cost = Cptr->edges[inc].cost;
			}

			Cptr->edges[inc].is_in_tree = 0;

			inc++;
		}
	}
	Cptr->max_edge_cost = max_edge_cost;
	//printf("Récup des arêtes\n");

	int degree;
	for( i = 0 ; i < num_nodes ; i++){
		iptr = &(Cptr->nodes[i]);
		degree = iptr->degree;
		iptr->edges = (edge**)malloc(sizeof(edge*)*degree);
	}
	//printf("Ini liste arêtes des noeuds\n");

	edge* eptr; int e;int idNode;
	int tmp;
	for( e = 0 ;  e < Cptr->m ; e++){
		eptr = &(Cptr->edges[e]);
		i = eptr->i;
		j = eptr->j;
		// printf("Récup arête\n");

		//tail node
		idNode = search_id_node(i,Cptr,num_nodes);
		iptr = &(Cptr->nodes[idNode]);
		tmp = iptr->tmp_degree;
		iptr->edges[tmp] = eptr;
		iptr->tmp_degree++;
		// printf("Ajout queue\n");


		//head node
		// printf("Valeur de j %d\n", j);
		idNode = search_id_node(j,Cptr,num_nodes);
		iptr = &(Cptr->nodes[idNode]);
		tmp = iptr->tmp_degree;
		iptr->edges[tmp] = eptr;
		iptr->tmp_degree++;
		// printf("Ajout tete\n");

	}
	// printf("Fin complete\n");
}

void eulerian_multigraph(edge** M,int len, tree * Tptr, graph * Eptr){
	// Adding the matching to the spanning tree

	// Number of nodes and edges
	Eptr->n = Tptr->n;
	Eptr->m = len+ Tptr->n -1;
	int n = Eptr->n;
	int m = Eptr->m;

	// Fill the nodes
	Eptr->nodes= (node*)malloc(sizeof(node)*n);
	node* iptr;
	for (int i= 0; i<n; i++){
		iptr = &(Eptr->nodes[i]);
		iptr->id = i;
		iptr->is_in_tree = 0;
		iptr->degree= 0;
		iptr->tmp_degree= 0;
	}

	//Fill the edges
	Eptr->edges = (edge*)malloc(sizeof(edge)*m);
	edge* eptr;edge* eptr2;
	double max_edge_cost = -LARGE_NUMBER;

	// Copy of the tree
	int e;
	for(e=0; e<(Tptr->n-1);e++){
		eptr = &(Eptr->edges[e]);
		eptr2 = Tptr->edges[e];
		eptr->id =e;
		eptr->i = eptr2->i;
		eptr->j = eptr2->j;
		eptr->cost = eptr2->cost;

		if(max_edge_cost < eptr2->cost)
			max_edge_cost = eptr2->cost;

		eptr->is_in_tree = 0;

		//We update the degree of the edge's endpoints
		//tail node
		iptr = &(Eptr->nodes[eptr->i]);
		iptr->degree++;
		//head node
		iptr = &(Eptr->nodes[eptr->j]);
		iptr->degree++;
	}

	int inc = 0;
	// Copy of the matching
	for(e=Tptr->n-1;e<m;e++){
		eptr = &(Eptr->edges[e]);
		eptr2 = M[inc];
		eptr->id =e;
		eptr->i = eptr2->i;
		eptr->j = eptr2->j;
		eptr->cost = eptr2->cost;

		if(max_edge_cost < eptr2->cost)
			max_edge_cost = eptr2->cost;

		eptr->is_in_tree = 0;

		//We update the degree of the edge's endpoints
		//tail node
		iptr = &(Eptr->nodes[eptr->i]);
		iptr->degree++;
		//head node
		iptr = &(Eptr->nodes[eptr->j]);
		iptr->degree++;

		inc++;
	}
	Eptr->max_edge_cost = max_edge_cost;

	// Creating the list of edges for every node
	// Calculation of the degree
	int degree; int i; int j;
	for( i = 0 ; i < n ; i++)
	{
		iptr = &(Eptr->nodes[i]);
		degree = iptr->degree;
		iptr->edges = (edge**)malloc(sizeof(edge*)*degree);
	}

	// Adding the edges to the nodes' lists
	int tmp;
	for( e = 0 ;  e < Eptr->m ; e++)
	{
		eptr = &(Eptr->edges[e]);
		i = eptr->i;
		j = eptr->j;

		//tail node
		iptr = &(Eptr->nodes[i]);
		tmp = iptr->tmp_degree;
		iptr->edges[tmp] = eptr;
		iptr->tmp_degree++;

		//head node
		iptr = &(Eptr->nodes[j]);
		tmp = iptr->tmp_degree;
		iptr->edges[tmp] = eptr;
		iptr->tmp_degree++;

	}
}

int inList(int x, int list[],int n){
	int res=0;
	for(int i=0;i<n;i++){
		if(x==list[i]){
			res = 1;
		}
	}
	return res;
}


int listNodes(edge** edges,int m, int nodes[]){
	// Return the number of nodes present in edges and a list of them
	int e,i=0;
	edge* eptr;
	// For each edge we had the node not already present in the list
	for(e=0;e<m;e++){
		eptr=edges[e];
		// printf("(%d,%d) ",eptr->i,eptr->j);
		if(!inList(eptr->i,nodes,i)){
			nodes[i]=eptr->i;
			i++;
		}
		if(!inList(eptr->j,nodes,i)){
			nodes[i]=eptr->j;
			i++;
		}
	}
	// printf("\n");
	// for(int k=0;k<i;k++){
	// 	printf("%d ",nodes[k]);
	// }
	// printf("\n");
	return i;
}


void createGraph(edge** edges,int m,graph* Gptr){
	// printf("Dans createGraph\n");
	Gptr->m = m;
	int nodes[2*m]; // 2*m the maximal number of nodes possible
	int n =listNodes(edges,m,nodes);
	// printf("createGraph recherche indice dans liste de Gptr\n");
	// display_edges(edges,m);
	// display_cycle(nodes,n);
	// printf("\n");
	// printf("Nodes : %d\n",n);
	Gptr->n = n;
	// printf("Fin ini des nombres\n");

	int i,e;
	//We malloc and fill the nodes in their structure
	Gptr->nodes = (node*) malloc(sizeof(node) * n);
	node* iptr;
	for( i = 0 ;  i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		// printf("%d nodes[%d] // ",nodes[i],i);
		iptr->id = nodes[i];
		//The node is initially not in the tree
		iptr->is_in_tree = 0;
		//The node has initially degree zero
		iptr->degree=0;
		iptr->tmp_degree= 0;
	}
	// printf("\n");
	// printf("Noeuds créés\n");
	//We malloc, read and fill the edges in their structure
	Gptr->edges = (edge*) malloc(sizeof(edge)*m);
	double cost;
	edge* eptr; int j;
	double max_edge_cost = -LARGE_NUMBER;
	for( e = 0 ;  e < m ; e++)
	{
		// Retrieve the values of the edge in egdes
		i = edges[e]->i;
		j = edges[e]->j;
		cost = edges[e]->cost;

		// Copy of the values in Gptr
		eptr = &(Gptr->edges[e]);
		eptr->id = e;
		eptr->i = i;
		eptr->j = j;
		eptr->cost = cost;

		//We update the largest edge cost
		if(max_edge_cost < cost)
			max_edge_cost = cost;
		//The edge is initially not in the tree
		eptr->is_in_tree = 0;

		//We update the degree of the edge's endpoints
		//tail node
		// indsea = search_id_node(i,Gptr,Gptr->n);
		// printf("Search ind %d -> index : %d\n",i,indsea);
		
		iptr = &(Gptr->nodes[search_id_node(i,Gptr,Gptr->n)]);
		iptr->degree++;
		//head node
		iptr = &(Gptr->nodes[search_id_node(j,Gptr,Gptr->n)]);
		iptr->degree++;

		// fprintf(stderr,"%d,%d,%d,%lf\n",eptr->id,eptr->i,eptr->j,eptr->cost);
	}
	//We store the largest edge cost in the graph structure
	Gptr->max_edge_cost = max_edge_cost;
	// printf("Création des arêtes done\n");

	//We fill in the neighbors for each node:
	//Note that we are working with undirected graphs
	//We malloc their arrays
	int degree;
	for( i = 0 ; i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		degree = iptr->degree;
		iptr->edges = (edge**) malloc(sizeof(edge*)*degree);
	}
	// printf("Ini list edges node done\n");
	//We loop over the edges to fill the neighbors arrays
	int tmp;
	for( e = 0 ;  e < m ; e++)
	{
		eptr = &(Gptr->edges[e]);
		i = eptr->i;
		j = eptr->j;
		// printf("On rajoute (%d,%d)\n",i,j);
		// printf("Find edge\n");

		//tail node
		iptr = &(Gptr->nodes[search_id_node(i,Gptr,Gptr->n)]);
		// printf("i : %d ",Gptr->nodes[tmp_id].id);
		tmp = iptr->tmp_degree;
		iptr->edges[tmp] = eptr;
		iptr->tmp_degree++;
		// printf("Tail node\n");

		//head node
		iptr = &(Gptr->nodes[search_id_node(j,Gptr,Gptr->n)]);
		// printf("j : %d \n",Gptr->nodes[tmp_id].id);
		// printf("j=%d\n",tmp_id);
		tmp = iptr->tmp_degree;
		iptr->edges[tmp] = eptr;
		iptr->tmp_degree++;
		// printf("Head node\n");
	}
	// printf("Création de Gptr\n");
}

void copyEdges(edge* destE,edge* sourceE){
	destE->id=sourceE->id;
	destE->i=sourceE->i;
	destE->j=sourceE->j;
	destE->cost=sourceE->cost;
	// printf("(%d,%d)\n",destE->i,destE->j);
}

int isABridge(edge** edges, int m, int n, int indE,int nextNode){
	// Determine if e is a bridge in the list of edges names edges
	// edges : table of edges still in the graph
	// m : number of edges
	// e : edge we want to check
	
	// int val =0;
	// int i,j; double cost;
	// int k;
	// int count=0;
	// edge * eptr;
	// // Count the occurences of the edge e
	// for (k=0;k<m;k++){
	// 	eptr = edges[k];
	// 	i=eptr->i;
	// 	j=eptr->j;
	// 	cost=eptr->cost;
	// 	if(e->i==i && e->j==j && e->cost==cost){
	// 		count ++;
	// 	}
	// }
	// // If there is more than one then it is not a bridge
	// if(count==1){
	// 	val = 1;
	// }
	// return val;
	int val = 0;
	// printf("\nisABridge check\nOn enlève (%d,%d)\n",edges[indE]->i,edges[indE]->j);
	// "Suppresion" from the list of edges remaining
	edge* eptr=(edge*)malloc(sizeof(edge*));
	if (edges[m-1]->id!=-1){
		copyEdges(eptr,edges[indE]);
		// printf("Valeur de eptr: %d (%d,%d) cost %f\n",eptr->id,eptr->i,eptr->j,eptr->cost);
		// printf("Valeur de e[indE]: %d (%d,%d) cost %f\n",edges[indE]->id,edges[indE]->i,edges[indE]->j,edges[indE]->cost);
		copyEdges(edges[indE],edges[m-1]);
		copyEdges(edges[m-1],eptr);
		free(eptr);
		//printf("Fin\n");
	}
	graph Gptr;


	// printf("Check of createGraph\n");
	// display_edges(edges,m-1);
	createGraph(edges,m-1,&Gptr);
	//printf("Graph inter edges\n");
	// display_graph(&Gptr);
	tree * Tptr = (tree *) malloc(sizeof(tree));
	initialize_tree(Tptr, n);

	int nodes[2*m]; // 2*m the maximal number of nodes possible
	int nb_noeud_start = listNodes(edges,m,nodes);

	// printf("DFS\n");
	int visited = DFS(&Gptr,nextNode,Tptr);
	// display_tree(Tptr);
	// printf("End DFS\n");
	// printf("visited = %d et Gptr->n = %d\n",visited,Gptr.n);
	// display_tree(Tptr);

	// int visited; node iptr;
	// for(int i=0;i<Tptr->n;i++){
	// 	iptr = Gptr.nodes[i];
	// 	visited = visited+ iptr.visited;
	// 	printf("%d\n",iptr.visited);
	// }
	// printf("visited = %d\n",visited);
	// printf("Nombre de noeuds au départ %d et à l'arrivée %d\n",nb_noeud_start,Gptr.n);
	if(nb_noeud_start!=Gptr.n || visited!=Gptr.n){val=1;}
	
	// if (val){
	// 	printf("Bridge\n");
	// } else {
	// 	printf("Not a bridge\n");
	// }

	if (edges[m-1]->id!=-1){
		copyEdges(eptr,edges[indE]);
		// printf("Valeur de eptr: %d (%d,%d) cost %f\n",eptr->id,eptr->i,eptr->j,eptr->cost);
		// printf("Valeur de e[indE]: %d (%d,%d) cost %f\n",edges[indE]->id,edges[indE]->i,edges[indE]->j,edges[indE]->cost);
		copyEdges(edges[indE],edges[m-1]);
		copyEdges(edges[m-1],eptr);
		//printf("Fin\n");
	}
	free_tree(Tptr);

	// if(prim(&Gptr,nextNode,Tptr)){
	// 	val = 1;
	// 	printf("Bridge!!\n");
	// } else{
	// 	printf("Not a bridge\n");
	// }
	// printf("Fin Bridge\n");
	return val;
}



void supprEdges(edge** e,int m, int k){
	edge* eptr=(edge*)malloc(sizeof(edge*));
	if (e[m-1]->id!=-1){
		e[k]->id=-1;
		copyEdges(eptr,e[k]);
		// printf("Copie 1\n");
		// printf("Valeur de eptr: %d (%d,%d) cost %f\n",eptr->id,eptr->i,eptr->j,eptr->cost);
		// printf("Valeur de e[k]: %d (%d,%d) cost %f\n",e[k]->id,e[k]->i,e[k]->j,e[k]->cost);
		copyEdges(e[k],e[m-1]);
		// printf("Copie 2\n");
		copyEdges(e[m-1],eptr);
		// printf("Copie 3\n");
		free(eptr);
	}

	// edge** edges=(edge**)malloc(sizeof(edge*)*(m-1));
	// if (e[k+1]!=NULL){
	// 	for(int i=0; i<k;k++){
	// 	edges[i] = e[i];
	// 	}
	// 	for(int i=k+1;i<m;i++){
	// 		edges[i-1]=e[i];
	// 	}
	// 	free(e);
	// 	e = edges;
	// } else {
	// 	e = NULL;
	// }
}

void eulerianCycle (graph* Gptr, int cycle[]){
	// Starting node
	int node = 0;

	// Cycle
	int indN = 0;
	cycle[indN] = node;
	indN++;

	// List of the edges that remain to be visited
	edge** edges = (edge**)malloc(sizeof(edge*)*Gptr->m);
	int remainM = Gptr->m;
	for (int i=0;i<remainM;i++){
		edges[i]=&(Gptr->edges[i]);
	}

	
	int e =0; int no;
	edge* eptr; int indE=0; int nextNode;
	//display_edges(edges,remainM);
	while(edges[0]->id!=-1){ // id = -1 when the edge has been deleted
		no = 1;
		e = 0;
		// display_edges(edges,remainM);
		while(no){ // While we don't have the best edge
			// Search for an edge that have the current node for an end
			while(e<remainM && !(edges[e]->i == node) && !(edges[e]->j == node)){
				e++;
			}
			// If we are not at the end of the list we take the current edge
			if (e!=remainM){
				indE = e;
			}
			// Save the edge and the endnode
			eptr=edges[indE];
			if (eptr->i == node){
				nextNode = eptr->j;
			} else {
				nextNode = eptr->i;
			}
			// We take the edge right away if it is not a bridge
			// printf("On regarde (%d,%d) avec noeud suivant %d\n",eptr->i,eptr->j,nextNode);
			if (!isABridge(edges,remainM,Gptr->n,indE,nextNode)){
				no = 0;
				// printf("Pas un bridge\n");
				// printf("Nope c'est pas isABridge\n");
			} else {
				// If the edge is a bridge but also is the only possibility, we take it
				// (it was the last one saved in indE)
				if(e==remainM){
					no = 0;
					// printf("Seul choix\n");
				} else { // We start looking for an other one
					e++;
				}
			}
		}
		// We remove the edge from the list
		// printf("Suppression décidée. On enlève id (%d,%d)\n",Gptr->edges[indE].i,Gptr->edges[indE].j);
		supprEdges(edges,remainM,indE);
		remainM--;
		
		// display_edges(edges,remainM);
		cycle[indN]=nextNode;
		indN++;
		node = nextNode;
	}
	free(edges);
}

void supprRep(int list[], int n, int sortie[]){
	if (n!=0){
		int leng=0;
	for(int i=0;i<n-1;i++){
		if(!inList(list[i],sortie,leng)){
			sortie[leng]=list[i];
			leng++;
		}
	}
	sortie[leng]=list[n-1];
	}
}

void saveFile(graph* Gptr,int sortie[]){

	// Cost calculation
	double cost=0; int id=0;
	for(int k=0; k<Gptr->n;k++){
		id = edge_search(sortie[k],sortie[k+1],Gptr);
		cost=cost+Gptr->edges[id].cost;
	}

	FILE* file=fopen("ouput.csv","w+");

	if(file){
		fprintf(file,"%d,%f",Gptr->n,cost);
		for(int l=0; l<Gptr->n+1;l++){
			fprintf(file,"\n%d",sortie[l]);
		}
		fprintf(file,"\n");
	}
	fclose(file);
}


void Christophides(graph * Gptr){
	int i; int starting_node = 0;

	// Creation tree
	tree * Tptr = (tree *) malloc(sizeof(tree));
	initialize_tree(Tptr, Gptr->n);
	
	printf("Minimum weight spanning tree (Prim)\n\n");
	prim(Gptr, starting_node, Tptr);
	display_tree(Tptr);	

	printf("Checking nodes parity\n");
	parite(Tptr);
	for(i = 0; i < Tptr->n; i++){
		printf("node n°%d : %d /",i, Tptr->parity[i]);
	}
	printf("\n\n");

	//Construction of the complete graph
	printf("Construction of the complete graph\n");
	graph * Cptr = (graph *) malloc(sizeof(graph));
	Complete_construction(Tptr,Gptr,Cptr);
	display_graph(Cptr);
	printf("\n");

	// Matching
	// We have half of the nodes attached in the matching
	edge** M  = (edge **) malloc(sizeof(edge)*(Cptr->n/2)); 
	printf("Matching\n");
	matching(Cptr,M);
	display_edges(M,Cptr->n/2);
	printf("\n");

	// Add spanning tree and matching to create eulerian graph
	graph * Eptr = (graph *) malloc(sizeof(graph));
	eulerian_multigraph(M,Cptr->n/2,Tptr,Eptr);
	printf("Eulerian graph constructed\n");
	display_graph(Eptr);
	printf("\n");

	// // Search of an eulerian walk
	int cycle[Eptr->m+3];
	eulerianCycle(Eptr,cycle);
	printf("Cycle eulérien\n");
	display_cycle(cycle,Eptr->m+1);
	printf("\n");

	// Hamiltonian Cycle
	int sortie[Eptr->n+1];
	supprRep(cycle,Eptr->m+1,sortie);
	printf("Cycle hamiltonien\n");
	display_cycle(sortie,Eptr->n+1);

	saveFile(Gptr,sortie);

	free_tree(Tptr);
	// free_graph(Cptr);
	// free_graph(Eptr);
	// free(M);

}

