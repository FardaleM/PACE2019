vector < vector <int> > connected_components(int ** graph, int nodes, int * degree, vector <int> subgraph);
vector <int> cut_vertices(int ** graph, int nodes, int * degree, vector <int> & subgraph);

void print_graph(int ** graph, int nodes, int * degree, vector <int> subgraph){
  // Print given subgraph
  // Input : graph, nodes, degree, and subgraph
  // output: nothing

  bool * alive = (bool *) malloc(sizeof(bool) * nodes);
  if (alive == NULL){
    cout << "No space for alive" << endl;
    exit(EXIT_FAILURE);
  }
  for(int i=0; i < nodes; i++){
    alive[i] = false;
  }

  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    alive[*it] = true;
  }
    
  for(vector <int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    cout << int(*it) << " : ";
    for(int i = 0; i < degree[*it]; i++){
      if (alive[graph[*it][i]] == true){
	cout << graph[*it][i] << " ";
      }
    }
    cout << endl;
  }
  free(alive);
}

void print_degree(int ** graph, int nodes, int * degree, vector <int> subgraph){
  // Print degree in given subgraph
  // Input : graph, nodes, degree, and subgraph
  // output: nothing

  bool * alive = (bool *) malloc(sizeof(bool) * nodes);
  int * deg_subgraph = (int *) malloc(sizeof(int) * nodes);
  if (alive == NULL){
    cout << "No space for alive" << endl;
    exit(EXIT_FAILURE);
  }
  for(int i=0; i < nodes; i++){
    alive[i] = false;
    deg_subgraph[i] = 0;
  }

  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    alive[*it] = true;
  }
    
  for(vector <int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    //cout << int(*it) << " : ";
    for(int i = 0; i < degree[*it]; i++){
      if (alive[graph[*it][i]] == true){
	deg_subgraph[*it] += 1;
      }
    }
    cout << deg_subgraph[*it] << ", ";
  }
  cout << endl;
  free(alive);
  free(deg_subgraph);
}


vector < vector <int> > connected_components(int ** graph, int nodes, int * degree, vector <int> subgraph){
  //Returns connected components of given graph
  //Input : Graph G
  //Ouput : Vector of vectors. Each vector contains a connected component of G. Same as number_of_connected_comp function

  vector< vector <int> > conn_componets;
  
  bool * is_explored = (bool *) malloc(sizeof(bool) * nodes);
  bool * alive = (bool *) malloc(sizeof(bool) * nodes);

  if ((is_explored == NULL) or (alive == NULL)){
    cout << "No space for is_explored " << endl;
    exit(EXIT_FAILURE);
    //return conn_componets; // effectively return 0 with error
  }

  for(int i=0; i < nodes; i++){
    is_explored[i] = false;
    alive[i] = false;
  }

  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    alive[*it] = true;
  }
  
  int xvertex;
  int unexplored_nbr;
  int yvertex;
  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    yvertex = *it; // converting subgraph into vertex
    if (is_explored[yvertex] == true) continue;

    vector <int> component; // Collect vertices in new component
    vector <int> dfs_stack;
    dfs_stack.push_back(yvertex);
    is_explored[yvertex] = true;
    component.push_back(yvertex); // Vertex is pushed in component right after it is first time marked as explored.
    while(!dfs_stack.empty()){
      xvertex = dfs_stack.back();
      unexplored_nbr = -1; // If there is unexplored nbr, this value will be change to something meaningful
      for(int j=0; j < degree[xvertex]; j++){
	if ((alive[graph[xvertex][j]] == true) and (is_explored[graph[xvertex][j]] == false)){
	  unexplored_nbr = graph[xvertex][j];
	  break;
	}
      }
      
      if(unexplored_nbr != -1){
	is_explored[unexplored_nbr] = true;
	component.push_back(unexplored_nbr); // add this vertex to component
	dfs_stack.push_back(unexplored_nbr);
      }else{
	dfs_stack.pop_back(); // remove xvertex from stack
      }
    }
    conn_componets.push_back(component);
  }

  free(is_explored);
  free(alive);
  
  return conn_componets;
}

vector <int> cut_vertices(int ** graph, int nodes, int * degree, vector <int> & subgraph){
  // Returns cut vertices in given graph
  // Input : graph, number of nodes, subgraph
  // Output : a vector containing cut vertices of given graph
  vector <int> articulation_points;

  bool * alive = (bool *) malloc(sizeof(bool *) * nodes);
  bool * is_explored = (bool *) malloc(sizeof(bool) * nodes);
  int * parent = (int *) malloc(sizeof(int) * nodes);
  int * time_visited = (int *) malloc(sizeof(int) * nodes);
  int * dfs_low = (int *) malloc(sizeof(int) * nodes);

  if ((alive == NULL) || (is_explored == NULL) || (parent == NULL) || (time_visited == NULL)){
    cout << "Not enough space for alive/is_explored_parent/time_visited array" << endl;
    exit(EXIT_FAILURE);
  }

  for(int i=0; i< nodes; i++){
    alive[i] = false;
    is_explored[i] = true;
    time_visited[i] = 0;
    dfs_low[i] = 0;
  }

  int xvertex;
  int time = 0; 
  bool is_unexplored_nbr_present;
  vector <int> dfs_stack;
  int unexplored_vertex;
  
  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    alive[*it] = true;
    is_explored[*it] = false;
  }

  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    
    if (is_explored[*it] == true) continue;
    dfs_stack.push_back(*it);
    is_explored[*it] = true;
    parent[*it] = -1; // to denote the root
    time_visited[*it] = time;
    time += 1;

    while(!dfs_stack.empty()){
      xvertex = dfs_stack.back();
      cout << "xvertex " << xvertex << endl;

      is_unexplored_nbr_present = false;
      unexplored_vertex = -1;// If there is an unexplored vertex, this is get some value

      for (int i = 0; i < degree[*it]; i++){
	if (alive[graph[*it][i]] == false){
	  continue;
	}
		
	if ((alive[graph[*it][i]] == true) && (is_explored[graph[*it][i]] == false)) {
	  is_unexplored_nbr_present = true;
	  unexplored_vertex = graph[*it][i];
	  break;
	}
      }    

      if (is_unexplored_nbr_present){
	cout << "Unexplored nbr " << unexplored_vertex << endl;
	dfs_stack.push_back(unexplored_vertex);
	is_explored[unexplored_vertex] = true;
	parent[unexplored_vertex] = *it;
	time_visited[unexplored_vertex] = time;
	time += 1;
      }
      else{
	dfs_stack.pop_back();
      } 
    }  
  }

  free(alive);
  free(is_explored);
  free(parent);
  free(time_visited);
  return articulation_points;
}
/*

void addEdge(int xvertex, int yvertex);
void print_degree(vector <int> & subgraph);

bool is_vertex_cover(vector <int> subgraph, vector <int> solution);
vector <edge_T> maximal_matching(vector <int> original_subgraph);



void print_degree(vector <int> & subgraph){
  // Print degree of all vertices which are in subgraph
  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    int deg_in_subgraph = 0;
    for(vector<int>::iterator it1 = subgraph.begin(); it1 != subgraph.end(); it1++){
      //cout << *it << " " << *it1 << " " << graph[*it][*it1] << endl;
      if (graph[*it][*it1] == 1){
	deg_in_subgraph++;
      }
    }
    cout << deg_in_subgraph << " ";
  }
  cout << endl;
}


bool is_vertex_cover(vector <int> subgraph, vector <int> solution){
  //Returns True if solution is a vertex cover fo graph. False otherwise.
  //Does a sanity check wherever every vertex in solution is present in graph or not.
  //If sanity check works then just deletes vertices in solution from graph and apply reduction rule which deletes deg0 vertices.
  //If graph is empty then return yes, otherwise it returns no.

  //cout << "\n Subgraph ";
  //print_vector(subgraph);
  //cout << "\n Solution ";
  //print_vector(solution);
  //cout << endl;
  
  bool is_element_present;
  for(vector<int>::iterator it = solution.begin(); it != solution.end(); it++){
    is_element_present = ( find(subgraph.begin(), subgraph.end(), int(*it)) != subgraph.end());//is_element_present is True if element is present
    if(!is_element_present){
      cout << "Vertex in solution which is not in graph " << int(*it) << endl;
      return false;
    }
    subgraph.erase(remove(subgraph.begin(), subgraph.end(), int(*it)), subgraph.end());
  }
  
  rr_remove_deg_zero(subgraph, solution); // Passing solution as placeholder. rr does nothing to second argument.

  cout << "\n Remaining Subgraph after removing solution \n";
  print_vector(subgraph);
  
  return subgraph.empty();
}

vector <edge_T> maximal_matching(vector <int> original_subgraph){
  // Returns a maximal matching in graph
  // Input : graph
  // Output : Vector of pair of integers contained in matching
  vector <int> subgraph = original_subgraph; // Local copy of subgraph
  vector <edge_T> matching;
  vector <int> temp_sol;
  int xvertex, yvertex;
  edge_T temp_edge;
  while(!subgraph.empty()){
    rr_remove_deg_zero(subgraph, temp_sol);
    
    xvertex = int(*subgraph.begin());
    for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
      if (graph[xvertex][int(*it)] == 1){
	temp_edge.x = xvertex;
	temp_edge.y = int(*it);
	matching.push_back(temp_edge);
	subgraph.erase(remove(subgraph.begin(), subgraph.end(), xvertex), subgraph.end());
	subgraph.erase(remove(subgraph.begin(), subgraph.end(), int(*it)), subgraph.end());
	break;
      }
    }   
  }
  return matching;
}
*/

