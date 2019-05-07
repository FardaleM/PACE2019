void rr_remove_self_loop(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution);
void rr_remove_deg_zero(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution);
bool rr_remove_deg_one(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution);
void rr_half_integrality(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution);
float * solve_lp(int ** graph, int nodes, int * degree, vector <int> & subgraph);
vector <int> _aug_path (int n, int ** bipartite_graph, int * degree, int * Matching);
int * maximum_bipartite_matching(int n, int ** bipartite_graph,  int * degree_bip_graph);
tuple<int **, int *> _make_aux_bipartite(int ** graph, int nodes, int * degree, vector <int> subgraph);

bool is_remove_deg01_applicable(int ** graph, int nodes, int * degree, vector <int> & subgraph){
  // Returns true if rr_remove_deg_0 or rr_remove_deg_1 is applicable. Returns false otherwise
  
  if (subgraph.empty()) {
    //cout << "Leaf node " << solution.size() << endl;
    return false;
  }
  
  int min_deg = nodes;
  int temp_deg;
  int xvertex;
  bool * alive = (bool *) malloc(sizeof(bool) * nodes);

  if (alive == NULL){
    cout << "No space for alive in br_highest_degree_vertex in vc_solver.cpp " << endl;
    exit(EXIT_FAILURE);
  }

  for(int i = 0; i < nodes; i++){
    alive[i] = false;
  }

  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    alive[*it] = true;
  }

  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    temp_deg = 0;
    for(int j = 0; j < degree[*it]; j++){
      if (alive[graph[*it][j]] == true){
	temp_deg += 1;
      }
    }
    
    if (temp_deg < min_deg) {
      min_deg = temp_deg;
    }
  }

  free(alive);
  if (min_deg < 2){
    return true;
  }

  return false;
}


void rr_remove_self_loop(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution){
  //Reduction Rule
  //Input : global_graph with number of vertices and degree vector, subgraph, partial solution 
  //Output : Nothing
  //If a vertex contains self-loop then include it in solution

  vector<int> self_loop_vertices;// to collect vertices with self-loop

  int * ptr_find; // pointer to find vertex x in adj_list of x

  for(int i = 0; i < nodes; i++){
    ptr_find = find(graph[i], graph[i] + degree[i], i);
    if (ptr_find != graph[i] + degree[i]){
      //vertex i has i in its adj-list
      self_loop_vertices.push_back(i);
    }
  }
  
  for(vector<int>::iterator it = self_loop_vertices.begin(); it != self_loop_vertices.end(); it++){
    solution.push_back(*it);
    subgraph.erase(find(subgraph.begin(), subgraph.end(), int(*it)));
  }
}

void rr_remove_deg_zero(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution){
  //Reduction Rule
  //Input : global_graph with number of vertices and degree vector, subgraph, partial solution
  //Output : Nothing
  //
  // If vertex of degree 0 is present in subgraph then reduction rule removes it.

  bool * alive = (bool *) malloc(sizeof(bool *) * nodes);

  if (alive == NULL){
    cout << "No Space for alive in function rr_remove_deg_zero in kernelization.cpp " << endl;
    exit(EXIT_FAILURE);
  }
  for(int i=0; i < nodes; i++){
    alive[i] = false;
  }
  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    alive[*it] = true;
  }

  int deg_in_subgraph = 0;
  
  vector<int> zero_deg_vertices;
  for (vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    deg_in_subgraph = 0;

    for (int i=0; i < degree[*it]; i++){
      if (alive[graph[*it][i]] == true){
	deg_in_subgraph += 1;
      }
    }    
    
    if (deg_in_subgraph == 0){
      zero_deg_vertices.push_back(*it);
    }
  }

  if (!zero_deg_vertices.empty()){
    for (vector<int>::iterator it = zero_deg_vertices.begin(); it != zero_deg_vertices.end(); it++){
      subgraph.erase(find(subgraph.begin(), subgraph.end(), int(*it))); 
    }
  }

  free(alive);
}


bool rr_remove_deg_one(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution){
  //Reduction Rule
  //Input : global_graph with number of vertices and degree vector, subgraph, partial solution
  //Output : True if reduction rule is applicable, False otherwise
  //
  //If the vertex has degree 1 in the subgraph then include its neighbor in solution and return True. Return False if no such vertex is found.
  bool * alive = (bool *) malloc(sizeof(bool *) * nodes);

  if (alive == NULL){
    cout << "No space for alive array" << endl;
    exit(EXIT_FAILURE);
  }
  
  for(int i=0; i < nodes; i++){
    alive[i] = false;
  }
    
  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    alive[int(*it)] = true;
  }
  
  int deg_in_subgraph = 0;
  int temp_nbr = 0;
  
  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    deg_in_subgraph = 0;
    for (int i=0; i < degree[*it]; i++){
      if (alive[graph[*it][i]] == true){
	deg_in_subgraph += 1;
	temp_nbr = graph[int(*it)][i];
	//cout << temp_nbr << " "; 
      }
      //cout << endl;
    }    
          
    if (deg_in_subgraph == 1){
      solution.push_back(temp_nbr); //Add its nbr to solution
      subgraph.erase(find(subgraph.begin(), subgraph.end(), int(*it)));
      subgraph.erase(find(subgraph.begin(), subgraph.end(), temp_nbr));
      free(alive);
      return true;
    }
  }

  free(alive);
  return false;
}

void rr_remove_deg_one_v2(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution){
  //Reduction Rule
  // Input : global_graph with number of vertices and degree vector, subgraph, partial solution
  // Output : nothing
  // v2 : Runs faster. Expected to run in O(n)
  // If the vertex has degree 1 in the subgraph then include its neighbor in solution and return True. Return False if no such vertex is found.

  bool * alive = (bool *) malloc(sizeof(bool) * nodes);
  int * degree_in_subgraph = (int *) malloc(sizeof(int) * nodes);
  
  if ((alive == NULL) or (degree_in_subgraph == NULL)){
    cout << "No space for alive/degree_in_subgraph array" << endl;
    exit(EXIT_FAILURE);
  }
  
  for(int i=0; i < nodes; i++){
    alive[i] = false;
    degree_in_subgraph[i] = 0;
  }
    
  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    alive[int(*it)] = true;
  }

  stack<int> s;
  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    //degree_in_subgraph[*it] = 0; // vertex is alive so initializing its deg to 0
    for (int i=0; i < degree[*it]; i++){
      if (alive[graph[*it][i]] == true){
	degree_in_subgraph[*it] += 1;
      }
    }
    if(degree_in_subgraph[*it] == 1)
        s.push(*it);
  }

  int temp_nbr;

  while(!s.empty()){
      int u = s.top();
      s.pop();
      if (degree_in_subgraph[u] != 1 || !alive[u])
          continue;

      for(int j=0; j < degree[u]; j++){
          if (alive[graph[u][j]] == true){
              temp_nbr = graph[u][j];
              break;
          }
      }

      solution.push_back(temp_nbr); //Add its nbr to solution
      subgraph.erase(find(subgraph.begin(), subgraph.end(), u));
      subgraph.erase(find(subgraph.begin(), subgraph.end(), temp_nbr));
      alive[u] = false;
      alive[temp_nbr] = false;
      for(int j=0; j < degree[temp_nbr]; j++){
          if (alive[graph[temp_nbr][j]] == true){
              degree_in_subgraph[graph[temp_nbr][j]] -= 1;
              if(degree_in_subgraph[graph[temp_nbr][j]] == 1)
                  s.push(graph[temp_nbr][j]);
          }
      }

  }
  
  free(alive);
  free(degree_in_subgraph);
}


void rr_remove_clique_nbrs(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution){
  // Reduction rule.
  // If N(v) is a clique than include N(v) in the solution.

  bool * alive = (bool *) malloc(sizeof(bool) * nodes);
  
  if ((alive == NULL)){
    cout << "No space for alive/degree_in_subgraph array" << endl;
    exit(EXIT_FAILURE);
  }
  
  for(int i=0; i < nodes; i++){
    alive[i] = false;
  }
    
  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    alive[int(*it)] = true;
  }

  int * ptr_find; // To check if vertex is present in an array or not
  int temp_nbr;
  bool solution_increased = true; // if a vertex is included then solution increases and we check once again whether degree is one or not.
  bool is_nbr_clique;
  bool is_vertex_adj_to_rest;
  int xvertex, yvertex;
  
  while(solution_increased == true){
    solution_increased = false;  

    for(int i=0; i < nodes; i++){
      if (alive[i] == false) continue;
      
      is_nbr_clique = true;
      
      for(int j=0; j < degree[i]; j++){
	xvertex = graph[i][j];
	if (alive[xvertex] == false) continue;
	
	is_vertex_adj_to_rest = true;
	
	for(int k = j + 1; k < degree[i]; k++){
	  yvertex = graph[i][k];
	  if (alive[yvertex] == false) continue;
	  ptr_find = find(graph[xvertex], graph[xvertex] + degree[xvertex], yvertex);
	  if (ptr_find == graph[xvertex] + degree[xvertex]){
	    // xvertex is not adj with yvertex
	    is_vertex_adj_to_rest = false;
	    break;
	  }  
	}

	if (is_vertex_adj_to_rest == false){
	  is_nbr_clique = false;
	  break;
	}
      }

      if (is_nbr_clique == false) continue;

      solution_increased = true;
      alive[i] = false;
      subgraph.erase(find(subgraph.begin(), subgraph.end(), i));
      
      vector <int> clique_nbrs;
      //cout << "Branching vector " << branching_vertex << endl;
  
      for(int j = 0; j < degree[i]; j++){
	if (alive[graph[i][j]] == true){
	  clique_nbrs.push_back(graph[i][j]);
	}
      }

      for(vector<int>::iterator it_clnr = clique_nbrs.begin(); it_clnr != clique_nbrs.end(); it_clnr++){
	xvertex = int(*it_clnr);
	alive[xvertex] = false;
	subgraph.erase(find(subgraph.begin(), subgraph.end(), xvertex));
	solution.push_back(xvertex);
      }      
    }
  }

  free(alive);
}


void rr_remove_dominator(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution){
  // Reduction rule.
  // If N[u] \subseteq N[v] then include v in solution. 
  // In this case, v is called dominator
  
  bool * alive = (bool *) malloc(sizeof(bool) * nodes);
  
  if ((alive == NULL)){
    cout << "No space for alive/degree_in_subgraph array" << endl;
    exit(EXIT_FAILURE);
  }
  
  for(int i=0; i < nodes; i++){
    alive[i] = false;
  }
    
  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    alive[int(*it)] = true;
  }

  int * ptr_find; // To check if vertex is present in an array or not
  int temp_nbr;
  bool solution_increased = true; // if a vertex is included then solution increases and we check once again whether degree is one or not.
  bool is_xvertex_dominator;
  int xvertex, yvertex;
  
  while(solution_increased == true){
    solution_increased = false;  

    for(int i = 0; i < nodes; i++){
      if (alive[i] == false) continue;
      xvertex = i;
      is_xvertex_dominator = false;
      vector <int> alive_cl_nbr_xvertex {xvertex};
      for(int j = 0; j < degree[xvertex]; j++){
	if (alive[graph[xvertex][j]] == true){
	  alive_cl_nbr_xvertex.push_back(graph[xvertex][j]);
	}
      }

      for(vector <int>::iterator it=alive_cl_nbr_xvertex.begin(); it != alive_cl_nbr_xvertex.end(); it++){
	yvertex = *it;
	if (yvertex == xvertex) continue;
	vector <int> alive_cl_nbr_yvertex {yvertex};
	for(int j=0; j < degree[yvertex]; j++){
	  if (alive[graph[yvertex][j]] == true){
	    alive_cl_nbr_yvertex.push_back(graph[yvertex][j]);
	  }
	}

	if (includes(alive_cl_nbr_xvertex.begin(), alive_cl_nbr_xvertex.end(), alive_cl_nbr_yvertex.begin(), alive_cl_nbr_yvertex.end()) == true){
	  // closed nbrhood of yvertex is contained in contained in closed neighborhood in xvertex
	  is_xvertex_dominator = true;
	  break;
	}
      }

      if (is_xvertex_dominator == true){
	solution_increased = true;
	alive[xvertex] = false;
	solution.push_back(xvertex);
	subgraph.erase(find(subgraph.begin(), subgraph.end(), xvertex));
	break;
      }
    }
  }
  
  free(alive);
}


void rr_half_integrality(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution){
  // Reduction Rule. Computes Half-Integral solution for given graph
  // and includes all vertices which got value 1 into the solution
  //Input : global_graph with number of vertices and degree vector, subgraph, partial solution
  //Output : Nothing. Modify subgraph and solution vectors
  
  float * lp_value = solve_lp(graph, nodes, degree, subgraph);
  
  //Vertices to be deleted subgraph. 
  vector <int> to_be_delete;
  // If lp_value[x] = 1 or 0  vertex x will be deleted.
  // If lp_value[x] = 1 then it is inclued in solution.
  for(int i=0; i < subgraph.size(); i++){
    if (lp_value[i] == 0){
      to_be_delete.push_back(subgraph[i]);
    }
    if (lp_value[i] == 1.0){
      to_be_delete.push_back(subgraph[i]);
      solution.push_back(subgraph[i]);
    }
  }

  for(vector <int>::iterator it = to_be_delete.begin(); it != to_be_delete.end(); it++){
    subgraph.erase(find(subgraph.begin(), subgraph.end(), int(*it)));
  }

  free(lp_value);
}

float * solve_lp(int ** graph, int nodes, int * degree, vector <int> & subgraph){
  // Given a graph, solve half-integral LP corresponding to VC.
  // Input : subgraph
  // Output : array of float values, each correponding to a value in subgraph

  //cout << "In solve lp" << endl;
  
  int ** bipartite_graph;
  int * deg_bipart;
  tie(bipartite_graph, deg_bipart) = _make_aux_bipartite(graph, nodes, degree, subgraph);
  // Size of bipartite graph
  int bn = 2 * subgraph.size();

  //for(int i=0; i < subgraph.size(); i++){
  //  if (deg_bipart[i] > degree[subgraph[i]]){
  //    cout << i << " : " << degree[subgraph[i]] << " " << deg_bipart[i] << endl;
  //  }
  //}
  
  int * Matching = maximum_bipartite_matching(bn, bipartite_graph, deg_bipart);
  bool * setZ = (bool *) malloc(sizeof(bool) * bn);
  // Set Z is used to find vertex cover from Matching.
  // Let L, R be two partitions of bipartite graph
  // For us, L = vertices with index < bn/2; R = rest
  // Z := unsaturated vertices in L (say U) union vertices which are rechable
  // from U via alternating path with resp to U.
  // Vertex cover = (L \setminus Z) union (R \cap Z)
  bool * is_explored = (bool *) malloc(sizeof(bool) * bn);
  int * layer = (int *) malloc(sizeof(int) * bn);
  //layer[i] can be odd (1) or even (0). It is initialized to -1

  if ((is_explored == NULL) or (setZ == NULL) or (layer == NULL)){
    cout << "Not enough space for is_explored or setZ or layer" << endl;
    exit(EXIT_FAILURE);
  }

  for (int i=0; i < bn; i++){
    setZ[i] = false;
    is_explored[i] = false;
    layer[i] = -1;
  }
  
  int xvertex, yvertex;  
  bool is_unexplored_nbr_present;
  vector <int> dfs_stack;
  int unexplored_vertex;
  
  for(int i = 0; i < bn/2; i++){
    //Notice change in range. 
    if (Matching[i] != -1) continue;
    //DFS starts from unsaturated vertex in L
    if (is_explored[i] == true) continue;
    
    dfs_stack.clear();
    dfs_stack.push_back(i);
    is_explored[i] = true;
    layer[i] = 0; //First vertex is at even level
    
    while (!dfs_stack.empty()){
      xvertex = dfs_stack.back();
      setZ[xvertex] = true;// Include this vertex in Z
      
      if (layer[xvertex] == 1){	
	// If there is a vertex at odd level which is saturated then we move
	// via Matching (if we can)
	// If Matching[a] = b then add b to stack (if it not already explored) and continue
	// Note that there is no vertex at odd level which is unsaturated and reachable from
	// unsaturated vertex from even level (i.e. L)
	yvertex = Matching[xvertex];
	if (is_explored[yvertex] == false){
	  //cout << "Via Matching, yvertex : " << yvertex << endl;
	  dfs_stack.push_back(yvertex);
	  is_explored[yvertex] = true;
	  layer[yvertex] = 1 - layer[xvertex];
	}else{
	  dfs_stack.pop_back();
	}
	continue;
      }

      is_unexplored_nbr_present = false;
      unexplored_vertex = -1;// If there is an unexplored vertex, this is get some value
      
      for(int j = 0; j < deg_bipart[xvertex]; j++){
	if (is_explored[bipartite_graph[xvertex][j]] == false){
	  is_unexplored_nbr_present = true;
	  unexplored_vertex = bipartite_graph[xvertex][j];
	  break;
	}
      }

      if (is_unexplored_nbr_present){
	//cout << "Unexplored nbr " << unexplored_vertex << endl;
	dfs_stack.push_back(unexplored_vertex);
	is_explored[unexplored_vertex] = true;
	layer[unexplored_vertex] = 1 - layer[xvertex]; // Chage layer from 1 to 0 or 0 to 1
      }
      else{
	dfs_stack.pop_back();
      }
    }
  }
  
  bool * vc_bipartite = (bool *) malloc(sizeof(bool) * bn);
  float * lp_value = (float *) malloc(sizeof(float) * subgraph.size());

  if((vc_bipartite == NULL) or (lp_value == NULL)){
    cout << "No space for vc_bipartite of lp value in kernelization.cpp " << endl;
    exit(EXIT_FAILURE);
  }
  
  for (int i=0; i < bn/2; i++){
    //Notice change is range.
    //For first half, it is negation of setZ[i] value
    vc_bipartite[i] = !setZ[i];
    vc_bipartite[bn/2 + i] = setZ[bn/2 +i];
  }

  for (int i=0; i < subgraph.size(); i++){
    lp_value[i] = 0;
    if (vc_bipartite[i] == true) lp_value[i] += 0.5;
    if (vc_bipartite[subgraph.size() + i] == true) lp_value[i] += 0.5;
  }

  for(int i=0; i < 2 * subgraph.size(); i++){
    free(bipartite_graph[i]);
  }
  free(bipartite_graph);
  free(deg_bipart);
  free(Matching);
  free(is_explored);
  free(setZ);
  free(layer);
  free(vc_bipartite);
  return lp_value;
}

tuple<int **, int *> _make_aux_bipartite(int ** graph, int nodes, int * degree, vector <int> subgraph){
  // Called in half_integrality
  // For a given subgraph, make its aux bipartite graph
  // For a graph G, its auxilary bipartite graph is H where
  // for H has two parts V1(G) and V2(G). For every edge uv
  // in G, add two edges u1v2 and v1u2 to H.
  // Input : vector <int> subgraph. 
  // Output : adj list (pointer of array consisting of pointers) and degree

  //cout << "In make_aux_bipartite" << endl;
  
  bool * alive = (bool *) malloc(sizeof(bool) * nodes);
  int * index_in_subgraph = (int *) malloc(sizeof(int) * nodes);
  //If index_in_subgraph[x] = 3 then x is 3rd element in vector subgraph.
  // -1 implies vertex is not present in subgraph
  if((alive == NULL) or (index_in_subgraph == NULL)){
    cout << "No space for alive or index_in_subgraph in function _make_aux_bipartite in kernelization.cpp " << endl;
    exit(EXIT_FAILURE);
  }
  //Initilization
  for(int i=0; i < nodes; i++){
    alive[i] = false;
    index_in_subgraph[i] = -1;
  }

  int index = 0;
  for(vector<int>::iterator it = subgraph.begin(); it != subgraph.end(); it++){
    alive[*it] = true;
    index_in_subgraph[*it] = index;
    index += 1;
  }

  int n = subgraph.size();
  int bn = 2 * n; //Number of vertices in bipartite graph
  
  int * deg_bipart = (int *) malloc(sizeof(int) * bn);
  if(deg_bipart == NULL){
    cout << "No space for deg_bipart in function _make_aux_bipartite in kernelization.cpp " << endl;
    exit(EXIT_FAILURE);
  }
  for(int i=0; i < bn; i++){
    deg_bipart[i] = 0;
  }
  
  int x, y;
  int deg_in_subgraph;
  
  for(int i=0; i < n; i++){
    x = subgraph[i];
    deg_in_subgraph = 0;
    for (int j=0; j < degree[x]; j++){
      if (alive[graph[x][j]]  == true){
	deg_in_subgraph += 1;
      }
    }

    deg_bipart[i] = deg_in_subgraph;
    deg_bipart[i + n] = deg_in_subgraph;
  }

  //cout << "Vertices in Aux Bipartitie " << bn << endl;
  int ** bipartite_graph;
  bipartite_graph = (int **) malloc(sizeof(int *) * bn);
  if (bipartite_graph == NULL){
    cout << "No space for bipartite_graph in function _make_aux_bipartite in file kernelization.cpp " << endl;
    exit(EXIT_FAILURE);
  }
  for (int i=0; i < bn; i++){
    bipartite_graph[i] = (int *) malloc(sizeof(int) * deg_bipart[i]);
    if (bipartite_graph[i] == NULL){
      cout << "No space for bipartite_graph[i] in function _make_aux_bipartite in file kernelization.cpp " << endl;
      exit(EXIT_FAILURE);
    }
  }

  int index_in_adj_list;
  for (int i=0; i < n; i++){
    x = subgraph[i];
    index_in_adj_list = 0;
    for (int j=0; j < degree[x]; j++){
      y = graph[x][j];
      if (alive[y] == true){
	//We need to store position of y in subgraph in adj_list
	bipartite_graph[i][index_in_adj_list] = index_in_subgraph[y] + n;
	bipartite_graph[i + n][index_in_adj_list] = index_in_subgraph[y];
	index_in_adj_list += 1;
      }
    }
  }

  //changing first half of vertex set. if edge (i, j) is present then in
  //first half, we need to make edge (i, n +j)
  //for(int i = 0; i < n; i++){
  //  for(int j = 0; j < deg_bipart[i]; j++){
  //    bipartite_graph[i][j] += n;
  //  }
  //}
  
  free(alive);
  free(index_in_subgraph);
  return make_pair(bipartite_graph, deg_bipart);
}

vector <int> _aug_path (int n, int ** bipartite_graph, int * deg_bipart, int * Matching){
  // Called in maximum_bipartite_batching
  // Returns aug_path. (It it doesn't exists then it returns empty vector)
  // Input : Size of graph, address of its adj matrix, matching
  // Ouput : a vector containing vertices in aug path

  //cout << "Aug path" << endl;
  
  bool * is_explored = (bool *) malloc(sizeof(bool) * n);
  bool * is_saturated = (bool *) malloc(sizeof(bool) * n);
  int * layer = (int *) malloc(sizeof(int) * n);
  //layer[i] can be odd (1) or even (0). It is initialized to -1
  
  if ((is_explored == NULL) or (is_saturated == NULL) or (layer == NULL)){
    cout << "Not enough space for is_explored or is_saturated or layer" << endl;
    exit(EXIT_FAILURE);
  }
  //
  for (int i=0; i < n; i++){
    is_explored[i] = false;
    layer[i] = -1;
    is_saturated[i] = !(Matching[i] == -1);
  }

  int xvertex, yvertex;  
  bool is_unexplored_nbr_present;
  vector <int> dfs_stack;
  int unexplored_vertex;

  for(int i = 0; i < n; i++){
    if (is_explored[i] == true) continue;

    if (is_saturated[i] == true) continue;
    // DFS should start from unsaturated vertex. If there is a connected component in
    // which all vertices are saturated then we don't need to do DFS on that to find
    // aug path

    dfs_stack.clear();
    dfs_stack.push_back(i);
    is_explored[i] = true;
    layer[i] = 0; //First vertex is at even level
    
    while (!dfs_stack.empty()){
      xvertex = dfs_stack.back();
      //cout << "xvertex : " << xvertex << endl;

      if ((layer[xvertex] == 1) and (is_saturated[xvertex] == false)){
	// If there is vertex at odd level which is not saturated, we have found
	// aug_path
	  free(is_explored);
	  free(is_saturated);
	  free(layer);
	return dfs_stack;
      }

      if ((layer[xvertex] == 1) and (is_saturated[xvertex] == true)){
	// If there is a vertex at odd level which is saturated then we move
	// via Matching (if we can)
	// If Matching[a] = b then add b to stack (if it not already explored) and continue
	yvertex = Matching[xvertex];
	if (is_explored[yvertex] == false){
	  //cout << "Via Matching, yvertex : " << yvertex << endl;
	  dfs_stack.push_back(yvertex);
	  is_explored[yvertex] = true;
	  layer[yvertex] = 1 - layer[xvertex];
	}else{
	  dfs_stack.pop_back();
	}
	continue;
      }

      // Condition layer[xvertex] == 0 and is_satureated[xvertex] == false
      // do not occure as DFS should come at this stage via layer[xvertex] = odd
      // first two conditions make sure that we don't need reach here.

      // For remaining part, layer[xvertex] = 0 and is_satureated[xvertex] = true
      is_unexplored_nbr_present = false;
      unexplored_vertex = -1;// If there is an unexplored vertex, this is get some value

      //cout << "Deg of xvertex " << deg_bipart[xvertex] << endl;
      for(int j = 0; j < deg_bipart[xvertex]; j++){
	//cout << "Nbr of xvertex " << bipartite_graph[xvertex][j] << endl;
	if (is_explored[bipartite_graph[xvertex][j]] == false){
	  is_unexplored_nbr_present = true;
	  unexplored_vertex = bipartite_graph[xvertex][j];
	  break;
	}
      }

      if (is_unexplored_nbr_present){
	//cout << "Unexplored nbr " << unexplored_vertex << endl;
	dfs_stack.push_back(unexplored_vertex);
	is_explored[unexplored_vertex] = true;
	layer[unexplored_vertex] = 1 - layer[xvertex]; // Chage layer from 1 to 0 or 0 to 1
      }
      else{
	dfs_stack.pop_back();
      }
    }
  }
  
  free(is_explored);
  free(is_saturated);
  free(layer);
  return dfs_stack;
}


int * maximum_bipartite_matching(int n, int ** bipartite_graph, int * deg_bipart){
  // Returns maximum matching in bipartite graph
  // Input: int (size of graph); address of adj matrix
  // Ouput : pointer to an array of size n. M[a] = b then a is mapped to b
  // If M[a] = -1 then a is unmatched.

  int * Matching;
  Matching = (int *) malloc(sizeof(int) * n);
  if(Matching == NULL){
    cout << "No space for Matching in function maximum_bipartite_matching in file kernelization.cpp " << endl;
    exit(EXIT_FAILURE);
  }
  for(int i=0; i < n; i++){
    Matching[i] = -1;
  }

  vector <int> path;
  path = _aug_path(n, bipartite_graph, deg_bipart, Matching);
  
  while(!path.empty()){
    for(int i=0; i < path.size(); i+=2){
      //Notice range and increment
      Matching[path[i]] = path[i + 1];
      Matching[path[i + 1]] = path[i];
    }
    //Find an aug_path
    path = _aug_path(n, bipartite_graph, deg_bipart, Matching);
  }

  // Check if matching is correct
  /*
  int * appred_mat;
  appred_mat = (int *) malloc(sizeof(int) * n);
  for(int i=0; i < n; i++){
    appred_mat[i] = 0;
  }
  for (int i=0; i < n; i++){
    if (Matching[i] != -1) appred_mat[Matching[i]] += 1;
  }
  for(int i=0; i < n; i++){
    if (appred_mat[i] > 1)
      cout << "Vertex appeared multiple time ";
  }
  cout << endl; */
  return Matching;
}
