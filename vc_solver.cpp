void solve_deg_two_v1(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution);
void vertex_cover_exact_v1(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution);
void br_highest_degree_vertex_v1(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution);


void solve_deg_two_v1(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution){
  //Solve instances which has maximum degree two
  //Input: Subgraph, Solution
  //Output: none
  //v1 : does not return anything
  //Apply reduction rule to remove deg-1 and deg-0 vertices.
  //Remaining graph is a cycle. Pick a vertex in solution and recall the function

  while(is_remove_deg01_applicable(graph, nodes, degree, subgraph)){
    rr_remove_deg_zero(graph, nodes, degree, subgraph, solution);
    rr_remove_deg_one_v2(graph, nodes, degree, subgraph, solution);
  }
  
  if (subgraph.empty()) return;
  //All vertices now has degree 2. Pick first vertex and put it in solution. 
  int xvertex;
  
  vector < vector <int> > conn_comp = connected_components(graph, nodes, degree, subgraph);
  
  for (vector<vector <int>>::iterator it = conn_comp.begin(); it != conn_comp.end(); it++){
    xvertex =  (*it)[0]; //Pick a vertex in each subgraph and include it in solution. Apply reduction rules after picking this vertex. 
    solution.push_back(xvertex);
    (*it).erase(find((*it).begin(), (*it).end(), xvertex));
    while(is_remove_deg01_applicable(graph, nodes, degree, *it)){
      rr_remove_deg_zero(graph, nodes, degree, *it, solution);
      rr_remove_deg_one_v2(graph, nodes, degree, *it, solution);
    }
  }
  subgraph.clear();
  return;
}

void vertex_cover_exact_v1(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution){
  //Main function to solve vertex cover
  //Input : Subgraph, solution
  //Output : none
  // v1 : returns none
  // Apply reduction rules. Calls br_highest_degree on each connected components. Joins the solution and returns it.

  int xvertex; //Temp vertex
  //print_degree(subgraph);
  
  if (subgraph.empty()) {
    return;
  }
  
  vector < vector <int> > conn_comp = connected_components(graph, nodes, degree, subgraph);
  
  for (vector<vector <int>>::iterator it = conn_comp.begin(); it != conn_comp.end(); it++){
    while(is_remove_deg01_applicable(graph, nodes, degree, *it)){
      rr_remove_deg_zero(graph, nodes, degree, *it, solution);
      rr_remove_deg_one_v2(graph, nodes, degree, *it, solution);
    }

    br_highest_degree_vertex_v1(graph, nodes, degree, *it, solution);
  }
  subgraph.clear();
  return;
}

void br_highest_degree_vertex_v1(int ** graph, int nodes, int * degree, vector <int> & subgraph, vector <int> & solution){
  // Branching Rule
  // Input: Subgraph, Solution
  // Output: none
  //
  // v1 : returns none. simple version. fewer return values
  // Pick highest degree vertex which is of degree 2 or more and makes two branches. In first, v is added in solution, in second N(v) is added in solution.
  //cout << "******* Call of branching function ******" << endl;
  //cout << subgraph.size() << endl;

  if (subgraph.empty()) {
    //cout << "Leaf node " << solution.size() << endl;
    return;
  }
  
  int branching_vertex = -1; // This will get some sensible value
  int max_deg = 0;
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
    
    if (temp_deg > max_deg) {
      branching_vertex = *it;
      max_deg = temp_deg;
    }
  }

  if (max_deg <= 2) {
    solve_deg_two_v1(graph, nodes, degree, subgraph, solution);
    free(alive);
    return;
  }
  
  vector <int> branching_nbrs;
  //cout << "Branching vertex " << branching_vertex << endl;
  //cout << "Deg of branching " << max_deg << endl;
  
  for(int j = 0; j < degree[branching_vertex]; j++){
    if (alive[ graph[branching_vertex][j] ] == true){
      branching_nbrs.push_back(graph[branching_vertex][j]);
    }
  }

  free(alive);

  vector <int> temp_sub_graph1, temp_solution1, temp_sub_graph2, temp_solution2;

  //cout << "*** **** " << endl;
  //cout << "Solution " << solution.size() << endl;
  //cout << "Subgraph  " << subgraph.size() << endl;

  // Creating first instance
  temp_sub_graph1 = subgraph;
  //temp_solution1 = solution;
  temp_sub_graph1.erase(find(temp_sub_graph1.begin(), temp_sub_graph1.end(), branching_vertex));
  temp_solution1.push_back(branching_vertex);

  // Running first instance
  vertex_cover_exact_v1(graph, nodes, degree, temp_sub_graph1, temp_solution1);

  
  //cout << "Temp Solution-1 " << temp_solution1.size() << endl;
  //cout << "Temp Subgraph-1 " << temp_sub_graph1.size() << endl;
  //cout << "Solution " << solution.size() << endl;
  //cout << "Subgraph  " << subgraph.size() << endl;

  
  // Creating second instance
  temp_sub_graph2 = subgraph;
  //temp_solution2 = solution;
  temp_sub_graph2.erase(find(temp_sub_graph2.begin(), temp_sub_graph2.end(), branching_vertex));
  //print_vector(temp_sub_graph2);
  //print_vector(temp_solution2);
  for(vector<int>::iterator it_brnr = branching_nbrs.begin(); it_brnr != branching_nbrs.end(); it_brnr++){
    xvertex = int(*it_brnr);
    //cout << xvertex << endl;
    temp_sub_graph2.erase(find(temp_sub_graph2.begin(), temp_sub_graph2.end(), xvertex));
    temp_solution2.push_back(xvertex);
  }

  // Running second instance
  vertex_cover_exact_v1(graph, nodes, degree, temp_sub_graph2, temp_solution2);

  
  //cout << "Temp Solution-2 " << temp_solution2.size() << endl;
  //cout << "Temp Subgraph-2 " << temp_sub_graph2.size() << endl;
  //cout << "### ####" << endl;
  
  if (temp_solution1.size() <= temp_solution2.size()){
    subgraph = temp_sub_graph1;
    for(vector <int>::iterator it1 = temp_solution1.begin(); it1 != temp_solution1.end(); it1++){
      solution.push_back(*it1);
    }
    //solution = temp_solution1;
    return;
  }

  subgraph = temp_sub_graph2;
  //solution = temp_solution2;
  for(vector <int>::iterator it2 = temp_solution2.begin(); it2 != temp_solution2.end(); it2++){
    solution.push_back(*it2);
  }
  return;
}
