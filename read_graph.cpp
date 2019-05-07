bool getInput(int *** ptr_graph, int * ptr_nodes, int * ptr_edges, int ** ptr_degree){
  struct edge_T{
    int xvertex;
    int yvertex;
  };
  edge_T * edge_list; //To store all edges
  edge_T temp_edge;
  temp_edge.xvertex = -1;
  temp_edge.yvertex = -1;
  int * raw_degree;
  string temp_input_line;
  size_t pos = 0;
  string token;
  int index_in_edge_list = 0;
  int x, y;
  
  while(getline(cin, temp_input_line)){
    if (temp_input_line[0] == 'c') continue;
    pos = 0;
    if (temp_input_line[0] == 'p'){
      //Line starting with p is of the form "p td #vertices #edges"
      pos = temp_input_line.find(delimiter);
      token = temp_input_line.substr(0, pos);
      temp_input_line.erase(0, pos + delimiter.length());
      //deleted p

      pos = temp_input_line.find(delimiter);
      token = temp_input_line.substr(0, pos);
      problem_descriptor = token;
      temp_input_line.erase(0, pos + delimiter.length());
      //deleted problem_descriptor

      pos = temp_input_line.find(delimiter);
      token = temp_input_line.substr(0, pos);
      *ptr_nodes = stoi(token);
      temp_input_line.erase(0, pos + delimiter.length());
      //deleted nodes

      *ptr_edges = stoi(temp_input_line);
      //only edges are remaining
      raw_degree = (int *) malloc(sizeof(int) * (*ptr_nodes));
      //raw_degree is used maininly to store vertices as they come. Used to remove duplicates
      edge_list = (edge_T *) malloc(sizeof(edge_T) * (*ptr_edges)); 

      if ((edge_list == NULL) or (raw_degree == NULL)){
	cout << "Not enough space for raw_degree or edge_list in function read_graph.cpp" << endl;
	exit(EXIT_FAILURE);
	//return false;
      }

      for(int i = 0; i < (*ptr_nodes); i++){
	raw_degree[i] = 0;//.xvertex = -1;// = temp_edge;
      }
      
      for(int i = 0; i < (*ptr_edges); i++){
	edge_list[i] = temp_edge;//.xvertex = -1;// = temp_edge;
      }

      continue;
    }

    pos = temp_input_line.find(delimiter);
    token = temp_input_line.substr(0, pos);
    x = stoi(token);
    temp_input_line.erase(0, pos + delimiter.length());
    y = stoi(temp_input_line);


    //Vertices start from 1 to nodes.
    raw_degree[x - 1] += 1;
    raw_degree[y - 1] += 1;
    edge_list[index_in_edge_list].xvertex = x - 1;
    edge_list[index_in_edge_list].yvertex = y - 1;
    index_in_edge_list += 1;
  }

  //Constructing graph from edge_list
  //To store position in adj_list
  int * adj_list_pos;
  //adj_list_pos = (int *) malloc(sizeof(int *) * (*ptr_nodes));
  *ptr_degree = (int *) malloc(sizeof(int) * (*ptr_nodes));
  *ptr_graph = (int **) malloc(sizeof(int *) * (*ptr_nodes));
  
  if ((*ptr_graph == NULL) or (*ptr_degree == NULL)){
    cout << "Not enough space to store graph or degree in read_graph " << endl;
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < (*ptr_nodes); i++){
    (*ptr_degree)[i] = 0; // setting degree 0
    (*ptr_graph)[i] = (int *) malloc(sizeof(int) * raw_degree[i]);
    if ((*ptr_graph)[i] == NULL){
      cout << "Not enough space to store graph at vertex " << i << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  int * ptr_find; // To check if vertex is present in an array or not

  for (int i = 0; i < *ptr_edges; i++){
    x = edge_list[i].xvertex;
    y = edge_list[i].yvertex;
    // if y is present in adj list of x graph then continue.
    ptr_find = find((*ptr_graph)[x], (*ptr_graph)[x] + (*ptr_degree)[x], y);
    if (ptr_find != (*ptr_graph)[x] + (*ptr_degree)[x]){
      // y is already present in adj list of x. Continue.
      continue;
    }
    // if x is present in adj list of y graph then continue.
    ptr_find = find((*ptr_graph)[y], (*ptr_graph)[y] + (*ptr_degree)[y], x);
    if (ptr_find != (*ptr_graph)[y] + (*ptr_degree)[y]){
      // x is already present in adj list of y. Continue.
      continue;
    }

    //(*ptr_graph)[x][adj_list_pos[x]] = y;
    (*ptr_graph)[x][(*ptr_degree)[x]] = y;
    //adj_list_pos[x] += 1;
    (*ptr_degree)[x] += 1;
    //(*ptr_graph)[y][adj_list_pos[y]] = x;
    (*ptr_graph)[y][(*ptr_degree)[y]] = x;
    //adj_list_pos[y] += 1;
    (*ptr_degree)[y] += 1;
  }

  free(raw_degree);
  free(edge_list);
  return true;
}
