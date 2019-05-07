#include <iostream>
#include <string>
#include <vector>
#include <stack> 
#include <algorithm> // find
#include <tuple>
#include <iterator> // distance
#include <bits/stdc++.h> 
#include <fstream>
#include <set>

using namespace std;
string problem_descriptor;
string delimiter = " "; // To parse input

#include "read_graph.cpp"
#include "kernelization.cpp"
#include "graph_function.cpp"
#include "vc_solver.cpp"


//int nodes, edges;
//int actual_edges = 0;

int main(){

  int ** graph;// Graph is store as pointer to list of n-pointers, each of which stores an adj-list of the particular vertex
  int nodes; //Number of vertices
  int edges; //Number of edges;
  int * degree; //Degree of each vertex. Useful while traversing graph
  getInput(&graph, &nodes, &edges, &degree);

  //cout << nodes << "\t"<< edges << endl;
  vector <int> subgraph;
  //vector <int> original_subgraph;
  
  for (int i = 0; i < nodes; i++){
    subgraph.push_back(i);
    //original_subgraph.push_back(i);
  }

  //Random shuffle
  srand(time(0));
  random_shuffle(subgraph.begin(), subgraph.end()); 

  
  /*
  //To print adj-list
  for(int i = 0; i < nodes; i++){
    cout << i << " : ";
    for(int j = 0; j < degree[i]; j++){
      cout << graph[i][j] << " ";
    }
    cout << endl;
    }*/

  vector <int> solution;
  rr_remove_self_loop(graph, nodes, degree, subgraph, solution);
  //cout << "After self loop " << solution.size() << endl;
  //for(vector <int>:: iterator it = solution.begin(); it != solution.end(); it++){
  //  cout << int(*it) << " ";
  //}
  //cout << "\n";  

  rr_remove_clique_nbrs(graph, nodes, degree, subgraph, solution);
  //cout << "After remove clique nbrs (solution) " << solution.size() << endl;
  //cout << "After remove clique nbrs (subgraph) " << subgraph.size() << endl;
    
  rr_half_integrality(graph, nodes, degree, subgraph, solution);
  //cout << "After half integrality (solution) " << solution.size() << endl;
  //cout << "After half integrality (subgraph) " << subgraph.size() << endl;

  rr_remove_dominator(graph, nodes, degree, subgraph, solution);
  //cout << "After remove_dominator (solution) " << solution.size() << endl;
  //cout << "After remove_dominator (subgraph) " << subgraph.size() << endl;

  
  while(is_remove_deg01_applicable(graph, nodes, degree, subgraph)){
  
    //while(rr_remove_deg_one(graph, nodes, degree, subgraph, solution));
    rr_remove_deg_one_v2(graph, nodes, degree, subgraph, solution);
    //cout << "After deg one (Solution) " << solution.size() << endl;
    //cout << "After deg one (Subgraph) " << subgraph.size() << endl;
  
  
    /*for(vector <int>:: iterator it = solution.begin(); it != solution.end(); it++){
      cout << int(*it) + 1 << " ";
      }
      cout << "\n"; */ 
    rr_remove_clique_nbrs(graph, nodes, degree, subgraph, solution);
    //cout << "After remove clique nbrs (solution) " << solution.size() << endl;
    //cout << "After remove clique nbrs (subgraph) " << subgraph.size() << endl;

    rr_remove_deg_zero(graph, nodes, degree, subgraph, solution);
    //cout << "After deg zero (solution) " << solution.size() << endl;
    //cout << "After deg zero (subgraph) " << subgraph.size() << endl;

  
    rr_half_integrality(graph, nodes, degree, subgraph, solution);
    //cout << "After half integrality (solution) " << solution.size() << endl;
    //cout << "After half integrality (subgraph) " << subgraph.size() << endl;
    rr_remove_dominator(graph, nodes, degree, subgraph, solution);
    //cout << "After remove_dominator (solution) " << solution.size() << endl;

  }

  vector<int> reduced_subgraph, increased_solution;
  //tie(reduced_subgraph, increased_solution) = vertex_cover_exact(graph, nodes, degree, subgraph, solution);
  vertex_cover_exact_v1(graph, nodes, degree, subgraph, solution);
  //cout << "solution " << solution.size() << endl;
  //cout << "subgraph " << subgraph.size() << endl;

  /*    vector< vector <int> > conn_componets;
  conn_componets = connected_components(graph, nodes, degree, subgraph);
  cout << conn_componets.size() << endl;

  for(vector< vector <int>>::iterator it = conn_componets.begin(); it != conn_componets.end(); it++){
    for(vector <int>::iterator it1 = (*it).begin(); it1 != (*it).end(); it1++){
      cout << *it1 << " ";
    }
    cout << endl;
    }*/
  
  
  //cut_vertices(graph, nodes, degree, subgraph);
  //cout << subgraph.size() << endl;
  //sort(solution.begin(), solution.end());
  //**** Writing output file starts ****//
  /*ofstream write_output_file;
  write_output_file.open("solution.vc");
  write_output_file << "s vc ";
  write_output_file << nodes << " ";
  write_output_file << solution.size();
  write_output_file << "\n";
  for(vector <int>:: iterator it = solution.begin(); it != solution.end(); it++){
    write_output_file << int(*it) + 1; // +1 to get vertices to 1 to n from 0 to n - 1.
    write_output_file << "\n";  
  }
  write_output_file.close();*/
  //**** Writing output file ends ****//
  // output solution starts //
  cout << "s vc ";
  cout << nodes << " ";
  cout << solution.size();
  cout << "\n";
  for(vector <int>:: iterator it = solution.begin(); it != solution.end(); it++){
    cout << int(*it) + 1; // +1 to get vertices to 1 to n from 0 to n - 1.
    cout << "\n";  
  }

  // output solution ends //
  
  
  for(int i=0; i < nodes; i++){
    free(graph[i]);
  }
  free(graph);
  free(degree);
  return 0;
}
