#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <boost/config.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <tuple>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <typeinfo>
using namespace std;


#pragma region 
// define an adjacency list with edge weights
typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty; // define edge weight property
typedef boost::property<boost::vertex_name_t, double> VertexWeightProperty; // define node weight property; note that: vertex_index_t is not mutable
typedef boost::adjacency_list<boost::setS, boost::vecS,
	boost::undirectedS, VertexWeightProperty, EdgeWeightProperty> graph; // define all the graph properties
typedef boost::graph_traits<graph>::adjacency_iterator AdjacencyIterator;
#pragma endregion define graph property  


#pragma region  
graph read_data_with_terminals(string file_name) {

	int V_num; // vertex number
	int P_num; // number of positive vertices
	int E_num; // edge number
	int v1;
	int v2;
	double weight;
	graph input_graph; // define the adjacency list of the input graph; there is no need to define the V_num
	string line_content;
	ifstream myfile(file_name); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		while (getline(myfile, line_content)) // read file line by line
		{
			// parse the sting：line_content
			list<string> Parsed_content;
			std::string delimiter = " "; // the delimiter
			size_t pos = 0;
			std::string token;
			while ((pos = line_content.find(delimiter)) != std::string::npos) {
				// find(const string& str, size_t pos = 0) function returns the position of the first occurrence of str in the string, or npos if the string is not found.
				token = line_content.substr(0, pos);
				// The substr(size_t pos = 0, size_t n = npos) function returns a substring of the object, starting at position pos and of length npos
				Parsed_content.push_back(token); // store the subtr to the list
				line_content.erase(0, pos + delimiter.length()); // remove the front substr and the first delimiter
			}
			Parsed_content.push_back(line_content); // store the subtr to the list
			if (!Parsed_content.front().compare("Nodes")) // when it's equal, compare returns 0
			{
				Parsed_content.pop_front();
				V_num = atoi(Parsed_content.front().c_str()); // convert string to int
				for (int i = 0; i < V_num; i++) {
					boost::add_vertex(i, input_graph);
					boost::put(boost::vertex_name_t(), input_graph, i, 0);
				}
			}
			else if (!Parsed_content.front().compare("Edges"))
			{
				Parsed_content.pop_front();
				E_num = atoi(Parsed_content.front().c_str());
			}
			else if (!Parsed_content.front().compare("E"))
			{
				Parsed_content.pop_front(); // remove E, expose v1
				v1 = atoi(Parsed_content.front().c_str()) - 1;
				Parsed_content.pop_front(); // remove v1, expose v2
				v2 = atoi(Parsed_content.front().c_str()) - 1;
				Parsed_content.pop_front(); // remove v2, expose weight
				weight = stof(Parsed_content.front().c_str());
				boost::add_edge(v1, v2, weight, input_graph); // add edge
			}
			else if (!Parsed_content.front().compare("Terminals"))
			{
				Parsed_content.pop_front();
				P_num = atoi(Parsed_content.front().c_str());
			}
			else if (!Parsed_content.front().compare("TP"))
			{
				Parsed_content.pop_front(); // remove TP, expose v1
				v1 = atoi(Parsed_content.front().c_str()) - 1;
				Parsed_content.pop_front(); // remove v1, expose weight
				boost::put(boost::vertex_name_t(), input_graph, v1, stof(Parsed_content.front().c_str()));
			}
			else if (!Parsed_content.front().compare("CTP"))
			{
				Parsed_content.pop_front(); // remove CTP, expose v1
				v1 = atoi(Parsed_content.front().c_str()) - 1;
				boost::put(boost::vertex_name_t(), input_graph, v1, 1e8);  // give terminals big prize
			}
		}

		// check number of vertices
		std::cout << "|V|= " << num_vertices(input_graph);
		std::cout << "  |P|= " << P_num;
		// check number of edges
		std::cout << "  |E|= " << num_edges(input_graph) << endl;
		// print errors
		if (V_num != num_vertices(input_graph)) {
			std::cout << "Error: the number of the input vertices is not right." << endl;
		}
		if (E_num != num_edges(input_graph)) {
			std::cout << "Error: the number of the input edges is not right." << endl;
		}
		// connectivity
		std::vector<int> component(num_vertices(input_graph)); // vertex i is in component[i]; No.component from 0
		int cpn_num = connected_components(input_graph, &component[0]); // the number of component; decrease component
		if (cpn_num > 1) {
			std::cout << "cpn_num: " << cpn_num << endl;
			std::cout << "Error: cpn_num>1" << endl;
		}
		return input_graph;

		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << file_name << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
}
#pragma endregion read_data_with_terminals 2016年12月8日


#pragma region 

graph read_data_staynerd_dimacs(graph input_graph, string staynerd_dimacs_file) {
	int v1;
	int v2;
	double weight;
	int N = num_vertices(input_graph);
	graph output_graph;
	for (int i = 0; i < N; i++) {
		boost::add_vertex(i, output_graph); // input node
		boost::put(boost::vertex_name_t(), output_graph, i, get(boost::vertex_name_t(), input_graph, i)); // input node weight
	}

	string line_content;
	ifstream myfile(staynerd_dimacs_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		while (getline(myfile, line_content)) // read file line by line
		{
			// parse the sting：line_content
			list<string> Parsed_content;
			std::string delimiter = " "; // the delimiter
			size_t pos = 0;
			std::string token;
			while ((pos = line_content.find(delimiter)) != std::string::npos) {
				// find(const string& str, size_t pos = 0) function returns the position of the first occurrence of str in the string, or npos if the string is not found.
				token = line_content.substr(0, pos);
				// The substr(size_t pos = 0, size_t n = npos) function returns a substring of the object, starting at position pos and of length npos
				Parsed_content.push_back(token); // store the subtr to the list
				//cout << token << std::endl; // print the front substr
				line_content.erase(0, pos + delimiter.length()); // remove the front substr and the first delimiter
			}
			//std::cout << line_content << std::endl; // this is the last substr in this line
			Parsed_content.push_back(line_content); // store the subtr to the list

			if (!Parsed_content.front().compare("E"))
			{
				Parsed_content.pop_front(); // remove E, expose v1
				v1 = atoi(Parsed_content.front().c_str()) - 1;
				Parsed_content.pop_front(); // remove v1, expose v2
				v2 = atoi(Parsed_content.front().c_str()) - 1;
				boost::add_edge(v1, v2, get(boost::edge_weight_t(), input_graph, boost::edge(v1, v2, input_graph).first), output_graph); // add edge
			}
		}

		myfile.close(); //close the file
		return output_graph;
		
	}
	else
	{
		std::cout << "Unable to open file " << staynerd_dimacs_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}

}

#pragma endregion read_data_staynerd_dimacs


#pragma region

double net_cost(graph input_graph) {

	double included_cost = 0;
	double missed_prize = 0;
	int N = num_vertices(input_graph); // number of vertices
	graph::out_edge_iterator eit1, eend1;

	if (num_edges(input_graph) == 0) {
		double sum_node_weight = get(boost::vertex_name_t(), input_graph, 0);
		double biggest_node_weight = get(boost::vertex_name_t(), input_graph, 0); // initialize as node weight 0
		for (int i = 1; i < N; i++) {
			double x = get(boost::vertex_name_t(), input_graph, i);
			sum_node_weight = sum_node_weight + x;
			if (biggest_node_weight < x) {
				biggest_node_weight = x; // find the maximal node weight
			}
		}
		missed_prize = sum_node_weight - biggest_node_weight;
	}
	else {
		std::vector<int> component(N); // vertex i is in component[i]; No.component from 0
		int cpn_num = connected_components(input_graph, &component[0]); // the number of component; decrease component
		if (cpn_num != N - num_edges(input_graph)) {
			cout << "Solution Error! This is not a simple tree!" << endl;
			getchar();
		}
		for (int i = 0; i < N; i++) {
			if (in_degree(i, input_graph) > 0) { // the included vertex
				tie(eit1, eend1) = boost::out_edges(i, input_graph); // adjacent_vertices of i
				for_each(eit1, eend1,
					[&input_graph, &i, &included_cost](graph::edge_descriptor it1)
				{
					int j = boost::target(it1, input_graph); // the adjacent vetex to i
					if (j > i) { // don't overcount an edge
						included_cost = included_cost + get(boost::edge_weight_t(), input_graph, boost::edge(i, j, input_graph).first); // the edge cost
					}
				});
			}
			else {
				missed_prize = missed_prize + get(boost::vertex_name_t(), input_graph, i);
			}
		}
	}

	return included_cost + missed_prize;
}

#pragma endregion net_cost  2017年4月17日


#pragma region
graph GPrA(graph input_graph, double& GPrA_time) {

	int N = num_vertices(input_graph); // number of vertices
	vector<double> nw1(N); // the nw value for finding R, it will be decreased after the check
	vector<double> nw2(N); // the nw value for pruning
	for (int i = 0; i < N; i++) {
		nw1[i] = get(boost::vertex_name_t(), input_graph, i);
		nw2[i] = get(boost::vertex_name_t(), input_graph, i);
	}
	vector<bool> pcheck1(N); // true means it has been checked 
	int num_check1 = 0; // the number of checked vertices to find R
	vector<bool> pcheck2(N); // true means it has been checked 
	int num_check2 = 0; // the number of checked vertices for pruning
	vector<int> pdegree1(N); // the processing degree to find R
	vector<int> pdegree2(N); // the processing degree for pruning
	vector<int> leaf1; // the leaves for finding R
	vector<int> leaf2; // the leaves for pruning

	for (int i = 0; i < N; i++) {
		pdegree1[i] = in_degree(i, input_graph); // decrease pdegree
		pdegree2[i] = pdegree1[i];
		if (pdegree1[i] == 0) {
			pcheck1[i] = true; // check disconnected vertices
			num_check1++;
			pcheck2[i] = true;
			num_check2++;
		}
		else if (pdegree1[i] == 1) {
			leaf1.insert(leaf1.end(), i);
			leaf2.insert(leaf2.end(), i);
		}
	}


	graph::out_edge_iterator eit, eend;
	AdjacencyIterator ai, a_end;
	int leaf_num1 = N - num_check1 - 1; // the number of leaves you need to process
	int leaf_num2 = N - num_check1 - 1; // the number of leaves you need to process

	auto begin_time = std::chrono::high_resolution_clock::now(); // start time

																 //this version is similar to the version below
	int k = 0;
	while (k < leaf_num1) {
		int i = leaf1[k];
		k++;
		tie(eit, eend) = boost::out_edges(i, input_graph); // adjacent_vertices of i
		for_each(eit, eend,
			[&input_graph, &pcheck1, &i, &nw1, &pdegree1, &num_check1, &leaf1](graph::edge_descriptor it)
		{
			int j = boost::target(it, input_graph);
			if (pcheck1[j] == false) {
				double cost = get(boost::edge_weight_t(), input_graph, boost::edge(i, j, input_graph).first);
				if (cost < nw1[i]) {
					nw1[j] = nw1[j] + nw1[i] - cost; // decrease nw1[j]
				}
				pcheck1[i] = true; // i has been checked, there is no need to delete it in this phase
				pdegree1[j]--;// decrease pdegree[j]
				if (pdegree1[j] == 1) {
					leaf1.insert(leaf1.end(), j); // it's fine to insert in the end, but not in the biginning
				}
				// break; // how to break a for_each???
			}
		});
	}

	//// find the Root, which is the mark of the optimal prunning result
	//while (num_check1 < N - 1) { 
	//	// there will be a single vertex left unchecked (but its nw value will be decreased)
	//	//// the version below is slower
	//	//while (leaf1.size() > 0) {
	//	//	int i = leaf1[0];
	//	//	leaf1.erase(leaf1.begin());
	//	//	tie(eit, eend) = boost::out_edges(i, input_graph); // adjacent_vertices of i
	//	//	for_each(eit, eend,
	//	//		[&input_graph, &pcheck1, &i, &nw1, &pdegree1, &num_check1, &leaf1](graph::edge_descriptor it)
	//	//	{
	//	//		int j = boost::target(it, input_graph);
	//	//		if (pcheck1[j] == false) {
	//	//			double cost = get(boost::edge_weight_t(), input_graph, boost::edge(i, j, input_graph).first);
	//	//			if (cost < nw1[i]) {
	//	//				nw1[j] = nw1[j] + nw1[i] - cost; // decrease nw1[j]
	//	//			}
	//	//			pcheck1[i] = true; // i has been checked, there is no need to delete it in this phase
	//	//			num_check1++;
	//	//			pdegree1[j]--;// decrease pdegree[j]
	//	//			if (pdegree1[j] == 1) {
	//	//				leaf1.insert(leaf1.end(), j); // it's fine to insert in the end, but not in the biginning; note for (int k = 0; k < leaf1.size(); k++)
	//	//			}
	//	//			// break; // how to break a for_each???
	//	//		}
	//	//	});
	//	//}
	//	//// the version below is fast
	//	for (int k = 0; k < leaf1.size(); k++)
	//	{
	//		int i = leaf1[k];
	//		if (pdegree1[i] == 1) {
	//			tie(eit, eend) = boost::out_edges(i, input_graph); // adjacent_vertices of i
	//			for_each(eit, eend,
	//				[&input_graph, &pcheck1, &i, &nw1, &pdegree1, &num_check1, &leaf1, &k](graph::edge_descriptor it)
	//			{
	//				int j = boost::target(it, input_graph);
	//				if (pcheck1[j] == false) {
	//					double cost = get(boost::edge_weight_t(), input_graph, boost::edge(i, j, input_graph).first);
	//					if (cost < nw1[i]) {
	//						nw1[j] = nw1[j] + nw1[i] - cost; // decrease nw1[j]
	//					}
	//					pcheck1[i] = true; // i has been checked, there is no need to delete it in this phase
	//					num_check1++;
	//					pdegree1[i] = 0; // it's not the leaf any more
	//					pdegree1[j]--; // decrease pdegree[j]
	//					if (pdegree1[j] == 1) {
	//						leaf1.insert(leaf1.end(), j); // it's fine to insert in the end, but not in the biginning; note for (int k = 0; k < leaf1.size(); k++)
	//					}
	//					// break; // how to break a for_each???
	//				}
	//			});
	//		}
	//		//// the version below is slower than that above
	//		//boost::tie(ai, a_end) = boost::adjacent_vertices(i, input_graph);
	//		//for (; ai != a_end; ai++) {
	//		//	int j = *ai;
	//		//	if (pcheck1[j] == false) {
	//		//		double cost = get(boost::edge_weight_t(), input_graph, boost::edge(i, j, input_graph).first);
	//		//		if (cost < nw1[i]) {
	//		//			nw1[j] = nw1[j] + nw1[i] - cost; // decrease nw1[j]
	//		//		}
	//		//		pcheck1[i] = true; // i has been checked, there is no need to delete it in this phase
	//		//		num_check1++;
	//		//		pdegree1[i] = 0; // it's not the leaf any more
	//		//		pdegree1[j]--;// decrease pdegree[j]
	//		//if (pdegree1[j] == 1) {
	//		//	leaf1.insert(leaf1.end(), j);
	//		//}
	//		//		break; // how to break a for_each???
	//		//	}
	//		//}
	//	}
	//}

	// R is the vertex with the biggest nw
	int R = 0;
	double max = nw1[0];
	for (int i = 1; i < N; i++) {
		if (nw1[i] > max) {
			max = nw1[i];
			R = i;
		}
	}

	// Strong pruning tree
	graph output_graph = input_graph; // the output graph

									  //this version is similar to the version below
	k = 0;
	while (k < leaf_num2 + 1) { // since R is ignored, it must be leaf_num2+1
		int i = leaf2[k];
		k++;
		if (i != R) {
			tie(eit, eend) = boost::out_edges(i, output_graph); // adjacent_vertices of i
			for_each(eit, eend,
				[&output_graph, &pcheck2, &i, &nw2, &pdegree2, &num_check2, &leaf2, &k](graph::edge_descriptor it)
			{
				int j = boost::target(it, output_graph);
				if (pcheck2[j] == false) {
					double cost = get(boost::edge_weight_t(), output_graph, boost::edge(i, j, output_graph).first);
					if (cost < nw2[i]) {
						nw2[j] = nw2[j] + nw2[i] - cost; // decrease nw2[j]
					}
					else {
						boost::remove_edge(i, j, output_graph); // remove edge(i,j)	
					}
					pcheck2[i] = true; // i has been checked
					pdegree2[j]--;// decrease pdegree[j]
					if (pdegree2[j] == 1) {
						leaf2.insert(leaf2.end(), j);
					}
					// break; // how to break a for_each???
				}
			});
		}
	}

	//while (num_check2 < N - 1) {
	//	for (int k = 0; k < leaf2.size(); k++)
	//	{
	//		int i = leaf2[k];
	//		if (pdegree2[i] == 1 && i != R) {
	//			tie(eit, eend) = boost::out_edges(i, output_graph); // adjacent_vertices of i
	//			for_each(eit, eend,
	//				[&output_graph, &pcheck2, &i, &nw2, &pdegree2, &num_check2, &leaf2, &k](graph::edge_descriptor it)
	//			{
	//				int j = boost::target(it, output_graph);
	//				if (pcheck2[j] == false) {
	//					double cost = get(boost::edge_weight_t(), output_graph, boost::edge(i, j, output_graph).first);
	//					if (cost < nw2[i]) {
	//						nw2[j] = nw2[j] + nw2[i] - cost; // decrease nw2[j]
	//					}
	//					else {
	//						boost::remove_edge(i, j, output_graph); // remove edge(i,j)	
	//					}
	//					pcheck2[i] = true; // i has been checked
	//					num_check2++;
	//					pdegree2[i] = 0; // it's not the leaf any more
	//					pdegree2[j]--;// decrease pdegree[j]
	//					if (pdegree2[j] == 1) {
	//						leaf2.insert(leaf2.end(), j);
	//					}
	//					// break; // how to break a for_each???
	//				}
	//			});
	//			// the version below is slower than that above
	//			//boost::tie(ai, a_end) = boost::adjacent_vertices(i, input_graph);
	//			//for (; ai != a_end; ai++) {
	//			//	int j = *ai;
	//			//	if (pcheck2[j] == false) {
	//			//		double cost = get(boost::edge_weight_t(), output_graph, boost::edge(i, j, output_graph).first);
	//			//		if (cost < nw2[i]) {
	//			//			nw2[j] = nw2[j] + nw2[i] - cost; // decrease nw2[j]
	//			//		}
	//			//		else {
	//			//			boost::remove_edge(i, j, output_graph); // remove edge(i,j)
	//			//			//// check
	//			//			//std::cout << "output_graph net_cost: " << net_cost(output_graph) << endl; // this line causes errors; becuase edge_descriptor? why?
	//			//		}
	//			//		pcheck2[i] = true; // i has been checked
	//			//		num_check2++;
	//			//		pdegree2[i] = 0; // it's not the leaf any more
	//			//		pdegree2[j]--;// decrease pdegree[j]
	//			//if (pdegree2[j] == 1) {
	//			//	leaf2.insert(leaf2.end(), j);
	//			//}
	//			//		break; 
	//			//	}
	//			//}
	//		}
	//	}
	//}

	// deleted disconnected parts
	std::vector<int> component(N); // vertex i is in component[i]; No.component from 0
	int cpn_num = connected_components(output_graph, &component[0]); // the number of component; decrease component
	for (int i = 0; i < N; i++) {
		if (component[i] != component[R]) { // disconnected vertex
			clear_vertex(i, output_graph); // clear_vertex removes adjacent vertices, but not node weight
		}
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count(); // Nanosecond
	GPrA_time = runningtime / 1e6;

	return output_graph;

}
#pragma endregion  GPrA   2016年11月26日


#pragma region 

graph find_MST_Yahui(graph input_graph, double& MST_time) {
	auto begin_time = std::chrono::high_resolution_clock::now(); // start time

	int N = num_vertices(input_graph);
	graph output_graph(N);
	for (int i = 0; i < N; i++) {
		boost::put(boost::vertex_name_t(), output_graph, i, get(boost::vertex_name_t(), input_graph, i)); // put inside the node weight
	}
	vector <boost::graph_traits<graph>::vertex_descriptor> p(N); // minimum_spanning_tree traits

	if (in_degree(0, input_graph) > 0) { // 0 is connected
										 // find minimum_spanning_tree
		prim_minimum_spanning_tree(input_graph, &p[0]); // it can only be &p[0], and 0 must be connected in MST;
														// print edges in minimum_spanning_tree
		for (int i = 1; i != p.size(); ++i) { // p[0]=0;
			if (p[i] != i) {
				boost::add_edge(i, p[i], get(boost::edge_weight_t(), input_graph, boost::edge(i, p[i], input_graph).first), output_graph); // add edge
			}
		}
	}
	else { // 0 is disconnected
		int v1;
		for (int i = 1; i < N; i++) {
			if (in_degree(i, input_graph) > 0) {
				v1 = i;
				break;
			}
		}
		boost::add_edge(0, v1, 1, input_graph); // add edge (0,v1)
												// find minimum_spanning_tree
		prim_minimum_spanning_tree(input_graph, &p[0]); // it can only be &p[0]; if 0 is disconnected, you need a fake edge to connect it
														// print edges in minimum_spanning_tree
		for (int i = 1; i != p.size(); ++i) { // p[0]=0;
			if (p[i] != i) {
				boost::add_edge(i, p[i], get(boost::edge_weight_t(), input_graph, boost::edge(i, p[i], input_graph).first), output_graph); // add edge
			}
		}
		boost::remove_edge(0, v1, output_graph); // remove edge (0,v1)
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count(); // Nanosecond
	MST_time = runningtime / 1e6;

	return output_graph;
}

#pragma endregion find_MST_Yahui 2016年12月7日


#pragma region 

void MST_P3_Yahui(graph input_graph, graph& solu_graph, double& MST_time) {

	if (num_edges(solu_graph) <= 1) {
		MST_time = 0;
	}
	else {
		auto begin_time = std::chrono::high_resolution_clock::now(); // start time

		int N = num_vertices(input_graph);
		//graph base_graph = input_graph;
		//for (int i = 0; i < N; i++) {
		//	if (in_degree(i, solu_graph) == 0) {
		//		clear_vertex(i, base_graph); // remove unrelated edges; this is much slower than the version below
		//	}
		//}
		graph base_graph;
		for (int i = 0; i < N; i++) {
			boost::add_vertex(i, base_graph);
			boost::put(boost::vertex_name_t(), base_graph, i, get(boost::vertex_name_t(), input_graph, i)); // input node weight
		}
		graph::out_edge_iterator eit, eend;
		for (int i = 0; i < N; i++) {
			if (in_degree(i, solu_graph) >= 1) { // i is in solu_graph
				tie(eit, eend) = boost::out_edges(i, input_graph); // adjacent_vertices of i in input_graph
				for_each(eit, eend,
					[&input_graph, &i, &base_graph, &solu_graph](graph::edge_descriptor it)
				{
					int j = boost::target(it, input_graph);
					if (i > j && in_degree(j, solu_graph) > 0) { // j is in solu_graph
						boost::add_edge(i, j, get(boost::edge_weight_t(), input_graph, boost::edge(i, j, input_graph).first), base_graph); // input edge 
					}
				});
			}
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count(); // Nanosecond
		double time1 = runningtime / 1e6;
		double time2;
		base_graph = find_MST_Yahui(base_graph, time2);
		MST_time = time1 + time2;
		solu_graph = base_graph;
	}
}

#pragma endregion MST_P3_Yahui 


#pragma region 
void Growing_P3(graph initial_graph, graph& solu_graph, double& Growing_time) {

	// the details of the growing algorithm for P3
	// basic veriables:
	// vector<int> UV: the included vertices that may be the root to a path candidate
	// vector<bool> included: the vertices included in solu_graph
	// process:
	// 1. construct vector<bool> included
	// 2. construct vector<int> UV
	// 3. check vertices in UV
	// 4. for each checked vertex, find a path candidate with at most n new vertices
	// 5. add the path candidate; add new vertices into vector<int> UV
	// 6. if there is no path candidate for a vertex in UV; remove it from UV

	vector<int> UV;
	int N = num_vertices(initial_graph);
	vector<bool> included(N);

	// initialize UV, included
	if (num_edges(solu_graph) > 0) { // there is at least an edge in solu_graph
		for (int i = 0; i < N; i++) {
			if (in_degree(i, solu_graph) > 0) {
				included[i] = true; // i is included
				if (in_degree(i, solu_graph) < in_degree(i, initial_graph)) {
					UV.insert(UV.end(), i); // i is connected to non-included vertices
				}
			}
		}
	}
	else { // there is no edge in solu_graph
		int max_v = 0;
		double max = get(boost::vertex_name_t(), initial_graph, 0);
		for (int i = 1; i < N; i++) {
			if (get(boost::vertex_name_t(), initial_graph, i) > max) {
				max = get(boost::vertex_name_t(), initial_graph, i);
				max_v = i;
			}
		}
		included[max_v] = true; // max_v is included
		UV.insert(UV.end(), max_v); // max_v is connected to non-included vertices
	}

	auto begin_time = std::chrono::high_resolution_clock::now(); // start time

	while (UV.size() > 0) { // while there is an unchecked vertex
		int target1 = UV[0]; // the root vertex
							 // find a path candidate with at most 2 new vertices
		bool end = false;
		typedef boost::graph_traits<graph>::adjacency_iterator AdjacencyIterator;
		AdjacencyIterator ai1, a_end1;
		boost::tie(ai1, a_end1) = boost::adjacent_vertices(target1, initial_graph);
		for (; ai1 != a_end1; ai1++) {
			int target2 = *ai1;
			if (included[target2] == false) {
				double net_prize1 = get(boost::vertex_name_t(), initial_graph, target2) -
					get(boost::edge_weight_t(), initial_graph, boost::edge(target1, target2, initial_graph).first);
				if (net_prize1 > 0) {
					end = true;
					// add the path candidate: target2
					//cout << "target2: " << target2 << endl;
					//cout << "target2: " << get(boost::vertex_name_t(), initial_graph, target2) << endl;
					//cout << "net_prize1: " << net_prize1 << endl;
					boost::add_edge(target1, target2, get(boost::edge_weight_t(), initial_graph, boost::edge(target1, target2, initial_graph).first), solu_graph);
					included[target2] = true;
					if (in_degree(target2, solu_graph) < in_degree(target2, initial_graph)) {
						UV.insert(UV.end(), target2);
					}
				}
				if (end == true) {
					break;
				}
				else { // find target3
					typedef boost::graph_traits<graph>::adjacency_iterator AdjacencyIterator;
					AdjacencyIterator ai2, a_end2;
					boost::tie(ai2, a_end2) = boost::adjacent_vertices(target2, initial_graph);
					for (; ai2 != a_end2; ai2++) {
						int target3 = *ai2;
						if (included[target3] == false) {
							double net_prize2 = get(boost::vertex_name_t(), initial_graph, target3) -
								get(boost::edge_weight_t(), initial_graph, boost::edge(target2, target3, initial_graph).first) + net_prize1;
							if (net_prize2 > 0) {
								end = true;
								// add the path candidate: target2, target3
								boost::add_edge(target1, target2, get(boost::edge_weight_t(), initial_graph, boost::edge(target1, target2, initial_graph).first), solu_graph);
								boost::add_edge(target2, target3, get(boost::edge_weight_t(), initial_graph, boost::edge(target2, target3, initial_graph).first), solu_graph);
								included[target2] = true;
								included[target3] = true;
								if (in_degree(target2, solu_graph) < in_degree(target2, initial_graph)) {
									UV.insert(UV.end(), target2);
								}
								if (in_degree(target3, solu_graph) < in_degree(target3, initial_graph)) {
									UV.insert(UV.end(), target3);
								}
							}
							if (end == true) {
								break;
							}
						}
					}
				}
			}
		}
		if (end == false) { // there is no path candidate for target1
			UV.erase(UV.begin()); // erase UV[0]
		}
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count(); // Nanosecond
	Growing_time = runningtime / 1e6;

}
#pragma endregion Growing_P3


#pragma region 

void P3(graph initial_graph, graph& solu_graph, double& P3_time) {

	bool changed = true;

	while (changed == true) {
		double old_cost = net_cost(solu_graph);
		//cout << "old_cost: " << net_cost(solu_graph) << endl;
		// TGA
		double Growing_time;
		Growing_P3(initial_graph, solu_graph, Growing_time);
		P3_time = P3_time + Growing_time;
		//cout << "TGA: " << net_cost(solu_graph) << endl;
		// MST
		double MST_time;
		MST_P3_Yahui(initial_graph, solu_graph, MST_time);
		P3_time = P3_time + MST_time;
		//cout << "MST: " << net_cost(solu_graph) << endl;
		// GPrA
		double GPrA_time;
		graph GPrA_graph = GPrA(solu_graph, GPrA_time);
		P3_time = P3_time + GPrA_time;
		//cout << "GPrA: " << net_cost(solu_graph) << endl;
		double new_cost = net_cost(solu_graph);
		if (new_cost == old_cost) {
			changed = false;
		}
	}
	cout << "P3_time: " << P3_time << "ms" << endl;
}

#pragma endregion P3 


#pragma region
void P3_staynerd_dimacs_hand() {
	string instance_name;
	ofstream outputFile;
	outputFile.precision(4);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);
	outputFile.open("P3_Result_staynerd_dimacs_hand.csv");
	outputFile << "P3_Result_staynerd_dimacs_hand.csv" << endl;
	outputFile << "Instance,Old Cost,New Cost,Time" << endl;

	for (int name = 1; name < 49; name++) {
		// instance_name
		if (name < 10) {
			instance_name = "handbd" + to_string(0) + to_string(name);
		}
		else if (name < 15) {
			instance_name = "handbd" + to_string(name);
		}
		else if (name < 24) {
			instance_name = "handbi" + to_string(0) + to_string(name - 14);
		}
		else if (name < 29) {
			instance_name = "handbi" + to_string(name - 14);
		}
		else if (name < 38) {
			instance_name = "handsd" + to_string(0) + to_string(name - 28);
		}
		else if (name == 38) {
			instance_name = "handsd" + to_string(name - 28);
		}
		else if (name < 48) {
			instance_name = "handsi" + to_string(0) + to_string(name - 38);
		}
		else if (name == 48) {
			instance_name = "handsi" + to_string(name - 38);
		}
		cout << instance_name << ":";

		graph input_graph = read_data_with_terminals(instance_name+".stp");
		graph staynerd_graph = read_data_staynerd_dimacs(input_graph, "stp." + instance_name + ".stpsolver.4.3600.sol");

		double old_cost = net_cost(staynerd_graph);
		cout << "old_cost= " << old_cost << endl;
		double P3_time = 0;
		P3(input_graph, staynerd_graph, P3_time);
		double new_cost = net_cost(staynerd_graph);
		cout << "new_cost= " << new_cost << endl;

		if (new_cost + 1e-5 < old_cost) {
			outputFile << instance_name << "," << old_cost << "," << new_cost << ","
				<< P3_time / 1000 << "s" << endl;
		}
	}
	outputFile.close();
	std::cout << "END" << endl;
}
#pragma endregion P3_staynerd_dimacs_hand


#pragma region
void P3_bbdualascent() {
	string instance_name;
	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);
	outputFile.open("P3_Result_bbdualascent_hand.csv");
	outputFile << "P3_Result_bbdualascent_hand.csv" << endl;
	outputFile << "Instance,Old Cost,New Cost,Time" << endl;

	for (int name = 1; name < 49; name++) {
		// instance_name
		if (name < 10) {
			instance_name = "handbd" + to_string(0) + to_string(name);
		}
		else if (name < 15) {
			instance_name = "handbd" + to_string(name);
		}
		else if (name < 24) {
			instance_name = "handbi" + to_string(0) + to_string(name - 14);
		}
		else if (name < 29) {
			instance_name = "handbi" + to_string(name - 14);
		}
		else if (name < 38) {
			instance_name = "handsd" + to_string(0) + to_string(name - 28);
		}
		else if (name == 38) {
			instance_name = "handsd" + to_string(name - 28);
		}
		else if (name < 48) {
			instance_name = "handsi" + to_string(0) + to_string(name - 38);
		}
		else if (name == 48) {
			instance_name = "handsi" + to_string(name - 38);
		}
		cout << instance_name << ":";

		graph input_graph = read_data_with_terminals(instance_name + ".stp");
		graph staynerd_graph = read_data_staynerd_dimacs(input_graph, "bbdualascent_" + instance_name + ".sol");

		double old_cost = net_cost(staynerd_graph);
		cout << "old_cost= " << old_cost << endl;
		double P3_time = 0;
		P3(input_graph, staynerd_graph, P3_time);
		double new_cost = net_cost(staynerd_graph);
		cout << "new_cost= " << new_cost << endl;

		if (new_cost + 1e-5 < old_cost) {
			outputFile << instance_name << "," << old_cost << "," << new_cost << ","
				<< P3_time / 1000 << "s" << endl;
		}
	}
	outputFile.close();
	std::cout << "END" << endl;
}
#pragma endregion P3_bbdualascent


int main()
{
	std::cout << std::setprecision(4) << std::fixed;
	P3_staynerd_dimacs_hand();
	P3_bbdualascent();
	getchar();
}

