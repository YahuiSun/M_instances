#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <chrono>
#include <boost/config.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
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
		int v1 = 0;
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
#pragma endregion GPrA     2016年11月26日18:06


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
void save_data_MST_GPrA(string instance_name, graph result_graph, double net_cost, double running_time) {

	string save_name = "Result_MST_GPrA_" + instance_name; // save_name
	ofstream outputFile;
	outputFile.precision(4);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);
	outputFile.open(save_name + ".stp"); // stp file

	outputFile << "33D32945 STP File, STP Format Version 1.0" << endl;
	outputFile << endl;

	// comments
	outputFile << "SECTION Comments" << endl;
	outputFile << "Name \"" << save_name << "\"" << endl;
	outputFile << "Net_cost " << net_cost << endl;
	outputFile << "Running_time " << running_time << "ms" << endl;
	outputFile << "Creator \"Yahui Sun\"" << endl;
	outputFile << "Problem \"Prize - Collecting Steiner Problem in Graphs\"" << endl;
	outputFile << "END" << endl;
	outputFile << endl;

	// graph
	outputFile << "SECTION Graph" << endl;
	outputFile << "Nodes " << num_vertices(result_graph) << endl;
	outputFile << "Edges " << num_edges(result_graph) << endl;
	graph::out_edge_iterator eit, eend;
	for (int i = 0; i < num_vertices(result_graph); i++) {
		tie(eit, eend) = boost::out_edges(i, result_graph); // adjacent_vertices of 2
		for_each(eit, eend,
			[&result_graph, &i, &outputFile](graph::edge_descriptor it)
		{
			int j = boost::target(it, result_graph);
			if (i < j) {
				outputFile << "E " << i + 1 << " " << j + 1 << " " << get(boost::edge_weight_t(), result_graph, boost::edge(i, j, result_graph).first) << endl;
			}
		});
	}
	outputFile << "END" << endl;
	outputFile << endl;

	// TP
	outputFile << "SECTION Non-Compulsory Terminals" << endl;
	int p = 0;
	for (int i = 0; i < num_vertices(result_graph); i++) {
		if (get(boost::vertex_name_t(), result_graph, i) > 0) {
			p++;
		}
	}
	outputFile << "Terminals " << p << endl;
	for (int i = 0; i < num_vertices(result_graph); i++) {
		if (get(boost::vertex_name_t(), result_graph, i) > 0) {
			outputFile << "TP " << i + 1 << " " << get(boost::vertex_name_t(), result_graph, i) << endl;
		}
	}
	outputFile << "END" << endl;
	outputFile << endl;

	outputFile << "EOF" << endl;

}
#pragma endregion save_data_MST_GPrA  2016年11月26日18:36


#pragma region
void solve_hand_MST_GPrA() {

	string instance_name; // instance_name
	string instance;
	ofstream outputFile;
	outputFile.precision(4);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);
	double MST_time;
	double GPrA_time;
	double Total_time;
	double GW_lb;
	double GAP;
	double min_time;
	double times = 10; // run each instance for x times
	graph input_graph;
	graph GPrA_graph;
	graph MST_graph;

	// write file
	outputFile.open("Result_MST_GPrA_hand.csv");
	outputFile << "Result_MST_GPrA_hand.csv" << endl;
	outputFile << "Instance,|V|,|E|,Time,Cost" << endl;

	for (int name = 1; name < 49; name++) {

		// instance_name
		if (name < 10) {
			instance = "handbd";
			instance_name = instance + to_string(0) + to_string(name);
		}
		else if (name < 15) {
			instance = "handbd";
			instance_name = instance + to_string(name);
		}
		else if (name < 24) {
			instance = "handbi";
			instance_name = instance + to_string(0) + to_string(name - 14);
		}
		else if (name < 29) {
			instance = "handbi";
			instance_name = instance + to_string(name - 14);
		}
		else if (name < 38) {
			instance = "handsd";
			instance_name = instance + to_string(0) + to_string(name - 28);
		}
		else if (name == 38) {
			instance = "handsd";
			instance_name = instance + to_string(name - 28);
		}
		else if (name < 48) {
			instance = "handsi";
			instance_name = instance + to_string(0) + to_string(name - 38);
		}
		else if (name == 48) {
			instance = "handsi";
			instance_name = instance + to_string(name - 38);
		}

		input_graph = read_data_with_terminals(instance_name + ".stp");
		min_time = 1e6;
		double max_time = 1e-6;
		double avg_time = 0;

		for (int i = 0; i < times; i++) {
			MST_graph = find_MST_Yahui(input_graph, MST_time);
			std::cout << "MST_time: " << MST_time << "ms" << std::endl;
			GPrA_graph = GPrA(MST_graph, GPrA_time);
			std::cout << "GPrA_time: " << GPrA_time << "ms" << std::endl;
			double Total_time = MST_time + GPrA_time;
			std::cout << "Total_time: " << Total_time << "ms" << std::endl;
			avg_time = avg_time + Total_time;
		}
		avg_time = avg_time / times;
		//save_data
		double cost = net_cost(GPrA_graph);
		cout << "cost= " << cost << endl;
		cout << "avg_time= " << avg_time << "ms" << endl;
		save_data_MST_GPrA(instance_name, GPrA_graph, cost, avg_time);

		outputFile << instance_name << "," << num_vertices(input_graph) << "," << num_edges(input_graph) << ","
			<< avg_time/1000 << "s" << "," << net_cost(GPrA_graph) << endl;
	}

	outputFile.close();
	std::cout << "END" << endl;
}
#pragma endregion solve_hand_MST_GPrA


#pragma region
void solve_PPI_MST_GPrA() {

	string instance_name; // instance_name
	string instance;
	ofstream outputFile;
	outputFile.precision(4);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);
	double MST_time;
	double GPrA_time;
	double Total_time;
	double GW_lb;
	double GAP;
	double min_time;
	double times = 10; // run each instance for x times
	graph input_graph;
	graph GPrA_graph;
	graph MST_graph;

	// write file
	outputFile.open("Result_MST_GPrA_PPI.csv");
	outputFile << "Result_MST_GPrA_PPI.csv" << endl;
	outputFile << "Instance,|V|,|E|,Time,Cost" << endl;

	for (int name = 1; name < 11; name++) {

		// instance_name
		instance = "PPI";
		instance_name = instance + "_" + to_string(name);

		input_graph = read_data_with_terminals(instance_name + ".stp");
		double distribution_ratio = 2; // ep slack = 1:(distribution_ratio-1)
		min_time = 1e6;
		double max_time = 1e-6;
		double avg_time = 0;

		for (int i = 0; i < times; i++) {
			MST_graph = find_MST_Yahui(input_graph, MST_time);
			std::cout << "MST_time: " << MST_time << "ms" << std::endl;
			GPrA_graph = GPrA(MST_graph, GPrA_time);
			std::cout << "GPrA_time: " << GPrA_time << "ms" << std::endl;
			double Total_time = MST_time + GPrA_time;
			std::cout << "Total_time: " << Total_time << "ms" << std::endl;
			if (min_time > Total_time) {
				min_time = Total_time;
			}
			if (max_time < Total_time) {
				max_time = Total_time;
			}
			avg_time = avg_time + Total_time;
		}
		avg_time = (avg_time - min_time - max_time) / 8;
		//save_data
		double cost = net_cost(GPrA_graph);
		cout << "cost= " << cost << endl;
		cout << "avg_time= " << avg_time << "ms" << endl;
		save_data_MST_GPrA(instance_name, GPrA_graph, cost, avg_time);

		outputFile << instance_name << "," << num_vertices(input_graph) << "," << num_edges(input_graph) << ","
			<< avg_time / 1000 << "s" << "," << net_cost(GPrA_graph) << endl;
	}

	outputFile.close();
	std::cout << "END" << endl;
}
#pragma endregion solve_PPI_MST_GPrA


#pragma region
void solve_DSN_MST_GPrA() {

	string instance_name; // instance_name
	string instance;
	ofstream outputFile;
	outputFile.precision(4);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);
	double MST_time;
	double GPrA_time;
	double Total_time;
	double GW_lb;
	double GAP;
	double min_time;
	double times = 10; // run each instance for x times
	graph input_graph;
	graph GPrA_graph;
	graph MST_graph;

	// write file
	outputFile.open("Result_MST_GPrA_DSN.csv");
	outputFile << "Result_MST_GPrA_DSN.csv" << endl;
	outputFile << "Instance,|V|,|E|,Time,Cost" << endl;

	for (int name = 1; name < 21; name++) {

		if (name != 6 && name != 16) {
			// instance_name
			instance = "D";
			if (name < 10) {
				instance_name = instance + "_0" + to_string(name) + "_a";
			}
			else if (name == 10) {
				instance_name = instance + "_" + to_string(name) + "_a";
			}
			else if (name < 20) {
				instance_name = instance + "_0" + to_string(name - 10) + "_b";
			}
			else if (name == 20) {
				instance_name = instance + "_" + to_string(name - 10) + "_b";
			}

			input_graph = read_data_with_terminals(instance_name + ".stp");
			min_time = 1e6;
			double max_time = 1e-6;
			double avg_time = 0;

			for (int i = 0; i < times; i++) {
				MST_graph = find_MST_Yahui(input_graph, MST_time);
				std::cout << "MST_time: " << MST_time << "ms" << std::endl;
				GPrA_graph = GPrA(MST_graph, GPrA_time);
				std::cout << "GPrA_time: " << GPrA_time << "ms" << std::endl;
				double Total_time = MST_time + GPrA_time;
				std::cout << "Total_time: " << Total_time << "ms" << std::endl;
				if (min_time > Total_time) {
					min_time = Total_time;
				}
				if (max_time < Total_time) {
					max_time = Total_time;
				}
				avg_time = avg_time + Total_time;
			}
			avg_time = (avg_time - min_time - max_time) / 8;
			//save_data
			double cost = net_cost(GPrA_graph);
			cout << "cost= " << cost << endl;
			cout << "avg_time= " << avg_time << "ms" << endl;
			save_data_MST_GPrA(instance_name, GPrA_graph, cost, avg_time);

			outputFile << instance_name << "," << num_vertices(input_graph) << "," << num_edges(input_graph) << ","
				<< avg_time / 1000 << "s" << "," << net_cost(GPrA_graph) << endl;
		}
	}

	outputFile.close();
	std::cout << "END" << endl;
}
#pragma endregion solve_DSN_MST_GPrA


#pragma region
void solve_M_MST_GPrA() {
	// write file
	string instance = "M";
	ofstream outputFile;
	outputFile.precision(4);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);
	outputFile.open("Result_MST_GPrA_M.csv");
	outputFile << "Result_MST_GPrA_M.csv" << endl;
	outputFile << "Instance,|V|,|E|,Time,Cost" << endl;

	for (int name = 21; name < 41; name++) {

		string instance_name; // instance_name
		if (name < 21) {
			instance_name = instance + to_string(name) + "A";
		}
		else {
			instance_name = instance + to_string(name - 20) + "B";
		}
		graph input_graph = read_data_with_terminals(instance_name + ".stp");

		double MST_time;
		double GPrA_time;
		double Total_time;
		double min_time = 1e6;
		graph MST_graph;
		graph GPrA_graph;
		double times = 1;
		double max_time = 1e-6;
		double avg_time = 0;

		for (int i = 0; i < times; i++) {
			MST_graph = find_MST_Yahui(input_graph, MST_time);
			std::cout << "MST_time: " << MST_time << "ms" << std::endl;
			GPrA_graph = GPrA(MST_graph, GPrA_time);
			std::cout << "GPrA_time: " << GPrA_time << "ms" << std::endl;
			double Total_time = MST_time + GPrA_time;
			std::cout << "Total_time: " << Total_time << "ms" << std::endl;
			avg_time = avg_time + Total_time;
		}
		avg_time = avg_time / times;
		//save_data
		double cost = net_cost(GPrA_graph);
		cout << "cost= " << cost << endl;
		cout << "avg_time= " << avg_time << "ms" << endl;
		save_data_MST_GPrA(instance_name, GPrA_graph, cost, avg_time);

		outputFile << instance_name << "," << num_vertices(input_graph) << "," << num_edges(input_graph) << ","
			<< avg_time / 1000 << "s" << "," << net_cost(GPrA_graph) << endl;
	}

	outputFile.close();
	std::cout << "END" << endl;
}
#pragma endregion solve_M_MST_GPrA


#pragma region 
void main()
{
	std::cout << std::setprecision(4) << std::fixed;

	solve_hand_MST_GPrA();

	getchar();
}
#pragma endregion main