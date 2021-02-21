// Copyright Michael Drexl 2005, 2006.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://boost.org/LICENSE_1_0.txt)

// Example use of the resource-constrained shortest paths algorithm.
#include <boost/config.hpp>

#ifdef BOOST_MSVC
#pragma warning(disable : 4267)
#endif

#include <boost/graph/adjacency_list.hpp>

#include <boost/graph/r_c_shortest_paths.hpp>
#include <iostream>
#include <map>

using namespace boost;

struct VertProp 
{
    VertProp(int num=0, double min_dist=0, double max_dist=0)
    : num(num), min_distance(min_dist), max_distance(max_dist) { }
    int num;				// id
    double min_distance;	// minimal traversal distance
    double max_distance;	// maximal traversal distance
};

struct EdgeProp
{
    EdgeProp(int n, double t, double d) : num(n), time(t), distance(d)  { }
				
    int num; // id
    double time; // traversal time
    double distance; // traversal distance
};

typedef adjacency_list<vecS, vecS, directedS, VertProp, EdgeProp> BoostGraph;
typedef graph_traits<BoostGraph>::edge_descriptor Edge;
typedef graph_traits<BoostGraph>::vertex_descriptor Vertex;

// data structures for shortest path problem with distance windows (spptw)
// ResourceContainer model 
struct ResourceContainer
{
    ResourceContainer(int c = 0, int t = 0) : time(c), distance(t) {}
    ResourceContainer& operator=(const ResourceContainer& other)
    {
        if (this == &other)
            return *this;
		
        this->time = other.time;
		this->distance = other.distance;

		return *this;
    }
    double time;
    double distance;
};

bool operator==(const ResourceContainer& res_cont_1, const ResourceContainer& res_cont_2)
{
    return (res_cont_1.time == res_cont_2.time && res_cont_1.distance == res_cont_2.distance);
}

bool operator<(const ResourceContainer& res_cont_1, const ResourceContainer& res_cont_2)
{
    if (res_cont_1.time > res_cont_2.time)
        return false;
    if (res_cont_1.time == res_cont_2.time)
        return res_cont_1.distance < res_cont_2.distance;
    return true;
}

// ResourceExtensionFunction model
class ResourceExtensionFunction
{
public:
    inline bool operator()(const BoostGraph& g, ResourceContainer& new_cont, const ResourceContainer& old_cont, Edge ed) const
    {
        const EdgeProp& arc_prop = get(edge_bundle, g)[ed];
        const VertProp& vert_prop = get(vertex_bundle, g)[target(ed, g)];
        new_cont.time = old_cont.time + arc_prop.time;
        double& i_distance = new_cont.distance;
        i_distance = old_cont.distance + arc_prop.distance;
        i_distance < vert_prop.min_distance ? i_distance = vert_prop.min_distance : 0;
        return i_distance <= vert_prop.max_distance ? true : false;
    }
};

// DominanceFunction model
class DominanceFunction
{
public:
    inline bool operator()(const ResourceContainer& res_cont_1,
        const ResourceContainer& res_cont_2) const
    {
        // must be "<=" here, not "<"
        return res_cont_1.time <= res_cont_2.time
            && res_cont_1.distance <= res_cont_2.distance;
    }
};
// end data structures for shortest path problem with distance windows (spptw)


int main()
{
    BoostGraph g;

	std::vector<Vertex> vertices;
	std::vector<Edge> edges;
	std::vector<EdgeProp> edge_props;
	std::vector<double> min_dist(5, 0.0);
	std::vector<double> max_dist(5, 100.0);
	std::vector<double> time(6, 1.0); 
	std::vector<double> distance(6, 1.0);
	std::map<std::pair<int,int>,int> edge_map;
	
    vertices.push_back(add_vertex(VertProp(0, min_dist[0], max_dist[0]), g));
    vertices.push_back(add_vertex(VertProp(1, min_dist[1], max_dist[1]), g));
    vertices.push_back(add_vertex(VertProp(2, min_dist[2], max_dist[2]), g));
    vertices.push_back(add_vertex(VertProp(3, min_dist[3], max_dist[3]), g));
	vertices.push_back(add_vertex(VertProp(4, min_dist[4], max_dist[4]), g));

	time[3] = 2.5;
	distance[3] = 0;

	time[5] = 3.5;
	distance[5] = 2.5;
	
	edge_props.push_back(EdgeProp(0, time[0], distance[0]));
	edge_props.push_back(EdgeProp(1, time[1], distance[1]));
	edge_props.push_back(EdgeProp(2, time[2], distance[2]));
	edge_props.push_back(EdgeProp(3, time[3], distance[3]));
	edge_props.push_back(EdgeProp(4, time[4], distance[4]));
	edge_props.push_back(EdgeProp(5, time[5], distance[5]));
	
    edges.push_back(add_edge(0, 1, edge_props[0], g).first);
    edges.push_back(add_edge(1, 2, edge_props[1], g).first);
    edges.push_back(add_edge(1, 3, edge_props[2], g).first);
    edges.push_back(add_edge(2, 4, edge_props[3], g).first);
    edges.push_back(add_edge(3, 4, edge_props[4], g).first);
	edges.push_back(add_edge(0, 4, edge_props[5], g).first);

	edge_map[std::make_pair(0,1)] = 0;
	edge_map[std::make_pair(1,2)] = 1;
	edge_map[std::make_pair(1,3)] = 2;
	edge_map[std::make_pair(2,4)] = 3;
	edge_map[std::make_pair(3,4)] = 4;
	edge_map[std::make_pair(0,4)] = 5;
	
    Vertex s = vertices[0];
    Vertex t = vertices[4];

    // spptw
    std::vector<std::vector<Edge>> opt_solutions;
    std::vector<ResourceContainer> pareto_opt_rcs;

    r_c_shortest_paths(g, get(&VertProp::num, g),
        get(&EdgeProp::num, g), s, t, opt_solutions,
        pareto_opt_rcs, ResourceContainer(0, 0), ResourceExtensionFunction(),
        DominanceFunction(),
        std::allocator<r_c_shortest_paths_label<BoostGraph, ResourceContainer>>(),
        default_r_c_shortest_paths_visitor());

    std::cout << "SPP with distance windows:" << std::endl;
    std::cout << "Number of optimal solutions: ";
    std::cout << static_cast<int>(opt_solutions.size()) << std::endl;
    for (int i = 0; i < static_cast< int >(opt_solutions.size()); ++i)
    {
        std::cout << "The " << i << "th shortest path from A to E is: ";
        std::cout << std::endl;
        for (int j = static_cast<int>(opt_solutions[i].size()) - 1;
             j >= 0; --j)
            std::cout << edge_map[std::make_pair(source(opt_solutions[i][j], g), target(opt_solutions[i][j], g))]
                      << std::endl;
        std::cout << "time: " << pareto_opt_rcs[i].time << std::endl;
        std::cout << "distance: " << pareto_opt_rcs[i].distance << std::endl;
    }

	opt_solutions.clear();
	pareto_opt_rcs.clear();

	g[4].max_distance = 2.5;
    return 0;
}
