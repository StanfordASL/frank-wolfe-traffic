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

using namespace boost;

struct VertProp 
{
    VertProp(int n = 0, int e = 0, int l = 0)
    : num(n), eat(e), lat(l) { }
    int num; // id
    int eat; // earliest arrival time
    int lat; // latest arrival time
};

struct ArcProp
{
    ArcProp(int n = 0, int c = 0, int t = 0)
    : num(n), cost(c), time(t)  { }
    int num; // id
    int cost; // traversal cost
    int time; // traversal time
};

typedef adjacency_list< vecS, vecS, directedS, VertProp, ArcProp >
    BoostGraph;

// data structures for shortest path problem with time windows (spptw)
// ResourceContainer model 
struct ResourceContainer
{
    ResourceContainer(int c = 0, int t = 0) : cost(c), time(t) {}
    ResourceContainer& operator=(const ResourceContainer& other)
    {
        if (this == &other)
            return *this;
		
        this->cost = other.cost;
		this->time = other.time;

		return *this;
    }
    int cost;
    int time;
};

bool operator==(const ResourceContainer& res_cont_1, const ResourceContainer& res_cont_2)
{
    return (res_cont_1.cost == res_cont_2.cost && res_cont_1.time == res_cont_2.time);
}

bool operator<(const ResourceContainer& res_cont_1, const ResourceContainer& res_cont_2)
{
    if (res_cont_1.cost > res_cont_2.cost)
        return false;
    if (res_cont_1.cost == res_cont_2.cost)
        return res_cont_1.time < res_cont_2.time;
    return true;
}

// ResourceExtensionFunction model
class ResourceExtensionFunction
{
public:
    inline bool operator()(const BoostGraph& g, ResourceContainer& new_cont, const ResourceContainer& old_cont, graph_traits< BoostGraph >::edge_descriptor ed) const
    {
        const ArcProp& arc_prop = get(edge_bundle, g)[ed];
        const VertProp& vert_prop = get(vertex_bundle, g)[target(ed, g)];
        new_cont.cost = old_cont.cost + arc_prop.cost;
        int& i_time = new_cont.time;
        i_time = old_cont.time + arc_prop.time;
        i_time < vert_prop.eat ? i_time = vert_prop.eat : 0;
        return i_time <= vert_prop.lat ? true : false;
    }
};

// DominanceFunction model
class DominanceFunction
{
public:
    inline bool operator()(const ResourceContainer& res_cont_1,
        const ResourceContainer& res_cont_2) const
    {
        // must be "<=" here!!!
        // must NOT be "<"!!!
        return res_cont_1.cost <= res_cont_2.cost
            && res_cont_1.time <= res_cont_2.time;
        // this is not a contradiction to the documentation
        // the documentation says:
        // "A label $l_1$ dominates a label $l_2$ if and only if both are
        // resident at the same vertex, and if, for each resource, the resource
        // consumption of $l_1$ is less than or equal to the resource
        // consumption of $l_2$, and if there is at least one resource where
        // $l_1$ has a lower resource consumption than $l_2$." one can think of
        // a new label with a resource consumption equal to that of an old label
        // as being dominated by that old label, because the new one will have a
        // higher number and is created at a later point in time, so one can
        // implicitly use the number or the creation time as a resource for
        // tie-breaking
    }
};
// end data structures for shortest path problem with time windows (spptw)

// example graph structure and cost from
// http://www.boost.org/libs/graph/example/dijkstra-example.cpp
enum nodes
{
    A,
    B,
    C,
    D,
    E
};
char name[] = "ABCDE";

int main()
{
    BoostGraph g;

    add_vertex(VertProp(A, 0, 0), g);
    add_vertex(VertProp(B, 5, 20), g);
    add_vertex(VertProp(C, 6, 10), g);
    add_vertex(VertProp(D, 3, 12), g);
    add_vertex(VertProp(E, 0, 100), g);

    add_edge(A, C, ArcProp(0, 1, 5), g);
    add_edge(B, B, ArcProp(1, 2, 5), g);
    add_edge(B, D, ArcProp(2, 1, 2), g);
    add_edge(B, E, ArcProp(3, 2, 7), g);
    add_edge(C, B, ArcProp(4, 7, 3), g);
    add_edge(C, D, ArcProp(5, 3, 8), g);
    add_edge(D, E, ArcProp(6, 1, 3), g);
    add_edge(E, A, ArcProp(7, 1, 5), g);
    add_edge(E, B, ArcProp(8, 1, 4), g);

    // the unique shortest path from A to E in the dijkstra-example.cpp is
    // A -> C -> D -> E
    // its length is 5
    // the following code also yields this result

    // with the above time windows, this path is infeasible
    // now, there are two shortest paths that are also feasible with respect to
    // the vertex time windows:
    // A -> C -> B -> D -> E and
    // A -> C -> B -> E
    // however, the latter has a longer total travel time and is therefore not
    // pareto-optimal, i.e., it is dominated by the former path
    // therefore, the code below returns only the former path

    // spp without resource constraints
    graph_traits< BoostGraph >::vertex_descriptor s = A;
    graph_traits< BoostGraph >::vertex_descriptor t = E;

    std::vector<
        std::vector< graph_traits< BoostGraph >::edge_descriptor > >
        opt_solutions;
    std::vector< spp_no_rc_res_cont > pareto_opt_rcs_no_rc;

    r_c_shortest_paths(g, get(&VertProp::num, g),
        get(&ArcProp::num, g), s, t, opt_solutions,
        pareto_opt_rcs_no_rc, spp_no_rc_res_cont(0), ref_no_res_cont(),
        dominance_no_res_cont(),
        std::allocator< r_c_shortest_paths_label< BoostGraph,
            spp_no_rc_res_cont > >(),
        default_r_c_shortest_paths_visitor());

    std::cout << "SPP without resource constraints:" << std::endl;
    std::cout << "Number of optimal solutions: ";
    std::cout << static_cast< int >(opt_solutions.size()) << std::endl;
    for (int i = 0; i < static_cast< int >(opt_solutions.size()); ++i)
    {
        std::cout << "The " << i << "th shortest path from A to E is: ";
        std::cout << std::endl;
        for (int j = static_cast< int >(opt_solutions[i].size()) - 1; j >= 0;
             --j)
            std::cout << name[source(opt_solutions[i][j], g)] << std::endl;
        std::cout << "E" << std::endl;
        std::cout << "Length: " << pareto_opt_rcs_no_rc[i].cost << std::endl;
    }
    std::cout << std::endl;

    // spptw
    std::vector<
        std::vector< graph_traits< BoostGraph >::edge_descriptor > >
        opt_solutions_spptw;
    std::vector< ResourceContainer > pareto_opt_rcs_spptw;

    r_c_shortest_paths(g, get(&VertProp::num, g),
        get(&ArcProp::num, g), s, t, opt_solutions_spptw,
        pareto_opt_rcs_spptw, ResourceContainer(0, 0), ResourceExtensionFunction(),
        DominanceFunction(),
        std::allocator< r_c_shortest_paths_label< BoostGraph,
            ResourceContainer > >(),
        default_r_c_shortest_paths_visitor());

    std::cout << "SPP with time windows:" << std::endl;
    std::cout << "Number of optimal solutions: ";
    std::cout << static_cast< int >(opt_solutions.size()) << std::endl;
    for (int i = 0; i < static_cast< int >(opt_solutions.size()); ++i)
    {
        std::cout << "The " << i << "th shortest path from A to E is: ";
        std::cout << std::endl;
        for (int j = static_cast< int >(opt_solutions_spptw[i].size()) - 1;
             j >= 0; --j)
            std::cout << name[source(opt_solutions_spptw[i][j], g)]
                      << std::endl;
        std::cout << "E" << std::endl;
        std::cout << "Length: " << pareto_opt_rcs_spptw[i].cost << std::endl;
        std::cout << "Time: " << pareto_opt_rcs_spptw[i].time << std::endl;
    }

    // utility function check_r_c_path example
    std::cout << std::endl;
    bool b_is_a_path_at_all = false;
    bool b_feasible = false;
    bool b_correctly_extended = false;
    ResourceContainer actual_final_resource_levels(0, 0);
    graph_traits< BoostGraph >::edge_descriptor ed_last_extended_arc;
    check_r_c_path(g, opt_solutions_spptw[0], ResourceContainer(0, 0), true,
        pareto_opt_rcs_spptw[0], actual_final_resource_levels, ResourceExtensionFunction(),
        b_is_a_path_at_all, b_feasible, b_correctly_extended,
        ed_last_extended_arc);
    if (!b_is_a_path_at_all)
        std::cout << "Not a path." << std::endl;
    if (!b_feasible)
        std::cout << "Not a feasible path." << std::endl;
    if (!b_correctly_extended)
        std::cout << "Not correctly extended." << std::endl;
    if (b_is_a_path_at_all && b_feasible && b_correctly_extended)
    {
        std::cout << "Actual final resource levels:" << std::endl;
        std::cout << "Length: " << actual_final_resource_levels.cost
                  << std::endl;
        std::cout << "Time: " << actual_final_resource_levels.time << std::endl;
        std::cout << "OK." << std::endl;
    }

    return 0;
}
