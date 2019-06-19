#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
using namespace std;

typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS,
        boost::no_property, boost::property < boost::edge_weight_t, double > > graph_t;
typedef boost::graph_traits < graph_t >::edge_descriptor edge_descriptor_t;
typedef pair<int,int> edge_t;
typedef pair<double,double> vertex_t;

vector< edge_t > vEdges;
map<int,vertex_t> mVertex;

void Read( const std::string& fname )
{
    ifstream f( fname );
    if( ! f.is_open() )
    {
        cout << "Cannot open " << fname << "\n";
        exit(1);
    }
    std::string l;
    while( getline( f, l ))
    {
        cout << "l " << l << "\n";

        vector< std::string > output;
        std::stringstream sst(l);
        std::string a;
        while( getline( sst, a, ' ' ) )
            output.push_back(a);

        if( output[0][0] == 'v' )
        {
            vertex_t v( stof( output[2]), stof( output[3]));
            mVertex.insert( pair<int,vertex_t>( stol( output[1]),v ));
        }
        else if ( output[0][0] == 'e' )
        {
            vEdges.push_back( edge_t( stol( output[1]), stol(output[2]) ) );
        }
    }

//    for( auto e : vEdges )
//        cout << e.first <<" "<< e.second << "\n";
//    cout << "edges\n";
//    for( auto v : mVertex )
//        cout << v.second.first <<" "<< v.second.second << "\n";
//    cout << "\n";
}


//void exMST()
//{
//    using namespace boost;
//    typedef adjacency_list < vecS, vecS, undirectedS,
//            no_property, property < edge_weight_t, int > > Graph;
//    typedef graph_traits < Graph >::edge_descriptor Edge;
//    typedef std::pair<int, int> E;
//
//    const int num_nodes = 5;
//    E edge_array[] = { E(0, 2), E(1, 3), E(1, 4), E(2, 1), E(2, 3),
//                       E(3, 4), E(4, 0), E(4, 1)
//                     };
//    int weights[] = { 1, 1, 2, 7, 3, 1, 1, 1 };
//    std::size_t num_edges = sizeof(edge_array) / sizeof(E);
//#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
//    Graph g(num_nodes);
//    property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
//    for (std::size_t j = 0; j < num_edges; ++j)
//    {
//        Edge e;
//        bool inserted;
//        boost::tie(e, inserted) = add_edge(edge_array[j].first, edge_array[j].second, g);
//        weightmap[e] = weights[j];
//    }
//#else
//    Graph g(edge_array, edge_array + num_edges, weights, num_nodes);
//#endif
//    property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
//    std::vector < Edge > spanning_tree;
//
//    kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));
//
//    std::cout << "Print the edges in the MST:" << std::endl;
//    for (std::vector < Edge >::iterator ei = spanning_tree.begin();
//            ei != spanning_tree.end(); ++ei)
//    {
//        std::cout << source(*ei, g) << " <--> " << target(*ei, g)
//                  << " with weight of " << weight[*ei]
//                  << std::endl;
//    }
//
//    std::ofstream fout("figs/kruskal-eg.dot");
//    fout << "graph A {\n"
//         << " rankdir=LR\n"
//         << " size=\"3,3\"\n"
//         << " ratio=\"filled\"\n"
//         << " edge[style=\"bold\"]\n" << " node[shape=\"circle\"]\n";
//    graph_traits<Graph>::edge_iterator eiter, eiter_end;
//    for (boost::tie(eiter, eiter_end) = edges(g); eiter != eiter_end; ++eiter)
//    {
//        fout << source(*eiter, g) << " -- " << target(*eiter, g);
//        if (std::find(spanning_tree.begin(), spanning_tree.end(), *eiter)
//                != spanning_tree.end())
//            fout << "[color=\"black\", label=\"" << get(edge_weight, g, *eiter)
//                 << "\"];\n";
//        else
//            fout << "[color=\"gray\", label=\"" << get(edge_weight, g, *eiter)
//                 << "\"];\n";
//    }
//    fout << "}\n";
//}

vector < edge_descriptor_t >
testMST( double& length,
         int num_nodes,
         vector<pair<int,int>>& vedge,
         vector<double>& vweight )
{
    using namespace boost;

    vector < edge_descriptor_t > spanning_tree;
    graph_t g(num_nodes);
    property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
    for (std::size_t j = 0; j < vedge.size(); ++j)
    {
        edge_descriptor_t e;
        bool inserted;
        boost::tie(e, inserted) = add_edge(vedge[j].first, vedge[j].second, g);
        weightmap[e] = vweight[j];
    }

    property_map < graph_t, edge_weight_t >::type weight = get(edge_weight, g);

    kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

//    std::cout << "Print the edges in the MST:" << std::endl;
//    for (std::vector < edge_descriptor_t >::iterator ei = spanning_tree.begin();
//            ei != spanning_tree.end(); ++ei)
//    {
//        std::cout << source(*ei, g) << " <--> " << target(*ei, g)
//                  << " with weight of " << weight[*ei]
//                  << std::endl;
//    }

    length = 0;
    for (std::vector < edge_descriptor_t >::iterator ei = spanning_tree.begin();
            ei != spanning_tree.end(); ++ei)
    {
        length += weight[*ei];
    }

    return spanning_tree;
}

void kMST( graph_t& G, int k )
{
    double smallest_length =  std::numeric_limits<double>::max();
    vector < edge_descriptor_t > smallest_mst;

    // loop over distinct pairs of points
    for( int si = 1; si <= mVertex.size(); si++ )
    {
        for( int sj = si+1; sj <= mVertex.size(); sj++ )
        {
            //cout << "s " << si <<" " << sj << "\n";
            /*
            Construct the circle C with diameter  =
            sqrt(3) *   d(i; j) centered at the midpoint of
            the line segment hsi ; sj i.
            */
            auto vi   = mVertex.find( si )->second;
            auto vj   = mVertex.find( sj )->second;
            double dx = vi.first - vj.first;
            double dy = vi.second - vj.second;
            double dC = 1.732 * sqrt( dx*dx+dy*dy);
            double rC = 0.5 * dC;
            double mx = ( vi.first + vj.first ) / 2;
            double my = ( vi.second + vj.second ) / 2;

            // Let SC be the subset of S contained in C .
            vector<int> vSC;
            vSC.push_back( si );
            vSC.push_back( sj );
            for( int i = 1; i <= mVertex.size(); i++ )
            {
                if( i == si )
                    continue;
                if( i == sj )
                    continue;
                auto t = mVertex.find( i )->second;
                double tx = mx - t.first;
                double ty = my - t.second;
                double td  = sqrt( tx*tx+ty*ty );
                if( td < dC )
                    vSC.push_back( i );
            }
            cout << vSC.size() << " Points in SC centered on "
                 << mx << " " << my << " d = "<< dC << "\n";

            /*
            if SC contains fewer than k points, skip
            to the next iteration of the loop (i.e., try the next pair of points).
            */
            if( vSC.size() < k )
                continue;


            for( auto p : vSC )
                cout << p << " ";
            cout << "\n";

            /*
            Let Q b e the square of side  circumscribing C .
            (4) Divide Q in to k square cells each with side = d / sqrt( K )
            */
            vector< vector< int > > v_pts_in_cell;
            double cellside = dC/sqrt( k );
            double blx = mx - rC;
            double bly = my - rC;

            //cout << " Q "<<blx<<" "<<bly<<" "<<blx+dC<<" "<<bly+dC<<" cellside " << cellside << "\n";
            for( ; ; )
            {
                vector< int > v;
                for( int i = 1; i <= mVertex.size(); i++ )
                {
                    auto t = mVertex.find( i )->second;
                    double tx = t.first;
                    double ty = t.second;
                    if( blx <= tx && tx <= blx+cellside
                            &&
                            bly <= ty && ty <= bly+cellside )
                    {
                        v.push_back( i );
                    }
                }
                v_pts_in_cell.push_back( v );

                //cout << "cell " << blx <<" "<< bly << " has " << v.size() << "\n";

                blx += cellside;
                if( blx >= mx + rC )
                {
                    blx = mx - rC;
                    bly += cellside;
                    if( bly >= my + rC )
                    {
                        break;
                    }
                }
            }
            // Sort the cells by the number of points from SC they contain
            sort(v_pts_in_cell.begin(),
                 v_pts_in_cell.end(),
                 []( vector<int> v1, vector<int> v2 )
            {
                return ( v1.size() > v2.size() );
            });

            /** choose the minimum number of cells so that the chosen cells together contain at least k
            points. If necessary , arbitrarily discard points from the last chosen cell so that
            the total number of points in all the cells is equal to k .
            */
            vector<int> v_pts_in_tree;
            int count = 0;
            for( auto& v : v_pts_in_cell )
            {
                for( int p : v )
                {
                    count++;
                    if( count <= k )
                        v_pts_in_tree.push_back( p );
                }
                if( count > k )
                    break;
            }
            cout << "points in test minimum spanning tree ";
            for( auto p : v_pts_in_tree )
                cout << p << " ";
            cout << "\n";

            vector< edge_t > v_edges_in_tree;
            vector<double> vdist;
            auto es = boost::edges(G);
            for (auto eit = es.first; eit != es.second; ++eit)
            {
                if( ( find( v_pts_in_tree.begin(), v_pts_in_tree.end(), boost::source(*eit, G) ) != v_pts_in_tree.end()  ) &&
                        ( find( v_pts_in_tree.begin(), v_pts_in_tree.end(), boost::target(*eit, G) ) != v_pts_in_tree.end()  ) )
                {
                    cout << "v"<<boost::source(*eit, G) << " - v" << boost::target(*eit, G) << "\n";
                    v_edges_in_tree.push_back( edge_t(boost::source(*eit, G), boost::target(*eit, G) ) );
                    auto v = mVertex.find( boost::source(*eit, G))->second;
                    //cout << " ( " << v.first <<" " << v.second << ") - ";
                    auto v2 = mVertex.find( boost::target(*eit, G))->second;
                    //cout << " ( " << v2.first <<" " << v2.second << ") ";
                    double dx = v.first - v2.first;
                    double dy = v.second - v2.second;
                    double d = sqrt( dx*dx+dy*dy);
                    vdist.push_back( d );
                }
            }

            double length;
            vector < edge_descriptor_t > mst = testMST(
                                                   length,
                                                   v_pts_in_tree.size(),
                                                   v_edges_in_tree,
                                                   vdist );
            if( length < smallest_length )
            {
                smallest_length = length;
                smallest_mst = mst;
            }

        }
    }
    cout << "\n\nSmallest " << k << " point MST has length " << smallest_length << "\n";
    for ( auto ei = smallest_mst.begin();
            ei != smallest_mst.end(); ++ei)
    {
        auto v = mVertex.find( ei->m_source)->second;
        auto v2 = mVertex.find( ei->m_target)->second;
        double dx = v.first - v2.first;
        double dy = v.second - v2.second;
        double d = sqrt( dx*dx+dy*dy);

        std::cout << ei->m_source << " <--> " << ei->m_target
                  << " with length of " << d
                  << std::endl;
    }
}

int main()
{
    Read( "s.txt" );


    graph_t G( vEdges.begin(), vEdges.end(), mVertex.size() );

    cout << "edges\n";
    auto es = boost::edges(G);
    for (auto eit = es.first; eit != es.second; ++eit)
    {
        cout << "v"<<boost::source(*eit, G) << " - v" << boost::target(*eit, G);
        auto v = mVertex.find( boost::source(*eit, G))->second;
        //cout << " ( " << v.first <<" " << v.second << ") - ";
        auto v2 = mVertex.find( boost::target(*eit, G))->second;
        //cout << " ( " << v2.first <<" " << v2.second << ") ";
        double dx = v.first - v2.first;
        double dy = v.second - v2.second;
        double d = sqrt( dx*dx+dy*dy);
        cout <<" "<< d << "\n";
    }



    //exMST();

//    typedef std::pair<int, int> E;
//    E edge_array[] = { E(0, 2), E(1, 3), E(1, 4), E(2, 1), E(2, 3),
//                       E(3, 4), E(4, 0), E(4, 1)
//                     };
//    vector<pair<int,int>> vedge;
//    for( int k = 0; k < 8; k++ )
//        vedge.push_back( edge_array[k] );
//    vector < double > weights { 1, 1, 2, 7, 3, 1, 1, 1 };
//    testMST( 5, vedge, weights );

    kMST( G, 3 );

    return 0;
}
