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

class cEdge
{
public:
    int v1;
    int v2;
    double length;

    cEdge( int i1, int i2 )
        : v1( i1 )
        , v2( i2 )
    {

    }
    void calc_length();
};

vector< cEdge > vEdges;
map<int,vertex_t> mVertex;

void cEdge::calc_length()
{
    auto p1 = mVertex.find( v1 )->second;
    auto p2 = mVertex.find( v2 )->second;
    double dx = p1.first -p2.first;
    double dy = p1.second - p2.second;
    length = sqrt( dx*dx+dy*dy);
}

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
            vEdges.push_back( cEdge( stol( output[1]), stol(output[2]) ) );
        }
    }

    for( auto& e : vEdges )
        e.calc_length();

//    for( auto e : vEdges )
//        cout << e.first <<" "<< e.second << "\n";
//    cout << "edges\n";
//    for( auto v : mVertex )
//        cout << v.second.first <<" "<< v.second.second << "\n";
//    cout << "\n";
}


void kMST(  int k )
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

            // construct graph including only the points selected
            graph_t g;
            boost::property_map<graph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, g);
            for( auto& e : vEdges )
            {
                if( ( find( v_pts_in_tree.begin(), v_pts_in_tree.end(), e.v1 ) != v_pts_in_tree.end()  ) &&
                        ( find( v_pts_in_tree.begin(), v_pts_in_tree.end(), e.v2 ) != v_pts_in_tree.end()  ) )
                {
                    edge_descriptor_t ed;
                    bool inserted;
                    boost::tie(ed, inserted) = add_edge( e.v1, e.v2, g);
                    weightmap[ed] = e.length;
                }
            }

            // calculate minimum spanning tree
            vector < edge_descriptor_t > mst;
            boost::kruskal_minimum_spanning_tree(g, std::back_inserter(mst));

            // calculate length of edges in spanning tree
            double length = 0;
            for (auto ei = mst.begin(); ei != mst.end(); ++ei)
            {
                length += weightmap[*ei];
            }

            // is this the smallest so far?
            if( length < smallest_length )
            {
                smallest_length = length;
                smallest_mst = mst;
            }

        }
    }

    // Display best result
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

    kMST( 3 );

    return 0;
}
