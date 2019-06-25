#include <iostream>
#include <fstream>
#include <map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
using namespace std;


class cVertex
{
public:
    cVertex()
        : x( numeric_limits<double>::lowest())
    {

    }
    double x;
    double y;
};

class cEdge
{
public:
    cEdge()
        : length( -1 )
    {

    }
    double length;
};

class cCircle
{
public:
    double x;
    double y;
    void Diameter( double d )
    {
        myD = d;
        myR = d/2;
    }
    double Diameter() const
    {
        return myD;
    }
    double Radius() const
    {
        return myR;
    }
private:
    double myD;
    double myR;
};

typedef boost::adjacency_list<
boost::listS, boost::vecS, boost::bidirectionalS,
      cVertex, cEdge > graph_t;


graph_t Read2( const std::string& fname )
{
    ifstream f( fname );
    if( ! f.is_open() )
    {
        cout << "Cannot open " << fname << "\n";
        exit(1);
    }
    struct sV
    {
        int index;
        double x;
        double y;
    };
    vector< sV > vVertex;
    vector<pair<int,int>> vEdge;
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
            sV V;
            V.index = stol( output[1]);
            V.x     = stof( output[2]);
            V.y     = stof( output[3]);
            vVertex.push_back( V );
        }
        else if ( output[0][0] == 'e' )
        {
            vEdge.push_back( make_pair( stol( output[1]), stol(output[2]) ) );
        }
    }

    graph_t g;

    // add edges to graph
    for( auto& e : vEdge )
        add_edge( e.first, e.second, g );

    // add properties
    double x1, y1, x2, y2;
    auto es = boost::edges(g);
    for (auto eit = es.first; eit != es.second; ++eit)
    {
        for( auto& v : vVertex )
        {
            int f = 0;
            int v1 =  source(*eit, g);
            int v2 =  target(*eit, g);
            if( v.index == v1 )
            {
                x1 = v.x;
                y1 = v.y;
                g[ v1 ].x = x1;
                g[ v1 ].y = y1;
                f++;
            }
            else if ( v.index == v2 )
            {
                x2 = v.x;
                y2 = v.y;
                g[ v2 ].x = x2;
                g[ v2 ].y = y2;
                f++;
            }
            if( f == 2 )
                break;
        }
        double dx = x1 - x2;
        double dy = y1 - y2;
        g[ *eit ].length = sqrt(dx*dx+dy*dy);
        cout << source(*eit, g) << " "
             << target(*eit, g) << " "
             << g[ *eit ].length << "\n";
    }
    return g;
}

void kMST( int k, graph_t& g )
{
    typedef boost::graph_traits < graph_t >::edge_descriptor edge_descriptor_t;
    double smallest_length = numeric_limits<double>::max();
    vector < edge_descriptor_t > smallest_mst;

    // loop over distinct pairs of points
    auto vs = boost::vertices(g);
    for (auto vit1 = vs.first; vit1 != vs.second; ++vit1 )
    {
        auto vit2 = vit1;
        if( g[*vit1].x < std::numeric_limits<int>::lowest() + 1 )
            continue;
        for( ++vit2; vit2 != vs.second; ++vit2 )
        {
            if( g[*vit2].x < std::numeric_limits<int>::lowest() + 1 )
                continue;
            cout << *vit1 <<" " << *vit2 << "\n";

            /* Construct the circle C with diameter  =
            sqrt(3) *   d(i; j) centered at the midpoint of
            the line segment */
            cCircle C;
            double x1 = g[*vit1].x;
            double y1 = g[*vit1].y;
            double x2 = g[*vit2].x;
            double y2 = g[*vit2].y;
            double dx = x1 - x2;
            double dy = y1 - y2;
            C.Diameter( 1.732 * sqrt( dx*dx+dy*dy) );
            C.x = ( x1 + x2 ) / 2;
            C.y = ( y1 + y2 ) / 2;

            // Let SC be the subset of S contained in C .
            vector<int> vSC;
            vSC.push_back( *vit1 );
            vSC.push_back( *vit2 );
            for( auto vit3 = vs.first; vit3 != vs.second; ++vit3  )
            {
                if( *vit3 == *vit1 )
                    continue;
                if( *vit3 == *vit2 )
                    continue;
                dx = g[*vit3].x - C.x;
                dy = g[*vit3].y - C.y;
                double td  = sqrt( dx*dx+dy*dy );
                if( td < C.Diameter() )
                    vSC.push_back( *vit3 );
            }
            // if SC contains fewer than k points, skip
            if( vSC.size() < k )
                continue;
            cout << "Points in SC centered on " << C.x << " " << C.y << "\n";
            for( int p : vSC )
                cout << p << " ";
            cout << "\n\n";

            /*
            Let Q b e the square of side  circumscribing C .
            (4) Divide Q in to k square cells each with side = d / sqrt( K )
            */
            vector< vector< int > > v_pts_in_cell;
            double cellside = C.Diameter()/sqrt( k );
            double blx = C.x - C.Radius();
            double bly = C.y - C.Radius();

            //cout << " Q "<<blx<<" "<<bly<<" "<<blx+dC<<" "<<bly+dC<<" cellside " << cellside << "\n";
            for( ; ; )
            {
                vector< int > v;
                for( int vt : vSC )
                {
                    double tx = g[vt].x;
                    double ty = g[vt].y;
                    if( blx <= tx && tx <= blx+cellside
                            &&
                            bly <= ty && ty <= bly+cellside )
                    {
                        v.push_back( vt );
                    }
                }
                v_pts_in_cell.push_back( v );

               // cout << "cell " << blx <<" "<< bly << " has " << v.size() << "\n";

                blx += cellside;
                if( blx >= C.x + C.Radius() )
                {
                    blx = C.x - C.Radius();
                    bly += cellside;
                    if( bly >= C.y + C.Radius() )
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

            // construct subgraph including only edges between the points selected
            graph_t sub;

            // loop over all edges
            auto es = boost::edges( g );
            for( auto eit = es.first; eit != es.second; eit++ )
            {
                //cout << source(*eit, g) <<" "<< target(*eit,g) << "\n";

                // check if both target and source are to be included
                if( ( find( v_pts_in_tree.begin(), v_pts_in_tree.end(), source(*eit, g) ) != v_pts_in_tree.end()  ) &&
                        ( find( v_pts_in_tree.begin(), v_pts_in_tree.end(), target(*eit,g) ) != v_pts_in_tree.end()  ) )
                {
                    // copy edge from main graph to subgraph
                    edge_descriptor_t ed;
                    bool inserted;
                    boost::tie(ed, inserted) = add_edge( source(*eit, g), target(*eit,g), sub);
                    sub[ed].length = g[*eit].length;
                }
            }

            // calculate minimum spanning tree
            vector <  boost::graph_traits < graph_t >::edge_descriptor > mst;
            boost::kruskal_minimum_spanning_tree(
                sub,
                std::back_inserter(mst),
                boost::weight_map(get(&cEdge::length, sub)));

            double mstLength = 0;
            for( auto e : mst )
            {
                std::cout << e.m_source << " <--> " << e.m_target <<" l=" << sub[e].length <<", ";
                mstLength += sub[e].length;
            }
            cout << "\n";

            // is this the smallest so far?
            if( mstLength < smallest_length )
            {
                smallest_length = mstLength;
                smallest_mst = mst;
            }
        }
    }

    // Display best result
    cout << "\n\nSmallest " << k << " point MST has length " << smallest_length << "\n";

    for( auto e : smallest_mst )
    {
        std::cout << e.m_source << " <--> " << e.m_target <<", ";
    }
    cout << "\n";
}
int main()
{
    graph_t g = Read2( "s.txt" );

    kMST( 3, g );

    return 0;
}
