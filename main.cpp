#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
using namespace std;

typedef boost::adjacency_list<
boost::listS, boost::vecS, boost::bidirectionalS > graph_t;
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

void kMST( graph_t& G, int k )
{
    // loop over distinct pairs of points
    for( int si = 1; si <= mVertex.size(); si++ )
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

                cout << "cell " << blx <<" "<< bly << " has " << v.size() << "\n";

                blx += cellside;
                if( blx >= mx + rC )
                {
                    blx = mx - rC;
                    bly += cellside;
                    if( bly >= my + rC )    {
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
            cout << "points in minimum spanning tree ";
            for( auto p : v_pts_in_tree )
                cout << p << " ";
            cout << "\n";
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

    kMST( G, 3 );

    return 0;
}
