#include <bits/stdc++.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
using namespace std;


vector<pair<int,int>> corresponding_path(int x, int y, vector<pair<int,int>> s);
int mu_cost(pair<int,int> x, int y);
int c_cost(pair<int,int> x, int y);
int _c(int x, int y);
int _mu(int x, int y);

int labeling_scheme(vector<pair<int,int>> &nodes_graph, vector<vector<int>> &directed_adjency_matrix, 
                vector<int> &topological_sort, vector<pair<int,int>> &initial_cost_function, 
                vector<vector<int>> &initial_cost_bound, vector<pair<int,int>> &initial_fuel_function, 
                vector<vector<int>> &initial_fuel_bound, int n)
{
    int s,t;

    vector<vector<pair<int,int>>> L(n+1);

    for(int i=1;i<=n;i++)
    {
        if(i!=t)
        {
            L[i].clear();
        }
    }

    L[t].push_back(make_pair(0,0));

    int Cout = INT_MAX;
    int Pout = 0;

    vector<pair<int,int>> P;

    for(int j=(n-1); j>=0; j--)
    {
        for(auto node:L[topological_sort[j]])
        {
            // M_j is the fuel label of l_j
            int M_j = node.first;

            // c_j is the cost label of l_j
            int c_j = node.second;
            
            int u_j, u;

            vector<pair<int,int>> P_j = corresponding_path(u_j, u, L[j]);

            for(auto vertex:directed_adjency_matrix[u_j])
            {
                int u_i = vertex;
                
                pair<int,int> a=make_pair(u_i, u_j);
                int M_i = mu_cost(a, M_j);
                int c_i = c_j + c_cost(a, M_j);

                if((c_i + _c(s,u_i)) >= Cout)
                {
                    continue;
                }

                int Mmax;
                
                if((M_i + _mu(s,u_i)) >= Mmax)
                {
                    continue;
                }                
                else
                {
                    P = corresponding_path( _c(s, u_i), a.first, P_j);
                    
                    int c_P;

                    if(P.size()>0 && c_P < Cout)
                    {
                        Cout = c_P;
                        Pout = P.size();
                    }
                }

                bool node_present = false;

                for(auto node:L[u_i])
                {
                    if(node.first == M_i && node.second == c_i)
                    {
                        node_present = true;
                    }
                }

                if(node_present)
                {
                    L[u_i].pop_back();
                    L[u_i].push_back(make_pair(M_i, c_i));
                }
            }
        }
    }

    for(auto value:L[s])
    {
        int M = value.first;
        int c = value.second;

        P = corresponding_path(s,t, L[s]);
        
        int c_P;
        
        if(c_P < Cout)
        {
            Cout = c_P;
            Pout = P.size();
        }
    }

    return Pout;
}


vector<pair<int,int>> corresponding_path(int x, int y, vector<pair<int,int>> &s)
{
    vector<pair<int,int>> var;
    
    for(auto node:s)
    {
        if(node.first==x && node.second==y)
        {
            continue;
        }
        else
        {
            var.push_back(node);
        }
    }
    
    return var;
}


int mu_cost(pair<int,int> x, int y)
{
    if(x.second == y)
    {
        return 0;
    }
    else
    {
        return mu_cost(make_pair(x.first+1, x.second), y);
    }
}

int c_cost(pair<int,int> x, int y)
{
    if(x.second == y)
    {
        return 0;
    }
    else
    {
        return c_cost(make_pair(x.first+1, x.second), y);
    }
}

int pricing_loop(vector<pair<int,int>> &nodes_graph, vector<vector<int>> &directed_adjency_matrix, 
        vector<vector<int>> &compound_variables, vector<int> &topological_sort, vector<pair<int,int>> &initial_cost_function, 
                vector<vector<int>> &initial_cost_bound, vector<pair<int,int>> &initial_fuel_function, 
                vector<vector<int>> &initial_fuel_bound, vector<pair<int,int>> &coordinates, int n)
{
    while(1>0)
    {
        int LP;
        vector<int> x(n+1);
        vector<int> y(n+1);
        vector<int> lambda(n+1);
        vector<int> delta(n+1);

        int Pout = labeling_scheme(nodes_graph, directed_adjency_matrix, topological_sort, 
            initial_cost_function, initial_cost_bound, initial_fuel_function, initial_fuel_bound, n);
        
        int c_Pout;
        if(c_Pout<0)
        {
            nodes_graph[LP].second += Pout;
            continue;
        }

        for(int i=0;i<n;i++)
        {
            int u=nodes_graph[i].first;
            int v=nodes_graph[i].second;
            int a=u;

            int A_fract;
            for(auto node:directed_adjency_matrix[u])
            {
                int a=node;
                if(y[a]!=0 && y[a]!=1)
                {
                    A_fract = node;
                    break;
                }
            }

            int A_comp;
            for(auto node:directed_adjency_matrix[u])
            {
                bool check = false;

                for(auto pre:compound_variables[LP])
                {
                    if(pre == y[a])
                    {
                        check = true;
                        break;
                    }
                }

                if(check)
                {
                    A_comp = node;
                    break;
                }
            }

            int index;
            bool present = false;
            for(auto node:directed_adjency_matrix[A_fract])
            {
                bool check = true;
                for(auto node1:directed_adjency_matrix[A_comp])
                {
                    for(auto node2:directed_adjency_matrix[A_fract])
                    {
                        if(node1 == node2)
                        {
                            check=false;
                            break;
                        }
                    }
                    if(!check)
                    {
                        break;
                    }
                }

                if(check)
                {
                    bool check1=false;
                    for(auto node1:directed_adjency_matrix[A_comp])
                    {
                        if(node1==node)
                        {
                            check1=true;
                            break;
                        }
                    }

                    bool check2=false;
                    for(auto node1:directed_adjency_matrix[A_fract])
                    {
                        if(node1==node)
                        {
                            check2=true;
                            break;
                        }
                    }

                    if(check1 && check2)
                    {
                        index = node;
                        break;
                    }
                }
            }

            if(present)
            {
                compound_variables[LP].push_back(y[index]);

                continue;
            }
        }
        return LP;
        break;
    }
}

int _c(int x, int y)
{
    if(x == 0)
    {
        return 0;
    }
    else
    {
        return _c(x-1, y);
    }
}




int _mu(int x, int y)
{
    if(x == 0)
    {
        return 0;
    }
    else
    {
        return _c(x-1, y);
    }
}


vector<pair<int,int>> corresponding_path(int x, int y, vector<pair<int,int>> s)
{
    vector<pair<int,int>> result;
    for(auto node:s)
    {
        if(node.first!=x && node.second!=y)
        {
            result.push_back(node);
        }
    }
    
    return result;
}



vector<int> topological_sort(vector<pair<int,int>> &nodes_graph, vector<vector<int>> &ls, int k)
{
    vector<int> arr(k);
    stack<int> s;
    
    set<int> st;
    int ind = k - 1;
    
    for (int i = k - 1; i >= 0; i--) 
    {
        if (st.find(i) == st.end()) 
        {
            s.push(i);
            st.insert(i);

            //check all the non visited nodes
            while (!s.empty()) 
            {
                int p = s.top();
                list<int>::iterator it;
                int temp = 0;

                //check its adjacent non visited nodes
                for (auto it = ls[p].begin(); it != ls[p].end(); it++) 
                {
                    if (st.find(*it) == st.end()) 
                    {
                        st.insert(*it);
                        s.push(*it);
                        temp = 1;
                    }
                }

                //if all adjaceny nodes are visited then pop that element from stack
                if (temp == 0) 
                {
                    arr[ind] = p;
                    ind--;
                    s.pop();
                }
            }
        }
    }
    return arr;
}



struct path
{
    int cost;
    vector<int> covered;
};


int greedy_approach_to_solve_set_covering(path feasible_paths[], set<int> &not_included, int n)
{
    int req = 10;
    int is_covered[10] = {0};

    int req_covered = 0;
    set<int> included_paths;

    while (req_covered != req)
    {
        int to_add;
        int found = 0;
        double value = 0;

        for (int p : not_included)
        {
            path fp = feasible_paths[p];
            int temp = 0;

            for (auto x : fp.covered)
            {
                if (is_covered[x] == 0)
                {
                    temp++;
                }
            }

            if (temp != 0 && value < ((double)temp) / fp.cost)
            {
                to_add = p;
                found = 1;
            }
        }
        
        if (found == 1)
        {
            path include = feasible_paths[to_add];
            not_included.erase(to_add);
            included_paths.insert(to_add);

            for (auto x : include.covered)
            {
                if (is_covered[x] == 0)
                {
                    req_covered++;
                }
                is_covered[x] = 1;
            }
        }
        else
        {
            break;
        }
    }

    int total_cost = 0;
    for (int p : included_paths)
    {
        total_cost += feasible_paths[p].cost;
    }
    
    return total_cost;
    
}



// Utility function for
// converting degrees to radians
long double toRadians(const long double degree)
{
    // cmath library in C++
    // defines the constant
    // M_PI as the value of
    // pi accurate to 1e-30
    long double one_deg = (M_PI) / 180;
    return (one_deg * degree);
}

long double distance(long double lat1, long double long1,
                    long double lat2, long double long2)
{
    // Convert the latitudes
    // and longitudes
    // from degree to radians.
    lat1 = toRadians(lat1);
    long1 = toRadians(long1);
    lat2 = toRadians(lat2);
    long2 = toRadians(long2);
    
    // Haversine Formula
    long double dlong = long2 - long1;
    long double dlat = lat2 - lat1;

    long double ans = pow(sin(dlat / 2), 2) +
                        cos(lat1) * cos(lat2) *
                        pow(sin(dlong / 2), 2);

    ans = 2 * asin(sqrt(ans));

    // Radius of Earth in
    // Kilometers, R = 6371
    // Use R = 3956 for miles
    long double R = 6371;
    
    // Calculate the result
    ans = ans * R;

    return ans;
}





int main()
{
    // Precomputation Section

    // no of bases in our model
    int n;  
    cin>>n;

    vector<pair<int,int>> nodes_graph(n+1);

    // Inputing the nodes of all the bases
    for(int i=0;i<n;i++)
    {
        int x,y;
        cin>>x>>y;
        nodes_graph[i].first = x;
        nodes_graph[i].second = y;
    }
    
    vector<pair<int,int>> coordinates(n+1);
    
    // Inputing the coordinates of the bases
    for(int i=1;i<=n;i++)
    {
        int x,y;
        cin>>x>>y;
        coordinates[i].first = x;
        coordinates[i].second = y;
    }
    
    // no of nodes in the graph
    int m;
    cin>>m;

    vector<vector<int>> directed_adjency_matrix(n+1);

    // Inputing the nodes rechechable out of a particular bases in the graph
    for(int i=1;i<=m;i++)
    {
        int x,y;
        cin>>x>>y;
        directed_adjency_matrix[x].push_back(y);
    }

    // Finding the topological ordering of the given adjency matrix of the graph
    vector<int> topological_graph_ordering;
    topological_graph_ordering = topological_sort(nodes_graph, directed_adjency_matrix, n);


    vector<vector<int>> initial_fuel_bound(n+1,vector<int>(n+1));
    vector<vector<int>> initial_cost_bound(n+1,vector<int>(n+1));
    vector<pair<int,int>> initial_fuel_function(n+1);
    vector<pair<int,int>> initial_cost_function(n+1);
    


    bool feasible_path_input;
    cin>>feasible_path_input;
    // Input 0

    path feasible_paths[n];
    set<int> not_included;
    
    if(feasible_path_input == 1)
    {
        for (int i = 0; i < n; i++)
        {
            int cost, size;
            vector<int> arcs;
            cin >> cost >> size;

            for (int j = 0; j < size; j++)
            {
                int a;
                cin >> a;
                arcs.push_back(a);
            }

            path f;
            f.cost = cost;
            f.covered = arcs;
            feasible_paths[i] = f;
            not_included.insert(i);
        }
    }


    // total no of requests to work upon
    int no_of_requests;   
    cin>>no_of_requests;

    while(no_of_requests>0)
    {
        vector<vector<int>> compound_variables(n+1);

        for(int i=1;i<=n;i++)
        {
            int size;
            cin>>size;

            for(int j=0;j<size;j++)
            {
                int x;
                cin>>x;
                compound_variables[i].push_back(x);
            }
        }

        long double start_longitudes, start_latitudes;
        cin>>start_longitudes>>start_latitudes;

        long double end_longitudes, end_latitudes;
        cin>>end_longitudes>>end_latitudes;

        int fuel_required;
        cin>>fuel_required;

        int start_time;
        cin>>start_time;

        double path_distance = distance(start_latitudes, start_longitudes, end_latitudes, end_longitudes);

        
        int result = pricing_loop(nodes_graph, directed_adjency_matrix, compound_variables, 
                    topological_graph_ordering, initial_cost_function, initial_cost_bound, 
                    initial_fuel_function, initial_fuel_bound, coordinates, n);

        cout<<result<<"\n";

        no_of_requests--;
    }

    return 0;
}



















































