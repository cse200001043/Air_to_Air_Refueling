#include <bits/stdc++.h>
using namespace std;

struct path
{
    int cost;
    vector<int> covered;
};

int main()
{
    int n;
    int req = 10;
    // cout << "enter Number of feasible paths: ";
    cin >> n;
    int is_covered[10] = {0};
    path feasible_paths[n];
    set<int> not_included;

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
        cout << p << " ";
        total_cost += feasible_paths[p].cost;
    }
    cout<<"\n";

    cout << "total cost is: " << total_cost << "\n";

    return 0;
}