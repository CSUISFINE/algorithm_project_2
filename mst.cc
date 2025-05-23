#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>

using namespace std;

struct City {
    int id;
    double x, y;
};

double euc_2D(const City& a, const City& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return round(sqrt(dx * dx + dy * dy));
}

vector<City> load_tsp(const string& filename) {
    ifstream infile(filename);
    string line;
    vector<City> cities;

    while (getline(infile, line)) {
        if (line.find("NODE_COORD_SECTION") != string::npos) break;
    }
    while (getline(infile, line)) {
        if (line == "EOF") break;
        istringstream iss(line);
        int id;
        double x, y;
        if (!(iss >> id >> x >> y)) continue;
        cities.push_back({id, x, y});
    }
    
    return cities;
}

double prim_mst(const vector<City>& cities, vector<vector<pair<int, double>>>& adj_list) {
    int n = cities.size();
    vector<bool> in_mst(n, false);
    vector<double> min_dist(n, numeric_limits<double>::max());
    vector<int> parent(n, -1);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

    min_dist[0] = 0.0;
    pq.push({0.0, 0});

    double mst_weight = 0.0;
    int edges = 0;
    
    adj_list.assign(n, vector<pair<int, double>>());

    while (!pq.empty()) {
        double d = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        if (in_mst[u])
            continue;

        in_mst[u] = true;
        mst_weight += d;
        edges++;

        if (edges % 1000 == 0) { 
            cout << "The " << edges << "th node (index: " << u << ", ID: " << cities[u].id 
             << ") was added to MST. Total nodes currently in MST: " << edges << endl;
        }

        if (parent[u] != -1) {
            adj_list[u].push_back({parent[u], d});
            adj_list[parent[u]].push_back({u, d});
        }

        for (int i = 0; i < n; i++) {
            if (u == i)
                continue;
            double w = euc_2D(cities[u], cities[i]);
            if (!in_mst[i] && w < min_dist[i]) {
                min_dist[i] = w;
                pq.push({min_dist[i], i});
                parent[i] = u;
            }
        }
    }

    return mst_weight;
}

pair<vector<pair<int, int>>, double> greedy_perfect_matching(const vector<int>& odd_d_idx, const vector<City>& cities) {
    std::vector<std::pair<int, int>> matching_edges;
    double pm_weight = 0.0;
    int n = odd_d_idx.size();
    vector<bool> matched(n, false);
    const double INF = std::numeric_limits<double>::max();

    for (int i = 0; i < n; i++) {
        if (matched[i])
            continue;
        int u_idx = odd_d_idx[i];
        int v_idx_idx = -1;
        double min_dist = INF;
        for (int j = i + 1; j < n; j++) {
            if (matched[j])
                continue;
            int temp_v = odd_d_idx[j];
            double cur_dist = euc_2D(cities[u_idx], cities[temp_v]);
            if (cur_dist < min_dist) {
                min_dist = cur_dist;
                v_idx_idx = j;
            }
        }
        if (v_idx_idx != -1) {
            int v_idx = odd_d_idx[v_idx_idx];
            matching_edges.push_back({u_idx, v_idx});
            pm_weight += min_dist;
            matched[i] = true;
            matched[v_idx_idx] = true;
        }
    }
    return {matching_edges, pm_weight};
}

void two_approximation(vector<vector<pair<int, double>>>& adj_list_2_aprx) {
    for (auto& row_vec : adj_list_2_aprx) {
        row_vec.insert(row_vec.end(), row_vec.begin(), row_vec.end());
    }
}

vector<int> euler_circuit(const vector<vector<pair<int, double>>>& adj_list) {
    int cities_num = adj_list.size();
    vector<int> circuit;
    vector<map<int, int>> adj_counts(cities_num);

    for (int u = 0; u < cities_num; u++) {
        for (auto &edge : adj_list[u]) {
            int v = edge.first;
            adj_counts[u][v]++;
        }
    }

    int start = 0;
    vector<int> current_path_stack;
    current_path_stack.push_back(start);

    while(!current_path_stack.empty()) {
        int u = current_path_stack.back();

        if (adj_counts[u].empty()) {
            circuit.push_back(u);
            current_path_stack.pop_back();
        } else {
            int v = adj_counts[u].begin()->first;
            adj_counts[u][v]--;
            if (adj_counts[u][v] == 0) {
                adj_counts[u].erase(v);
            }
            adj_counts[v][u]--;
            if (adj_counts[v][u] == 0) {
                adj_counts[v].erase(u);
            }
            current_path_stack.push_back(v);
        }
    }
    reverse(circuit.begin(), circuit.end());
    return circuit;
}

vector<int> hamilton_circuit(const vector<int>& euler_circuit, int num_cities) {
    vector<int> tsp;
    vector<bool> visited(num_cities, false);
    for (int city_idx : euler_circuit) {
        if (!visited[city_idx]) {
            tsp.push_back(city_idx);
            visited[city_idx] = true;
        }
    }
    if (tsp.front() != tsp.back() || tsp.size() == 1) {
        tsp.push_back(tsp.front());
    }
    return tsp;
}

double calculate_cost(const vector<int>& tour, const vector<City>& cities) {
    double cost = 0.0;
    for (int i = 0; i < tour.size() - 1; i++) {
        cost += euc_2D(cities[tour[i]], cities[tour[i + 1]]);
    }
    return cost;
}

int main() {
    string filename = "mona-lisa100K.tsp";
    vector<City> cities = load_tsp(filename);

    if (cities.empty()) {
        cerr << "Failed to load cities or file is empty." << endl;
        return 1;
    }

    cout << "Loaded " << cities.size() << " cities." << endl;

    vector<vector<pair<int, double>>> adj_list_mst;
    double mst_w = prim_mst(cities, adj_list_mst);

    cout << "MST Weight: " << fixed << setprecision(2) << mst_w << endl;

    vector<int> odd_degree_vertices_indices;
    for(int i=0; i < cities.size(); ++i) {
        if(adj_list_mst[i].size() % 2 != 0) {
            odd_degree_vertices_indices.push_back(i);
        }
    }

    cout << "Odd degree vertices (0-based city index): ";
    for(int city_idx : odd_degree_vertices_indices) {
        cout << city_idx << "(ID:" << cities[city_idx].id << ") ";
    }
    cout << endl;
    cout << "Number of odd degree vertices: " << odd_degree_vertices_indices.size() << endl;

    // perfect matching by greedy
    pair<vector<pair<int, int>>, double> matching_result = greedy_perfect_matching(odd_degree_vertices_indices, cities);

    vector<pair<int, int>> matched_pairs = matching_result.first;
    double total_weight = matching_result.second;

    cout << "\n greedy perfect matching result:" << endl;
    for (const auto& pair : matched_pairs) {
        cout << "  matched edge: (" << pair.first << ", " << pair.second << ") "
                  << "  weight: " << euc_2D(cities[pair.first], cities[pair.second]) << endl;
    }
    cout << "matching edges number: " << matched_pairs.size() << endl;
    cout << "total matching cost: " << fixed << setprecision(2) << total_weight << endl;
    
    vector<vector<pair<int, double>>> adj_list_pm;
    adj_list_pm = adj_list_mst;
    for (const auto& pm_edge : matched_pairs) {
        int u = pm_edge.first;
        int v = pm_edge.second;
        double weight = euc_2D(cities[u], cities[v]);
        adj_list_pm[u].push_back({v, weight});
        adj_list_pm[v].push_back({u, weight});
    }
    // 2-approximation algorithm
    vector<vector<pair<int, double>>> adj_list_2_aprx;
    adj_list_2_aprx = adj_list_mst;
    two_approximation(adj_list_2_aprx);

    vector<int> euler_circ_pm = euler_circuit(adj_list_pm); 
    vector<int> euler_circ_2_aprx = euler_circuit(adj_list_2_aprx); 

    // make hamilton circuit using by shortcutting
    vector<int> hamil_circ_pm = hamilton_circuit(euler_circ_pm, cities.size());
    vector<int> hamil_circ_2_aprx = hamilton_circuit(euler_circ_2_aprx, cities.size());
    cout << "\nHamiltonian Circuit by greedy perfect matching (TSP Tour - city ID):" << endl;
    if (hamil_circ_pm.empty()) {
        cout << "  No Hamiltonian circuit generated." << endl;
    } else {
        cout << "  ";
        for (size_t i = 0; i < hamil_circ_pm.size(); i++) {
            cout << cities[hamil_circ_pm[i]].id << (i == hamil_circ_pm.size() - 1 ? "" : " -> ");
        }
        cout << endl;
    }
    cout << "\nHamiltonian Circuit by 2-approximation (TSP Tour - city ID):" << endl;
    if (hamil_circ_2_aprx.empty()) {
        cout << "  No Hamiltonian circuit generated." << endl;
    } else {
        cout << "  ";
        for (size_t i = 0; i < hamil_circ_2_aprx.size(); i++) {
            cout << cities[hamil_circ_2_aprx[i]].id << (i == hamil_circ_2_aprx.size() - 1 ? "" : " -> ");
        }
        cout << endl;
    }

    // total TSP cost
    double final_cost_pm = calculate_cost(hamil_circ_pm, cities);
    double final_cost_2_aprx = calculate_cost(hamil_circ_2_aprx, cities);
    double optimal_tour_cost = 2579.0;
    cout << "\nTotal TSP tour cost (Christofides with Greedy Matching): " << fixed << setprecision(2) << final_cost_pm << endl;
    cout << "Total TSP tour cost (MST-based 2-approximation algorithm): " << fixed << setprecision(2) << final_cost_2_aprx << endl;
    cout << "\nOptimal TSP tour cost: " << fixed << setprecision(2) << optimal_tour_cost << ", ratio: " << fixed << setprecision(3) << optimal_tour_cost / mst_w << endl;
    // Approximation ratio
    cout << "Approximation Ratio (Greedy Matching Tour Cost / optimal value): " << fixed << setprecision(3) << final_cost_pm / optimal_tour_cost << endl;
    cout << "Approximation Ratio (2-approximation Tour Cost / optimal value): " << fixed << setprecision(3) << final_cost_2_aprx / optimal_tour_cost << endl;

    return 0;
}
