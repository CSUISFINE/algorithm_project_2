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
#include <random>
#include <chrono>
#include <thread>

using namespace std;
using namespace std::chrono;

struct City {
    int id;
    double x, y;
};

double euc_2D(const City& a, const City& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
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

        // if (edges % 1000 == 0) { 
        //     cout << "The " << edges << "th node (index: " << u << ", ID: " << cities[u].id 
        //      << ") was added to MST. Total nodes currently in MST: " << edges << endl;
        // }

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

vector<int> euler_circuit(const vector<vector<pair<int, double>>>& adj_list) {
    int num_cities = adj_list.size();
    vector<int> circuit;
    vector<map<int, int>> adj_counts(num_cities);

    for (int u = 0; u < num_cities; u++) {
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

pair<vector<int>, double> double_bridge_move(const vector<int>& tour, const vector<City>& cities, mt19937& gen) {
    int num_cities = tour.size() - 1;

    uniform_int_distribution<> dis(1, num_cities - 1);
    int points[4];
    while(true){
        for (int i = 0; i < 4; i++){
            points[i] = dis(gen);
        }
        sort(points, points + 4);
        if (points[0] == points[1] || points[1] == points[2] || points[2] == points[3]) {
            continue;
        } else {
            break;
        }
    }

    vector<int> new_tour;
    new_tour.insert(new_tour.end(), tour.begin(), tour.begin() + points[0]);
    new_tour.insert(new_tour.end(), tour.begin() + points[2], tour.begin() + points[3]);
    new_tour.insert(new_tour.end(), tour.begin() + points[1], tour.begin() + points[2]);
    new_tour.insert(new_tour.end(), tour.begin() + points[0], tour.begin() + points[1]);
    new_tour.insert(new_tour.end(), tour.begin() + points[3], tour.begin() + num_cities);
    new_tour.push_back(new_tour[0]);

    double cost_added = euc_2D(cities[tour[points[0]-1]], cities[tour[points[2]]]) + 
                        euc_2D(cities[tour[points[3]-1]], cities[tour[points[1]]]) + 
                        euc_2D(cities[tour[points[2]-1]], cities[tour[points[0]]]) + 
                        euc_2D(cities[tour[points[1]-1]], cities[tour[points[3]]]);
    double cost_removed = euc_2D(cities[tour[points[0]-1]], cities[tour[points[0]]]) +
                            euc_2D(cities[tour[points[1]-1]], cities[tour[points[1]]]) +
                            euc_2D(cities[tour[points[2]-1]], cities[tour[points[2]]]) +
                            euc_2D(cities[tour[points[3]-1]], cities[tour[points[3]]]);
    return {new_tour, cost_added - cost_removed};
}

pair<vector<int>, double> three_opt(const vector<int>& tour, const vector<City>& cities, mt19937& gen) {
    int num_cities = tour.size() - 1;
    uniform_int_distribution<> dis_idx(1, num_cities - 1);
    int i, j, k;
    int p[3];
    while(true) {
        for (int i = 0; i < 3; i++) {
            p[i] = dis_idx(gen);
        }
        sort(p, p + 3);
        i = p[0]; j = p[1]; k = p[2];
        if (abs(i - j) < 2 || abs(j - k) < 2 || abs(k - i) < 2){
            continue;
        } else {
            break;
        }
    }
    vector<int> S_opt = tour;
    double cost_removed = euc_2D(cities[tour[i - 1]], cities[tour[i]]) +
                            euc_2D(cities[tour[j - 1]], cities[tour[j]]) +
                            euc_2D(cities[tour[k - 1]], cities[tour[k]]);
    double cost_added = 0;
    double cost_change[8];
    double best_cost_change = numeric_limits<double>::infinity();
    int best_opt_case = 0;

    for (int opt_case = 1; opt_case < 8; opt_case++) {
        switch (opt_case) {
        case 0: // original
            cost_change[0] = 0;
            break;
        case 1: // reverse i ... j - 1
            cost_change[1] = -(euc_2D(cities[tour[i - 1]], cities[tour[i]]) +
                            euc_2D(cities[tour[j - 1]], cities[tour[j]])) + 
                            euc_2D(cities[tour[i - 1]], cities[tour[j - 1]]) + 
                            euc_2D(cities[tour[i]], cities[tour[j]]);
            break;
        case 2: // reverse j ... k - 1
            cost_change[2] = -(euc_2D(cities[tour[j - 1]], cities[tour[j]]) +
                            euc_2D(cities[tour[k - 1]], cities[tour[k]])) + 
                            euc_2D(cities[tour[j - 1]], cities[tour[k - 1]]) + 
                            euc_2D(cities[tour[j]], cities[tour[k]]);
            break;
        case 3: // reverse i ... k - 1
            cost_change[3] = -(euc_2D(cities[tour[i - 1]], cities[tour[i]]) +
                            euc_2D(cities[tour[k - 1]], cities[tour[k]])) + 
                            euc_2D(cities[tour[i - 1]], cities[tour[k - 1]]) + 
                            euc_2D(cities[tour[i]], cities[tour[k]]);
            break;
        case 4: // i - 1, j ... k - 1, i ... j - 1, k
            cost_change[4] = -cost_removed + 
                            euc_2D(cities[tour[i - 1]], cities[tour[j]]) + 
                            euc_2D(cities[tour[k - 1]], cities[tour[i]]) + 
                            euc_2D(cities[tour[j - 1]], cities[tour[k]]);
            break;
        case 5: // reverse i ... j - 1 after that, make i - 1, j ... k - 1, j - 1 ... i, k
            cost_change[5] = -cost_removed + 
                            euc_2D(cities[tour[i - 1]], cities[tour[j]]) + 
                            euc_2D(cities[tour[k - 1]], cities[tour[j - 1]]) + 
                            euc_2D(cities[tour[i]], cities[tour[k]]);
            break;
        case 6: // reverse j ... k - 1 after that, make i - 1, k - 1 ... j, i ... j - 1, k
            cost_change[6] = -cost_removed + 
                            euc_2D(cities[tour[i - 1]], cities[tour[k - 1]]) + 
                            euc_2D(cities[tour[j]], cities[tour[i]]) + 
                            euc_2D(cities[tour[j - 1]], cities[tour[k]]);
            break;
        case 7: // reverse i ... j - 1 and j ... k - 1, make i - 1, j - 1 ... i, k - 1 ... j, k
            cost_change[7] = -cost_removed + 
                            euc_2D(cities[tour[i - 1]], cities[tour[j - 1]]) + 
                            euc_2D(cities[tour[i]], cities[tour[k - 1]]) + 
                            euc_2D(cities[tour[j]], cities[tour[k]]);
            break;
        default:
            return {tour, 0};
        }
        if(best_cost_change > cost_change[opt_case]){
            best_cost_change = cost_change[opt_case];
            best_opt_case = opt_case;
        }
    }

    switch(best_opt_case) {
        case 0:
            break;
        case 1:
            reverse(S_opt.begin() + i, S_opt.begin() + j);
            break;
        case 2:
            reverse(S_opt.begin() + j, S_opt.begin() + k);
            break;
        case 3:
            reverse(S_opt.begin() + i, S_opt.begin() + k);
            break;
        case 4:
            {
            vector<int> new_tour;
            new_tour.insert(new_tour.end(), tour.begin(), tour.begin() + i);
            new_tour.insert(new_tour.end(), tour.begin() + j, tour.begin() + k);
            new_tour.insert(new_tour.end(), tour.begin() + i, tour.begin() + j);
            new_tour.insert(new_tour.end(), tour.begin() + k, tour.end());
            S_opt = new_tour;
            break;
            }
        case 5:
            {
            vector<int> new_tour;
            new_tour.insert(new_tour.end(), tour.begin(), tour.begin() + i);
            new_tour.insert(new_tour.end(), tour.begin() + j, tour.begin() + k);
            vector<int> rev(tour.begin() + i, tour.begin() + j);
            reverse(rev.begin(), rev.end());
            new_tour.insert(new_tour.end(), rev.begin(), rev.end());
            new_tour.insert(new_tour.end(), tour.begin() + k, tour.end());
            S_opt = new_tour;
            break;
            }
        case 6:
            {
            vector<int> new_tour;
            new_tour.insert(new_tour.end(), tour.begin(), tour.begin() + i);
            vector<int> rev(tour.begin() + j, tour.begin() + k);
            reverse(rev.begin(), rev.end());
            new_tour.insert(new_tour.end(), rev.begin(), rev.end());
            new_tour.insert(new_tour.end(), tour.begin() + i, tour.begin() + j);
            new_tour.insert(new_tour.end(), tour.begin() + k, tour.end());
            S_opt = new_tour;
            break;
            }
        case 7:
            {
            vector<int> new_tour;
            new_tour.insert(new_tour.end(), tour.begin(), tour.begin() + i);
            vector<int> rev_B(tour.begin() + i, tour.begin() + j);
            reverse(rev_B.begin(), rev_B.end());
            new_tour.insert(new_tour.end(), rev_B.begin(), rev_B.end());
            vector<int> rev_C(tour.begin() + j, tour.begin() + k);
            reverse(rev_C.begin(), rev_C.end());
            new_tour.insert(new_tour.end(), rev_C.begin(), rev_C.end());
            new_tour.insert(new_tour.end(), tour.begin() + k, tour.end());
            S_opt = new_tour;
            break;
            }
        default:
            break;
    }
    return {S_opt, best_cost_change};
}

struct SAParameters {
    double P_init_accept;
    double T_final;
    double alpha;
    double K;
    double L_multiplier;
    double T_transition_factor;
};

// circuit: 1, 2, 3, ..., n, 1
vector<int> SA_3opt(const vector<int>& initial_circuit, const vector<City>& cities, double avg_delta_E_pos, const SAParameters& params) {
    int num_cities = initial_circuit.size() - 1;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis_prob(0, 1.0);
    uniform_int_distribution<> dis_idx(1, num_cities - 1);
    vector<int> S_cur = initial_circuit;
    vector<int> S_best = S_cur;
    double cost_cur = calculate_cost(S_cur, cities);
    double cost_best = cost_cur;

    double T_init = -avg_delta_E_pos / log(params.P_init_accept);
    double T_final = params.T_final;
    double alpha = params.alpha;
    double K = params.K;
    int L = num_cities * params.L_multiplier;
    double T_transition = T_init * params.T_transition_factor;
    int iter = 0;
    double T = T_init;

    while (T > T_final) {
        iter++;
        for (int m = 0; m < L; m++) {
            vector<int> S_opt;
            double cost_change;
            double cost_opt;

            if (T > T_transition) { // With high T, use double bridge.
                pair<vector<int>, double> result = double_bridge_move(S_cur, cities, gen);
                S_opt = result.first;
                cost_change = result.second;
                cost_opt = cost_cur + cost_change;
            } else {
                pair<vector<int>, double> result = three_opt(S_cur, cities, gen);
                S_opt = result.first;
                cost_change = result.second;
                cost_opt = cost_cur + cost_change;
            }
            // 2-opt part
            // else if (T <= T_low_transition) {
            //     int i = dis_idx(gen);
            //     int j = dis_idx(gen);

            //     if (i > j)
            //         swap(i, j);
            //     if (i == j || (i == 1 && j == num_cities - 1))
            //         continue;
            //     // reverse S_opt[i ... j]
            //     S_opt = S_cur;
            //     reverse(S_opt.begin() + i, S_opt.begin() + j + 1);
            //     S_opt[num_cities] = S_opt[0];

            //     double added_cost = euc_2D(cities[S_cur[i - 1]], cities[S_opt[i]]) + euc_2D(cities[S_opt[j]], cities[S_cur[j + 1]]);
            //     double removed_cost = euc_2D(cities[S_cur[i - 1]], cities[S_cur[i]]) + euc_2D(cities[S_cur[j]], cities[S_cur[j + 1]]);
            //     cost_change = added_cost - removed_cost;
            //     cost_opt = cost_cur + cost_change;
            // }
            if (cost_change < 0) {
                S_cur = S_opt;
                cost_cur = cost_opt;
                if (cost_cur < cost_best){
                    S_best = S_cur;
                    cost_best = cost_cur;
                    // cout << " best: " << fixed << setprecision(2) << cost_best << " at T=" << T << endl;
                }
            } else {
                if (dis_prob(gen) < exp(-cost_change / (K * T))){
                    S_cur = S_opt;
                    cost_cur = cost_opt;
                }
            }
        }
        T = T * alpha;
        // if (iter % 500 == 0) {
        //      cout << "Iter: " << iter << " T: " << fixed << setprecision(4) << T
        //           << ", Cost_cur: " << cost_cur << ", Cost_best: " << cost_best << endl;
        // }
    }
    return S_best;
}

double calculate_avg_delta_E_positive(const vector<int>& initial_circuit, const vector<City>& cities, int num_iterations = 1000) {
    int num_cities = initial_circuit.size() - 1;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis_idx(1, num_cities - 1);

    double cost_cur = calculate_cost(initial_circuit, cities);
    double sum_delta_E = 0.0;
    int count_positive = 0;

    for (int k = 0; k < num_iterations; ++k) {
        int i = dis_idx(gen);
        int j = dis_idx(gen);

        if (i > j) swap(i, j);

        if (i == j || (i == 1 && j == num_cities - 1)) continue;

        vector<int> S_opt = initial_circuit;
        reverse(S_opt.begin() + i, S_opt.begin() + j + 1);
        S_opt[num_cities] = S_opt[0];

        double cost_opt = calculate_cost(S_opt, cities);
        double Delta_E = cost_opt - cost_cur;

        if (Delta_E > 0) {
            sum_delta_E += Delta_E;
            count_positive++;
        }
    }

    double avg_delta_E = (count_positive > 0) ? (sum_delta_E / count_positive) : 0.0;
    return avg_delta_E;
}

void two_approximation(vector<vector<pair<int, double>>>& adj_list_2_aprx) {
    for (auto& row_vec : adj_list_2_aprx) {
        row_vec.insert(row_vec.end(), row_vec.begin(), row_vec.end());
    }
}

pair<double, double> measure_2approx_runtime(const vector<City>& cities) {
    auto start = high_resolution_clock::now();
    
    vector<vector<pair<int, double>>> adj_list_mst;
    double mst_w = prim_mst(cities, adj_list_mst);
    
    vector<vector<pair<int, double>>> adj_list_2_aprx = adj_list_mst;
    two_approximation(adj_list_2_aprx);
    vector<int> euler_circ_2_aprx = euler_circuit(adj_list_2_aprx);
    vector<int> hamil_circ_2_aprx = hamilton_circuit(euler_circ_2_aprx, cities.size());
    double final_cost_2_aprx = calculate_cost(hamil_circ_2_aprx, cities);
    
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    
    return {final_cost_2_aprx, duration.count() / 1000000.0};
}

pair<double, double> measure_christofides_runtime(const vector<City>& cities) {
    auto start = high_resolution_clock::now();
    
    vector<vector<pair<int, double>>> adj_list_mst;
    double mst_w = prim_mst(cities, adj_list_mst);

    vector<int> odd_degree_vertices_indices;
    for(int i=0; i < cities.size(); ++i) {
        if(adj_list_mst[i].size() % 2 != 0) {
            odd_degree_vertices_indices.push_back(i);
        }
    }

    pair<vector<pair<int, int>>, double> matching_result = greedy_perfect_matching(odd_degree_vertices_indices, cities);
    vector<pair<int, int>> matched_pairs = matching_result.first;
    
    vector<vector<pair<int, double>>> adj_list_pm = adj_list_mst;
    for (const auto& pm_edge : matched_pairs) {
        int u = pm_edge.first;
        int v = pm_edge.second;
        double weight = euc_2D(cities[u], cities[v]);
        adj_list_pm[u].push_back({v, weight});
        adj_list_pm[v].push_back({u, weight});
    }
    
    vector<int> euler_circ_pm = euler_circuit(adj_list_pm);
    vector<int> hamil_circ_pm = hamilton_circuit(euler_circ_pm, cities.size());
    double final_cost_pm = calculate_cost(hamil_circ_pm, cities);
    
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    
    return {final_cost_pm, duration.count() / 1000000.0}; // Return time in seconds
}

pair<double, double> measure_SA_runtime(const vector<City>& cities, const SAParameters& params) {
    auto start = high_resolution_clock::now();
    
    vector<vector<pair<int, double>>> adj_list_mst;
    double mst_w = prim_mst(cities, adj_list_mst);
    
    vector<int> odd_degree_vertices_indices;
    for(int i=0; i < cities.size(); ++i) {
        if(adj_list_mst[i].size() % 2 != 0) {
            odd_degree_vertices_indices.push_back(i);
        }
    }
    
    pair<vector<pair<int, int>>, double> matching_result = greedy_perfect_matching(odd_degree_vertices_indices, cities);
    vector<pair<int, int>> matched_pairs = matching_result.first;
    
    vector<vector<pair<int, double>>> adj_list_pm = adj_list_mst;
    for (const auto& pm_edge : matched_pairs) {
        int u = pm_edge.first;
        int v = pm_edge.second;
        double weight = euc_2D(cities[u], cities[v]);
        adj_list_pm[u].push_back({v, weight});
        adj_list_pm[v].push_back({u, weight});
    }
    
    vector<int> euler_circ_pm = euler_circuit(adj_list_pm);
    vector<int> hamil_circ_pm = hamilton_circuit(euler_circ_pm, cities.size());
    
    double avg_delta_E_pos = calculate_avg_delta_E_positive(hamil_circ_pm, cities, 1000);
    vector<int> SA3opt_circ = SA_3opt(hamil_circ_pm, cities, avg_delta_E_pos, params);
    double final_cost_SA = calculate_cost(SA3opt_circ, cities);
    
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    
    return {final_cost_SA, duration.count() / 1000000.0}; // Return time in seconds
}

void measure_algorithm_runtimes(const string& filename, int num_runs, const map<string, SAParameters>& all_sa_params, const map<string, double>& optimal_solutions) {
    vector<City> cities = load_tsp(filename);
    
    if (cities.empty()) {
        cerr << "Failed to load cities or file is empty." << endl;
        return;
    }

    double optimal_cost = -1.0;
    if (optimal_solutions.count(filename)) {
        optimal_cost = optimal_solutions.at(filename);
    } else {
        cout << "Warning: Optimal solution cost not found for " << filename 
             << ". Approximation ratios will not be calculated." << endl;
    }
    
    cout << "\n=== Runtime Measurement for " << filename << " ===" << endl;
    cout << "Dataset size: " << cities.size() << " cities" << endl;
    cout << "Known Optimal Cost: " << fixed << setprecision(2) << optimal_cost << endl;
    cout << "Number of runs per algorithm: " << num_runs << endl << endl;
    
    cout << "\n--- 2-Approximation Algorithm ---" << endl;
    double total_2approx_time = 0.0;
    double total_2approx_cost = 0.0;
    double best_2approx_cost = numeric_limits<double>::max();

    for (int i = 0; i < num_runs; i++) {
        pair<double, double> result = measure_2approx_runtime(cities);
        double current_cost = result.first;
        double runtime = result.second;
        
        total_2approx_cost += current_cost;
        best_2approx_cost = min(best_2approx_cost, current_cost);
        total_2approx_time += runtime;
        cout << "  Run " << i + 1 << ": Cost = " << fixed << setprecision(2) << current_cost 
             << ", Time = " << setprecision(6) << runtime << " s" << endl;
    }
    double avg_2approx_time = total_2approx_time / num_runs;
    double avg_2approx_cost = total_2approx_cost / num_runs;
    
    cout << "\n--- Christofides (Greedy Matching) ---" << endl;
    double total_christofides_time = 0.0;
    double total_christofides_cost = 0.0;
    double best_christofides_cost = numeric_limits<double>::max();

    for (int i = 0; i < num_runs; i++) {
        pair<double, double> result = measure_christofides_runtime(cities);
        double current_cost = result.first;
        double runtime = result.second;

        total_christofides_cost += current_cost;
        best_christofides_cost = min(best_christofides_cost, current_cost);
        total_christofides_time += runtime;
        cout << "  Run " << i + 1 << ": Cost = " << fixed << setprecision(2) << current_cost 
             << ", Time = " << setprecision(6) << runtime << " s" << endl;
    }
    double avg_christofides_time = total_christofides_time / num_runs;
    double avg_christofides_cost = total_christofides_cost / num_runs;
    
    cout << "\n--- SA with 3-Opt (initialized by Christofides) ---" << endl;
    SAParameters current_sa_params;
    auto it_sa = all_sa_params.find(filename);
    if (it_sa != all_sa_params.end()) {
        current_sa_params = it_sa->second;
    } else {
        current_sa_params = {0.5, 0.1, 0.9999, 0.6, 10, 0.6}; 
        cout << "Warning: Using default SA parameters for " << filename << endl;
    }

    double total_SA_time = 0.0;
    double total_SA_cost = 0.0;
    double best_SA_cost = numeric_limits<double>::max();

    for (int i = 0; i < num_runs; i++) {
        cout << "  SA Run " << i + 1 << "/" << num_runs << " progress (T_init/best updates from SA_3opt):" << endl;
        pair<double, double> result = measure_SA_runtime(cities, current_sa_params);
        double current_cost = result.first;
        double runtime = result.second;

        total_SA_cost += current_cost;
        best_SA_cost = min(best_SA_cost, current_cost);
        total_SA_time += runtime;
        cout << "  -> SA Run " << i + 1 << " completed: Cost = " << fixed << setprecision(2) << current_cost
             << ", Time = " << setprecision(3) << runtime << "s" << endl;
    }
    double avg_SA_time = total_SA_time / num_runs;
    double avg_SA_cost = total_SA_cost / num_runs;
    
    // --- Display summary ---
    cout << "\n\n====== SUMMARY for " << filename << " ======" << endl;
    if (optimal_cost > 0) {
        cout << "Optimal Cost: " << fixed << setprecision(2) << optimal_cost << endl;
    }
    cout << fixed << setprecision(4);
    
    cout << "\n--- 2-Approximation Algorithm ---" << endl;
    cout << "Avg Cost:   " << avg_2approx_cost;
    if (optimal_cost > 0) cout << " (Ratio: " << (avg_2approx_cost / optimal_cost) << ")";
    cout << endl;
    cout << "Best Cost:  " << best_2approx_cost;
    if (optimal_cost > 0) cout << " (Ratio: " << (best_2approx_cost / optimal_cost) << ")";
    cout << endl;
    cout << "Avg Time:   " << setprecision(6) << avg_2approx_time << " s" << endl;

    cout << "\n--- Christofides (Greedy Matching) ---" << endl;
    cout << "Avg Cost:   " << avg_christofides_cost;
    if (optimal_cost > 0) cout << " (Ratio: " << (avg_christofides_cost / optimal_cost) << ")";
    cout << endl;
    cout << "Best Cost:  " << best_christofides_cost;
    if (optimal_cost > 0) cout << " (Ratio: " << (best_christofides_cost / optimal_cost) << ")";
    cout << endl;
    cout << "Avg Time:   " << setprecision(6) << avg_christofides_time << " s" << endl;

    cout << "\n--- SA with 3-Opt ---" << endl;
    cout << "Avg Cost:   " << avg_SA_cost;
    if (optimal_cost > 0) cout << " (Ratio: " << (avg_SA_cost / optimal_cost) << ")";
    cout << endl;
    cout << "Best Cost:  " << best_SA_cost;
    if (optimal_cost > 0) cout << " (Ratio: " << (best_SA_cost / optimal_cost) << ")";
    cout << endl;
    cout << "Avg Time:   " << setprecision(6) << avg_SA_time << " s" << endl;
    cout << "====================================" << endl;
}

int main() {
    // struct SAParameters {
    //     double P_init_accept;
    //     double T_final;
    //     double alpha;
    //     double K;
    //     int L_multiplier;
    //     double T_transition_factor;
    // };
    // Parameters for a280.tsp

    map<string, SAParameters> dataset_sa_params;
    dataset_sa_params["a280.tsp"] = {
        0.5,
        0.1,
        0.9995,
        0.5,
        10, 
        0.5
    };

    // Parameters for xql662.tsp
    dataset_sa_params["xql662.tsp"] = {
        0.5,
        0.1,
        0.9995,
        1,
        10, 
        0.5
    };

    // Parameters for kz9976.tsp
    dataset_sa_params["kz9976.tsp"] = {
        0.5,
        0.1,
        0.99,
        0.003,
        5,
        0.5
    };
    
    // Parameters for mona-lisa100K.tsp
    dataset_sa_params["mona-lisa100K.tsp"] = {
        0.5,
        0.1,
        0.95,
        0.01,
        0.5,
        0.5
    };

    map<string, double> optimal_solutions;
    optimal_solutions["a280.tsp"] = 2579;
    optimal_solutions["xql662.tsp"] = 2513;
    optimal_solutions["kz9976.tsp"] = 1061881;
    optimal_solutions["mona-lisa100K.tsp"] = 5757191;
    
    vector<string> test_files = {"a280.tsp", "xql662.tsp", "kz9976.tsp", "mona-lisa100K.tsp"};
    int num_runs = 10;
    
    for (const string& filename : test_files) {
        measure_algorithm_runtimes(filename, num_runs, dataset_sa_params, optimal_solutions);
        cout << "\n" << string(50, '-') << "\n" << endl;
    }
    
    return 0;
}