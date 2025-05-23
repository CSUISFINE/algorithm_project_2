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

pair<vector<int>, double> held_karp (const vector<City>& cities) {
    const double INF = numeric_limits<double>::infinity();
    int num_cities = cities.size();
    int start = 0;
    vector<vector<double>> dp((1 << num_cities), vector<double>(num_cities, INF));
    dp[1 << start][start] = 0;

    for (int mask = 1; mask < (1 << num_cities); mask++) {
        if ((mask & (1 << start)) == 0)
            continue;
        for (int j = 0; j < num_cities; j++) {
            if (((mask & (1 << j)) != 0)) {
                if (mask == (1 << start) && j == start)
                    continue;
                int prev_mask = mask ^ (1 << j);
                if (prev_mask == 0)
                    continue;
                for (int k = 0; k < num_cities; k++) {
                    if ((prev_mask & (1 << k)) != 0 && dp[prev_mask][k] < INF){
                        dp[mask][j] = min(dp[mask][j], dp[prev_mask][k] + euc_2D(cities[k], cities[j]));
                    }
                }
            }
        }
    }
    int final_mask = (1 << num_cities) - 1;
    double tsp_cost = INF;
    int last_city;

    for (int j = 0; j < num_cities; j++) {
        if (j == start)
            continue;
        if (dp[final_mask][j] < INF) {
            double cost_js = euc_2D(cities[j], cities[start]);
            if (dp[final_mask][j] + cost_js < tsp_cost) {
                tsp_cost = dp[final_mask][j] + cost_js;
                last_city = j;
            }
        }
    }

    vector<int> tsp_path;
    int cur_city = last_city;
    int cur_mask = final_mask;
    int prev_city;
    const double EPSILON = 1e-9;
    tsp_path.push_back(cur_city);
    while(cur_mask != (1 << start)){
        int prev_mask = cur_mask ^ (1 << cur_city);
        int prev_city = -1;
        for (int k = 0; k < num_cities; k++) {
            if ((prev_mask & (1 << k)) && dp[prev_mask][k] < INF) {
                if (abs(dp[cur_mask][cur_city] - (dp[prev_mask][k] + euc_2D(cities[k], cities[cur_city]))) < EPSILON) {
                    prev_city = k;
                    break;
                }
            }
        }
        cur_city = prev_city;
        cur_mask = prev_mask;
        tsp_path.push_back(cur_city);
    }
    reverse(tsp_path.begin(), tsp_path.end());
    return {tsp_path, tsp_cost};
}

int main() {
    string filename = "test.tsp";
    vector<City> cities = load_tsp(filename);

    if (cities.empty()) {
        cerr << "Failed to load cities or file is empty." << endl;
        return 1;
    }

    cout << "Loaded " << cities.size() << " cities." << endl;
    pair<vector<int>, double> tsp = held_karp(cities);
    cout << "tsp path: ";
    for (auto city : tsp.first) {
        cout << city << " ";
    }
    cout << endl;
    cout << "tsp cost: " << tsp.second << endl;
}