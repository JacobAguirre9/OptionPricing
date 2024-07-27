// This script is an implementation of the European binomial option pricing model

// The main difference between the American and European binomial option pricing models is that the American model allows for early exercise of the option, while the European model does not.

// Include the necessary libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

// Function to calculate the binomial coefficient
double binomial_coefficient(int n, int k) {
    double result = 1.0;
    for (int i = 1; i <= k; i++) {
        result = result * (n - k + i) / i;
    }
    return result;
}

// Function to calculate the European option price
double european_option_price(double S, double K, double r, double sigma, double T, int n, char type) {
    double dt = T / n;
    double u = exp(sigma * sqrt(dt));
    double d = 1.0 / u;
    double p = (exp(r * dt) - d) / (u - d);
    double q = 1.0 - p;
    vector<double> prices(n + 1);
    vector<double> option_values(n + 1);
    for (int i = 0; i <= n; i++) {
        prices[i] = S * pow(u, n - i) * pow(d, i);
        if (type == 'C') {
            option_values[i] = max(0.0, prices[i] - K);
        } else if (type == 'P') {
            option_values[i] = max(0.0, K - prices[i]);
        }
    }
    for (int j = n - 1; j >= 0; j--) {
        for (int i = 0; i <= j; i++) {
            option_values[i] = (p * option_values[i] + q * option_values[i + 1]) * exp(-r * dt);
        }
    }
    return option_values[0];
}

// Compute statistics
double mean(vector<double> values) {
    double sum = 0.0;
    for (int i = 0; i < values.size(); i++) {
        sum += values[i];
    }
    return sum / values.size();
}

double standard_deviation(vector<double> values) {
    double avg = mean(values);
    double sum = 0.0;
    for (int i = 0; i < values.size(); i++) {
        sum += pow(values[i] - avg, 2);
    }
    return sqrt(sum / values.size());
}

int main() {
    double S = 100.0; // Stock price
    double K = 100.0; // Strike price
    double r = 0.05; // Risk-free rate
    double sigma = 0.2; // Volatility
    double T = 1.0; // Time to maturity
    int n = 100; // Number of time steps
    char type = 'C'; // Option type (C for call, P for put)
    int num_simulations = 1000; // Number of simulations
    vector<double> option_prices;
    for (int i = 0; i < num_simulations; i++) {
        double option_price = european_option_price(S, K, r, sigma, T, n, type);
        option_prices.push_back(option_price);
    }
    double option_price_mean = mean(option_prices);
    double option_price_stddev = standard_deviation(option_prices);
    cout << "European option price: " << fixed << setprecision(2) << option_price_mean << endl;
    cout << "Standard deviation: " << fixed << setprecision(2) << option_price_stddev << endl;
    return 0;
}

// Export results to a CSV file and generate plots
#include <fstream>
#include <sstream>
#include <string>

void export_to_csv(vector<double> values, string filename) {
    ofstream file(filename);
    for (int i = 0; i < values.size(); i++) {
        file << values[i] << endl;
    }
    file.close();
}

// Full example
int main() {
    double S = 100.0; // Stock price
    double K = 100.0; // Strike price
    double r = 0.05; // Risk-free rate
    double sigma = 0.2; // Volatility
    double T = 1.0; // Time to maturity
    int n = 100; // Number of time steps
    char type = 'C'; // Option type (C for call, P for put)
    int num_simulations = 1000; // Number of simulations
    vector<double> option_prices;
    for (int i = 0; i < num_simulations; i++) {
        double option_price = european_option_price(S, K, r, sigma, T, n, type);
        option_prices.push_back(option_price);
    }
    double option_price_mean = mean(option_prices);
    double option_price_stddev = standard_deviation(option_prices);
    cout << "European option price: " << fixed << setprecision(2) << option_price_mean << endl;
    cout << "Standard deviation: " << fixed << setprecision(2) << option_price_stddev << endl;
    export_to_csv(option_prices, "option_prices.csv");
    return 0;
}