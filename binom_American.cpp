// This script is an implementation of the American binomial option pricing model

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

// Function to calculate the American option price
double american_option_price(double S, double K, double r, double sigma, double T, int n, char type) {
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
            if (type == 'C') {
                prices[i] = S * pow(u, j - i) * pow(d, i);
                option_values[i] = max(option_values[i], prices[i] - K);
            } else if (type == 'P') {
                prices[i] = S * pow(u, j - i) * pow(d, i);
                option_values[i] = max(option_values[i], K - prices[i]);
            }
        }
    }
    return option_values[0];
}

int main() {
    double S = 100.0; // Stock price
    double K = 100.0; // Strike price
    double r = 0.05; // Risk-free rate
    double sigma = 0.2; // Volatility
    double T = 1.0; // Time to maturity
    int n = 100; // Number of time steps
    char type = 'C'; // Option type (C for call, P for put)
    double option_price = american_option_price(S, K, r, sigma, T, n, type);
    cout << "American option price: " << fixed << setprecision(2) << option_price << endl;
    return 0;
}