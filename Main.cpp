#include <iostream>
#include <vector>
#include <cmath>
#include <random>

// Define constants
const double PI = 3.14159265358979323846;

// Function to generate a random normal variable using Box-Muller transform
double randomNormal(double mean = 0.0, double stddev = 1.0) {
    static std::mt19937 gen{std::random_device{}()};
    static std::normal_distribution<> dist(mean, stddev);
    return dist(gen);
}

// Pricing function for a simple European call option using Black-Scholes
double callOptionPrice(double S, double K, double T, double r, double sigma) {
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    double N_d1 = 0.5 * (1.0 + erf(d1 / sqrt(2.0)));
    double N_d2 = 0.5 * (1.0 + erf(d2 / sqrt(2.0)));
    return S * N_d1 - K * exp(-r * T) * N_d2;
}

// Monte Carlo simulation to estimate the option price
double monteCarloCallOptionPrice(double S1, double K1, double S2, double K2,
                                  double T, double r, double sigma1, double sigma2,
                                  double correlation, int numSimulations) {
    std::vector<double> prices(numSimulations);
    double sqrtT = sqrt(T);

    for (int i = 0; i < numSimulations; ++i) {
        double z1 = randomNormal();
        double z2 = correlation * z1 + sqrt(1 - correlation * correlation) * randomNormal();

        double S1_T = S1 * exp((r - 0.5 * sigma1 * sigma1) * T + sigma1 * sqrtT * z1);
        double S2_T = S2 * exp((r - 0.5 * sigma2 * sigma2) * T + sigma2 * sqrtT * z2);

        double payoff = std::max(S1_T - K1, 0.0) + std::max(S2_T - K2, 0.0);
        prices[i] = exp(-r * T) * payoff;
    }

    double meanPrice = 0.0;
    for (double price : prices) {
        meanPrice += price;
    }
    meanPrice /= numSimulations;

    return meanPrice;
}

// Objective function to minimize (squared difference between model and market prices)
double objectiveFunction(double correlation, double marketPrice,
                         double S1, double K1, double S2, double K2,
                         double T, double r, double sigma1, double sigma2,
                         int numSimulations) {
    double modelPrice = monteCarloCallOptionPrice(S1, K1, S2, K2, T, r, sigma1, sigma2, correlation, numSimulations);
    return (modelPrice - marketPrice) * (modelPrice - marketPrice);
}

// Simple optimization using grid search
double optimizeCorrelation(double marketPrice,
                           double S1, double K1, double S2, double K2,
                           double T, double r, double sigma1, double sigma2,
                           int numSimulations) {
    double bestCorrelation = 0.0;
    double minError = std::numeric_limits<double>::max();
    
    for (double correlation = -1.0; correlation <= 1.0; correlation += 0.01) {
        double error = objectiveFunction(correlation, marketPrice, S1, K1, S2, K2, T, r, sigma1, sigma2, numSimulations);
        if (error < minError) {
            minError = error;
            bestCorrelation = correlation;
        }
    }

    return bestCorrelation;
}

int main() {
    // Example market data
    double marketPrice = 10.0; // Example market option price
    double S1 = 100.0; // Asset 1 price
    double K1 = 95.0;  // Asset 1 strike price
    double S2 = 200.0; // Asset 2 price
    double K2 = 190.0; // Asset 2 strike price

    double T = 1.0; // Time to maturity
    double r = 0.05; // Risk-free rate
    double sigma1 = 0.20; // Volatility of asset 1
    double sigma2 = 0.25; // Volatility of asset 2

    int numSimulations = 10000; // Number of Monte Carlo simulations

    double impliedCorrelation = optimizeCorrelation(marketPrice, S1, K1, S2, K2, T, r, sigma1, sigma2, numSimulations);

    std::cout << "Implied Correlation: " << impliedCorrelation << std::endl;

    return 0;
}
