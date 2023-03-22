#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <ctime>
#include <cstdlib>
#define PI 3.14159265

using namespace std;

double hiperElipsoide(vector<double> x, int dim) {
    double sum = 0.0;
    for(int i = 1; i <= dim; i++) {
        sum += i * (x[i-1] * x[i-1]);
    }
    return sum;
}

double rastring(vector<double> x, int dim) {
    double sum = 0.0;
    for(int i = 1; i <= dim; i++) {
        sum += (x[i-1] * x[i-1]) - 10 * cos(PI * 2 * x[i-1]);
    }

    return 10 * dim + sum;
}

double schwefel(vector<double> x, int dim) {
    double sum = 0.0;
    for(int i = 1; i <= dim; i++) {
        sum += -x[i-1] * sin(sqrt(fabs(x[i-1])));
    }

    return -sum;
}


vector<double> generateRandomSolution(int dim, double lowerBound, double upperBound) {
    mt19937 random(time(NULL));
    uniform_real_distribution<> bounds(lowerBound, upperBound);
    vector<double> solution(dim);

    for (int i = 0; i < dim; i++) {
        solution[i] = bounds(random);
    }

    return solution;
}


vector<double> findGlobalOptimum(vector<vector<double>> cockroaches, double (*testFunction)(vector<double>, int)) {
    vector<double> globalOptimum;
    int dim = cockroaches[0].size();
    globalOptimum = cockroaches[0];

    for(int i = 1; i < cockroaches.size(); i++) {
        if(testFunction(cockroaches[i], dim) < testFunction(globalOptimum, dim)) {
            globalOptimum = cockroaches[i];
        }
    }

    return globalOptimum;
}


vector<double> cockroachAlgorithm(int numCockroach, int dim, int maxIter, double lowerBound, double upperBound, double visual, double eps, double w, double (*testFunction)(vector<double>, int)) {
    vector<double> globalOptimum;
    vector<double> localOptimum;
    vector<vector<double>> cockroaches(numCockroach);

    for (int i = 0; i < numCockroach; i++) {
        cockroaches[i] = generateRandomSolution(dim, lowerBound, upperBound);
    }

    globalOptimum = findGlobalOptimum(cockroaches, testFunction);

    for(int t = 0; t < maxIter; t++) {
        if(testFunction(globalOptimum, dim) <= eps) {
            return globalOptimum;
        }
        for(int i = 0; i < numCockroach; i++) {
            for(int j = 0; j < numCockroach; j++) {
                if((fabs(cockroaches[t][i] - cockroaches[t][j]) < visual) && (testFunction(cockroaches[i], dim) < testFunction(cockroaches[j], dim))) {
                    localOptimum = cockroaches[j];
                }
            }

            if(localOptimum[i] == cockroaches[t][i]) {
                cockroaches[t][i] = w * cockroaches[t][i] + i + rand() % 2 * (globalOptimum[i] - cockroaches[t][i]);
            } else {
                cockroaches[t][i] = w * cockroaches[t][i] + i + rand() % 2 * (localOptimum[i] - cockroaches[t][i]);
            }

            if(testFunction(cockroaches[i], dim) < testFunction(globalOptimum, dim)) {
                globalOptimum = cockroaches[i];
            }
        }

        if(t == 4 || t == 10 || t == 16 || t == 50 || t == 89) {

            for(int i = 0; i < numCockroach; i++) {
                cockroaches[t][i] = cockroaches[t][i] + rand() % 2 + 1;

                if(testFunction(cockroaches[i], dim) < testFunction(globalOptimum, dim)) {
                    globalOptimum = cockroaches[i];
                }
            }
        }

        int k = rand() % numCockroach + 1;

        if(cockroaches[k] != globalOptimum) {
            cockroaches[k] = globalOptimum;
        }
    }

    return globalOptimum;


}



int main()
{

    const int NUM_OF_CROCKROACHES = 500;
    const int MAX_ITER = 1000;
    const double EPS = 0.001;
    const double VISUAL = 0.5;
    const int DIM = 2;
    const double W = 0.1;
    const double LOWER_BOUND = -5.12;
    const double UPPER_BOUND = 5.12;

    vector<double> solution = cockroachAlgorithm(NUM_OF_CROCKROACHES, DIM, MAX_ITER, LOWER_BOUND, UPPER_BOUND, VISUAL, EPS, W, &schwefel);


    cout << "Znaleziono rozwiazanie: " << endl;
    cout << "F(x) = " << hiperElipsoide(solution, DIM) << " dla x { ";
    for(auto x : solution) {
        cout << x << " ";
    }

    cout << "}" << endl;
    ;





}
