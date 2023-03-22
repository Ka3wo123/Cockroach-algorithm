#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <ctime>

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

double schwefel(double x, int dim) {
    double sum = 0.0;
    for(int i = 1; i <= dim; i++) {
        sum += -x * sin(sqrt(fabs(x)));
    }

    return -sum;
}

vector<double> updatePosition() {

}

vector<double> findGlobalBest(vector<double> cockroaches, int dim, double (*func)(double)) {
    vector<double> globalBest;
    globalBest[0] = cockroaches[0];
    for(int i = 1; i < dim; i++) {
        if(func(cockroaches[i]) < func(globalBest[0])) {
            globalBest[0] = cockroaches[i];
        }
    }

    return globalBest;
}

vector<double> cockroachAlgorithm(int numCockroaches, int maxIter, double eps, double lowerBound, double upperBound) {
    vector<double> cockroachesPosition(numCockroaches);
}



int main()
{



}
