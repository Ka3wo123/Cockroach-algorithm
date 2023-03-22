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

double schwefel(vector<double> x, int dim) {
    double sum = 0.0;
    for(int i = 1; i <= dim; i++) {
        sum += -x[i-1] * sin(sqrt(fabs(x[i-1])));
    }

    return -sum;
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

vector<double> updatePosition() {

}

vector<double> cockroachAlgorithm() {


}



int main()
{

    const int NUM_OF_CROCKROACHES = 500;
    const int MAX_ITER = 1000;
    const double EPS = 0.001;
    const double VISIBILITY = 0.5;
    const double W = 0.1;






}
