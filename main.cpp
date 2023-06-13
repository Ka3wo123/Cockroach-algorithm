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
        sum += x[i-1] * sin(sqrt(fabs(x[i-1])));
    }

    return 418.9829 * dim - sum;
}


vector<double> generateRandomSolution(int dim, double lowerBound, double upperBound) {
    random_device rd;
    mt19937 random(rd());
    uniform_real_distribution<> bounds{lowerBound, upperBound};
    vector<double> solution;


    for (int i = 0; i < dim; i++) {
        solution.push_back(bounds(random));
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

vector<double> updatePosition(vector<double> position, vector<double> optimumPosition, double stepSize) {
    random_device rd;
    mt19937 random(rd());
    uniform_real_distribution<> dis{0, 1};
    vector<double> newPosition(position.size());

    for (int i = 0; i < position.size(); i++) {
        double rand = dis(random);
        double delta = rand * (optimumPosition[i] - position[i]);
        newPosition[i] = position[i] + stepSize * delta;
    }
    return newPosition;
}

bool isLocalOptimum(vector<double> localOpt, vector<double> pos) {
    for(int i = 0; i < localOpt.size(); i++) {
        if (localOpt[i] != pos[i]) return false;
    }

    return true;
}

double diffCocks(vector<double> cockI, vector<double> cockJ) {
    double sum = 0;

    for(int i = 0; i < cockI.size(); i++) {
        sum += fabs(cockI[i] - cockJ[i]) * fabs(cockI[i] - cockJ[i]);
    }

    return sqrt(sum);
}
vector<double> updateInLightIter(vector<double> pos) {
    random_device rd;
    mt19937 random(rd());
    uniform_real_distribution<> dis(0, 1);

    vector<double> newPos;

    for(int i = 0; i < pos.size(); i++) {
        newPos.push_back(pos[i] + dis(random));
    }

    return newPos;
}


vector<double> cockroachAlgorithm(int numCockroach, int dim, int maxIter, double lowerBound, double upperBound, double visual, double eps, double w, double (*testFunction)(vector<double>, int)) {
    vector<double> globalOptimum;
    vector<double> localOptimum;
    vector<vector<double>> cockroaches;
    random_device rd;
    mt19937 random(rd());
    uniform_real_distribution<> dis(1, (double) numCockroach);
    uniform_int_distribution<> light(1, maxIter);

    int lightIteration;

    for (int i = 0; i < numCockroach; i++) {
        cockroaches.push_back(generateRandomSolution(dim, lowerBound, upperBound));
    }

    globalOptimum = findGlobalOptimum(cockroaches, testFunction);


    for(int t = 0; t < maxIter; t++) {        
        if(testFunction(globalOptimum, dim) <= eps) {
            return globalOptimum;
        }

        lightIteration = light(random);

        for(int i = 0; i < numCockroach; i++) {            
            for(int j = 0; j < numCockroach; j++) {
                if(diffCocks(cockroaches[i], cockroaches[j]) < visual && (testFunction(cockroaches[j], dim) < testFunction(cockroaches[i], dim))) {
                    localOptimum = cockroaches[j];
                }
            }

            if(isLocalOptimum(localOptimum, cockroaches[i])) {
                cockroaches[i] = updatePosition(cockroaches[i], globalOptimum, w);
            } else {                
                cockroaches[i] = updatePosition(cockroaches[i], localOptimum, w);
            }

            if(testFunction(cockroaches[i], dim) < testFunction(globalOptimum, dim)) {
                globalOptimum = cockroaches[i];
            }           
        }                

        if(t == lightIteration) {

            for(int i = 0; i < numCockroach; i++) {
                cockroaches[i] = updateInLightIter(cockroaches[i]);


                if(testFunction(cockroaches[i], dim) < testFunction(globalOptimum, dim)) {
                    globalOptimum = cockroaches[i];
                }
            }
        }

        int k = dis(random);

        if(!isLocalOptimum(globalOptimum, cockroaches[k])) {
            cockroaches[k] = globalOptimum;
        }

        cout << "Iteration " << t << " ";
        for(int i = 0; i < globalOptimum.size(); i++) {
            cout << globalOptimum[i] << " ";
        }
        cout << endl;
    }

    return globalOptimum;

}


int main()
{

    const int NUM_OF_CROCKROACHES = 260;
    const int MAX_ITER = 200;
    const double EPS = 0.001;
    const double VISUAL = 0.1;
    const int DIM = 2;
    const double W = 0.1;
    const double LOWER_BOUND = -500;
    const double UPPER_BOUND = 500;

    vector<double> solution = cockroachAlgorithm(NUM_OF_CROCKROACHES, DIM, MAX_ITER, LOWER_BOUND, UPPER_BOUND, VISUAL, EPS, W, &schwefel);

    cout << "test";
    cout << "Znaleziono rozwiazanie: " << endl;
    cout << "F(x) = " << schwefel(solution, DIM) << " dla x { ";
    for(auto x : solution) {
        cout << x << " ";
    }

    cout << "}" << endl;










}
