#include <iostream>
#include <vector>
#include <cmath>

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

int main()
{

}
