#include <vector>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <CALAtmSim.hpp>
#include <math_rng.hpp>
#include <sys_utils.hpp>
#include <assert.h>

using namespace std;

int main(){

    int dim = 10;

	cal::AlignedVector<double> *ts = new cal::AlignedVector<double> (dim);
	cal::AlignedVector<double> *az = new cal::AlignedVector<double> (dim);
	cal::AlignedVector<double> *el = new cal::AlignedVector<double> (dim);
    cal::AlignedVector<double> *diff = new cal::AlignedVector<double> (dim);

    cal::rng_dist_normal(dim, 1, dim, 0, dim, ts->data());
	cal::rng_dist_normal(dim, int(dim/2), dim+int(dim/2), 0, dim, az->data());
    cal::rng_dist_normal(dim, 1, dim, 0, dim, el->data());

    set_difference(ts->begin(), ts->end(), el->begin(), el->end(), inserter(*diff, diff->begin()));

	for(int i=0; i<dim; i++){
        cout << ts->at(i) + 5 << "\t" << az->at(i) + 5 << "\t" << el->at(i) + 5<< "\t" << diff->at(i) <<  endl;
        try{
            assert (diff->at(i) == 0.0);
        } catch (int e) {
            cout << "An exception occurred. Exception Nr. " << e << endl;
            return 1;
        }
    }

    free(ts);
    free(az);
    free(el);

	return 0;
}
