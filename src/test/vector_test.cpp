#include <vector>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <CALAtmSim.hpp>
#include <math_rng.hpp>
#include <sys_utils.hpp>
#include <assert.h>
#include <limits>
#include <atmosphere.hpp>
using namespace std;

int main(){

    int dim = 10;

    cal::AlignedVector<double> *a = new cal::AlignedVector<double> (dim);
	cal::AlignedVector<double> *b = new cal::AlignedVector<double> (dim);
	cal::AlignedVector<double> *c = new cal::AlignedVector<double> (dim);
    cal::AlignedVector<double> *diff = new cal::AlignedVector<double> (dim);

    cal::rng_dist_normal(dim, 1, dim, 0, dim, a->data());
	cal::rng_dist_normal(dim, int(dim/2), dim+int(dim/2), 0, dim, b->data());
    cal::rng_dist_normal(dim, 1, dim, 0, dim, c->data());

    set_difference(a->begin(), a->end(), c->begin(), c->end(), inserter(*diff, diff->begin()));

    cout << "a_vec" << "\t" << "b_vec." << "\t" << "c_vec." << "\t" << "diff" << endl;
	for(int i=0; i<dim; i++){
        cout << a->at(i) + 5 << "\t" << b->at(i) + 5 << "\t" << c->at(i) + 5<< "\t" << diff->at(i) <<  endl;
        assert (diff->at(i) <= numeric_limits<double>::denorm_min());
    }

    free(a);
    free(b);
    free(c);

    cal::atm_sim *atm_strip = new cal::atm_sim(azmin, azmax, elmin, elmax,
                                               tmin,  tmax,
                                               lmin_center, lmin_sigma,
                                               lmax_center, lmax_sigma,
                                               w_center, w_sigma, wdir_center, wdir_sigma,
                                               z0_center, z0_sigma,
                                               T0_center, T0_sigma,
                                               zatm, zmax, xstep, ystep, zstep, nelem_sim_max, verbosity,
                                               key1, key2,
                                               counterval1, counterval2,
                                               cachedir,
                                               rmin, rmax);

    double fs_hz  = 20.0;                // Sample Freq. [Hz]
    double dt_sec = 1.0 / fs_hz;         // Time between two samples [sec.]
    double t_start = 0.0;                // Start time [sec.]
    double t_stop  = 60.0 + dt_sec;      // 1h of observation [sec.]

    int samples = int((t_stop - t_start) / dt_sec);

    double az_min = 0.0;
    double az_max = M_PI / 6.0;
    double ces_el = (70.0/180.0) * M_PI;
    double ss = (1.0 / 180.0) * M_PI;
    bool direction_lr = true;

    vector<double> ts;
	vector<double> az;
	vector<double> el;
    double *tod;
    tod = (double*)malloc(samples*sizeof(double));

    double t=0;
    for(int i=0; i<samples; i++){
        ts.push_back(t+(dt_sec*i));
        el.push_back(ces_el);
    }

    double d_az = ss / fs_hz;
    double az_now = az_min;
    for(int i=0; i<samples; i++){
        if(direction_lr){
            az_now += d_az;
            az.push_back(az_now);
            if(az_now >= az_max) direction_lr != direction_lr;
        } else {
            az_now -= d_az;
            az.push_back(az_now);
            if(az_now <= 0) direction_lr != direction_lr;
        }
    }

    atm_strip->simulate(false);

    //cout << ts.size() << " " << az.size() << " " << el.size() << " " <<tod.size() << " " << samples << endl;

    atm_strip->observe(ts.data(), az.data(), el.data(), tod, samples);

    for(int i=0; i<samples; i++)
        cout << ts.at(i) << '\t' << tod[i] << endl;

	return 0;
}
