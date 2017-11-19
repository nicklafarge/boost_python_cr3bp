#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/math/special_functions/sign.hpp>

using namespace std;
using namespace boost::numeric::odeint;
using namespace boost::python;

typedef std::vector<double> state_type;

// CR3BP Equations of Motion
class eomCR3BP {
    double m_mu;

public:
    eomCR3BP(double mu) : m_mu(mu) {}

    void operator()(const state_type &x, state_type &dxdt, const double /* t */) {

        // EOMs written in other boost integration entry on wiki
        double oneminusmu = 1 - m_mu;
        double d = sqrt(((x[0] + m_mu) * (x[0] + m_mu)) + (x[1] * x[1]) + (x[2] * x[2]));
        double r = sqrt(((x[0] - oneminusmu) * (x[0] - oneminusmu)) + (x[1] * x[1]) + (x[2] * x[2]));
        double r3 = r * r * r;
        double d3 = d * d * d;

        dxdt[0] = x[3];
        dxdt[1] = x[4];
        dxdt[2] = x[5];
        dxdt[3] = -(oneminusmu * (x[0] + m_mu) / d3) - (m_mu * (x[0] - oneminusmu) / r3) + (2 * x[4]) + x[0];
        dxdt[4] = -(oneminusmu * (x[1]) / d3) - (m_mu * (x[1]) / r3) - (2 * x[3]) + x[1];
        dxdt[5] = -(oneminusmu * (x[2]) / d3) - (m_mu * (x[2]) / r3);

        state_type pos(3);
        pos[0] = x[0];
        pos[1] = x[1];
        pos[2] = x[2];

    }

};

// Used to pull states out of the integration
struct getStateAndTime {
    std::vector<state_type> &m_states;
    std::vector<double> &m_times;

    getStateAndTime(std::vector<state_type> &states, std::vector<double> &times) : m_states(states), m_times(times) {}

    void operator()(const state_type &x, double t) {
        m_states.push_back(x);
        m_times.push_back(t);

    }

};

struct CR3BP {
    // Put function into format python can understand and the propagate
    boost::python::list propPy(boost::python::list &ic, boost::python::list &in_times, double mu, int state_dim,
                               int t_dim, double tol, double step_size) {

        typedef std::vector<double> state_type;
        std::vector<state_type> statesAndTimes;

        state_type IC(state_dim, 0);
        state_type t(t_dim, 0);

        // Transform inputs
        for (int i = 0; i < len(ic); ++i) {
            IC[i] = boost::python::extract<double>(ic[i]);
        }

        for (int i = 0; i < len(in_times); ++i) {
            t[i] = boost::python::extract<double>(in_times[i]);
        }

        // Propagate
        statesAndTimes = prop(IC, t, mu, state_dim, t_dim, tol, step_size);

        // Create python list from data to return
        return toTwoDimPythonList(statesAndTimes);
    }

    // Propagation function
    std::vector<vector<double >> prop(vector<double> ic, vector<double> t, double mu, int state_dim, int t_dim,
                                      double tol, double step_size) {
        using namespace std;
        using namespace boost::numeric::odeint;

        typedef std::vector<double> state_type;
        std::vector<state_type> statesAndTimes;

        // Set vectors intermediate steps during integration
        std::vector<double> tOut;
        std::vector<state_type> statesOut;

        // Determine step size (forward or backward) and set initial step size
        double h = t[1] > t[0] ? step_size : -step_size;

        // Set integrator type -> Currently set at rk78
        double relTol = tol;
        double absTol = tol;
        typedef runge_kutta_fehlberg78<state_type> rk78;
        auto stepper = make_controlled<rk78>(absTol, relTol);

        // Create eom to integrate
        eomCR3BP eom(mu);
        size_t steps = integrate_adaptive(stepper, eom, ic, t[0], t[1], h, getStateAndTime(statesOut, tOut));


        // Insert IC into list of state vectors
        statesAndTimes.resize(statesOut.size());


        for (int i = 0; i < statesAndTimes.size(); i++) {
            statesAndTimes[i].resize(ic.size() + 1);
            for (int j = 0; j < statesAndTimes[i].size(); j++) {
                if (j == 0) {
                    statesAndTimes[i][j] = tOut[i];
                } else {
                    statesAndTimes[i][j] = statesOut[i][j - 1];
                }
            }
        }

        return statesAndTimes;
    }


    template<class T>
    boost::python::list toPythonList(std::vector<T> vector) {
        typename std::vector<T>::iterator iter;
        boost::python::list list;
        for (iter = vector.begin(); iter != vector.end(); ++iter) {
            list.append(*iter);
        }
        return list;
    }

    template<class T>
    boost::python::list toTwoDimPythonList(std::vector<std::vector<T> > vector) {
        typename std::vector<std::vector<T> >::iterator iter;

        boost::python::list list;
        for (iter = vector.begin(); iter != vector.end(); ++iter) {
            list.append(toPythonList(*iter));
        }
        return list;
    }
};


BOOST_PYTHON_MODULE (boost_cr3bp) {
    class_<CR3BP>("CR3BP")
            .def("prop", &CR3BP::propPy);
};

int main() {
    cout << "main" << endl;

    double mu = 0.0122;

    vector<double> IC(6, 0);
    vector<double> t(2, 0);

    double stateDim = 6;
    double tDim = 2;
    std::vector<vector<double >> outTest;

    CR3BP cr3bp;

    // Initilal state and time span of integration
    IC[0] = 0.788;
    IC[1] = 0.200;
    IC[2] = 0.0;
    IC[3] = -0.88;
    IC[4] = 0.20;
    IC[5] = 0.0;
    t[0] = 0;
    t[1] = 0.5;

    // Integrate initial state
    outTest = cr3bp.prop(IC, t, mu, 6, 2, 1e-12, 1e-5);

    // Print time and states to a Console Window
    for (int i = 0; i < outTest.size(); i++) {
        std::cout << outTest[i][0] << " " << outTest[i][1] << " " << outTest[i][2] << " " << outTest[i][3] << " "
                  << outTest[i][4] << " " << outTest[i][5] << " " << outTest[i][6] << "\n";
    }


}