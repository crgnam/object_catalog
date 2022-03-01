#include <math.h>
#include <vector>
#include <iostream>

using Scalar = double;

// Function for converting orbital elements to states:
template <typename Scalar>
std::vector<Scalar> kepler_to_state(Scalar mu, Scalar a, Scalar e, Scalar i, Scalar peri, Scalar node, Scalar M0, Scalar tsince_epoch){

    // Calculate Mean Anomaly:
    Scalar n = std::sqrt(mu/(a*a*a));
    Scalar M = M0 + n*tsince_epoch;

    // Calculate eccentric anomaly:
    Scalar convergencePercentage = 0.05;
    Scalar relativeDifference;
    Scalar E = 1.0;
    Scalar f;
    Scalar df;
    Scalar E_new;
    for (int ii = 0; ii < 1000; ++ii) {
        // Set up Newton-Raphson method
        f = E - e*std::sin(E) - M;
        df = 1 - e*std::cos(E);
        E_new = E - f/df;

        // Check for convergence
        relativeDifference = abs(E_new - E) / E * 100;
        if (relativeDifference < convergencePercentage) {
            break;
        }
        E = E_new;
    }

    // Calculate the true anomaly:
    Scalar theta = std::atan2(std::sqrt(1- e*e)*std::sin(E), std::cos(E) - e);

    // Calculate the orbital momentum and radial position:
    Scalar h = std::sqrt(mu*a*(1-e*e));
    Scalar r_mag = ((h*h)/mu)*(1/(1+(e*std::cos(theta))));

    // Calculate the states in peri-focal coordinates:
    std::vector<Scalar> R_pqw = {r_mag*std::cos(theta), r_mag*std::sin(theta), 0.0};

    // Calculate vectorized rotations:
    Scalar a11 =  std::cos(node)*std::cos(peri) - std::sin(node)*std::sin(peri)*std::cos(i);
    Scalar a12 =  std::sin(node)*std::cos(peri) + std::cos(node)*std::sin(peri)*std::cos(i);
    Scalar a13 =  std::sin(peri)*std::sin(i);
    Scalar a21 = -std::cos(node)*std::sin(peri) - std::sin(node)*std::cos(peri)*std::cos(i);
    Scalar a22 = -std::sin(node)*std::sin(peri) + std::cos(node)*std::cos(peri)*std::cos(i);
    Scalar a23 =  std::cos(peri)*std::sin(i);
    Scalar a31 =  std::sin(node)*std::sin(i);
    Scalar a32 = -std::cos(node)*std::sin(i);
    Scalar a33 =  std::cos(i);

    // Apply rotations to obtain position in inertial frame:
    std::vector<Scalar> r = {a11*R_pqw[0] + a12*R_pqw[1], a21*R_pqw[0] + a22*R_pqw[1], a31*R_pqw[0] + a32*R_pqw[1]};

    return r;
}

// Main execution loop:
int main(){
    // Constants for conversion.  These will be removed as ephem table will already be converted:
    Scalar AU = 149597870.6907;
    Scalar deg2rad = 3.14159265/180.0;

    // Define placeholder values:
    Scalar mu = 132712440041.93938;
    Scalar e = 0.07850100198908602;
    Scalar a = 2.766043062222408*AU;
    Scalar i = 10.58769305845201*deg2rad;
    Scalar node = 80.26859547732911*deg2rad;
    Scalar peri = 73.63703979153577*deg2rad;
    Scalar M = 291.3755993017663*deg2rad;
    Scalar epoch = 695995269.1844931;

    // Calculate time since epoch:
    Scalar tsince_epoch = epoch;

    // Evaluate:
    std::vector<Scalar> r = kepler_to_state<Scalar>(mu, a, e, i, peri, node, M, tsince_epoch);

    // Display the result:
    std::cout << r[0] << ", " << r[1] << ", " << r[2] <<"\n";
}