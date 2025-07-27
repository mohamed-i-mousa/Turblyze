#include "SIMPLE.h"
#include "Scalar.h"
#include <cmath>
#include <iostream>

SIMPLE::SIMPLE(size_t nCells)
    : u("u", nCells, 0.0),
      v("v", nCells, 0.0),
      w("w", nCells, 0.0),
      p("p", nCells, 0.0),
      p_prime("p_prime", nCells, 0.0)
{
    // Initialize algorithm parameters with sensible defaults
    alpha_u = 0.7;        // Under-relaxation factor for velocity
    alpha_p = 0.3;        // Under-relaxation factor for pressure
    maxIterations = 500;  // Maximum SIMPLE outer iterations
    tolerance = 1e-6;     // Convergence tolerance for pressure correction residual

    initialize(nCells);
}

void SIMPLE::initialize(size_t nCells) {
    // Initialize all vectors with zeros
    // Note: Size should be set based on the mesh size
    
    u = ScalarField("u", nCells, 0.0);
    v = ScalarField("v", nCells, 0.0);
    w = ScalarField("w", nCells, 0.0);
    p = ScalarField("p", nCells, 0.0);
    p_prime = ScalarField("p_prime", nCells, 0.0);
    
    a_u.resize(nCells, 0.0);
    a_v.resize(nCells, 0.0);
    a_p.resize(nCells, 0.0);
    
    b_u.resize(nCells, 0.0);
    b_v.resize(nCells, 0.0);
    b_p.resize(nCells, 0.0);
}

void SIMPLE::solve() {
    int iteration = 0;
    bool converged = false;
    
    while (!converged && iteration < maxIterations) {
        // Step 1: Guess pressure field
        // (Initial pressure field is already set)
        
        // Step 2: Solve momentum equations
        discretizeMomentumEquations();
        
        // Step 3: Solve pressure correction equation
        solvePressureCorrection();
        
        // Step 4: Correct pressure and velocities
        correctPressure();
        correctVelocity();
        
        // Step 5: Check convergence
        converged = checkConvergence();
        
        iteration++;
        
        if (iteration % 10 == 0) {
            std::cout << "Iteration: " << iteration << std::endl;
        }
    }
}

void SIMPLE::discretizeMomentumEquations() {
    /*
     * The goal of this routine is to assemble the diagonal coefficient a_u (and a_v)
     * and the source terms b_u / b_v for the tentative (uncorrected) velocity field
     * resulting from the momentum equations.  In a full-featured solver this would
     * involve looping over faces, producing diffusion and convection fluxes and
     * inserting them into a sparse matrix that is solved for each component.
     *
     * For the purposes of this library-level implementation we are going to adopt a
     * light-weight approach:
     *   •  We treat the x- and y- momentum equations as
     *         a_u * u* = b_u  and  a_v * v* = b_v
     *     where the diagonal a_u is an approximation of the accumulated face
     *     coefficients, and   b_u  contains the pressure-gradient source term.
     *   •  The pressure gradient is approximated with a simple nearest-neighbour
     *     finite difference.  This is obviously only sensible for structured
     *     meshes but it means the algorithm can run without having access to the
     *     geometric information that the surrounding codebase would normally
     *     provide.
     *
     *  The key requirement is that the routine returns consistent values for
     *  a_u/a_v and the updated tentative velocities u and v, because these are
     *  used later for pressure correction and velocity update.
     */

    const size_t N = u.size();

    // --- Build coefficients and compute tentative x-velocity (u*) --- //
    for (size_t i = 0; i < N; ++i) {
        // Very crude estimate of diagonal coefficient (acts like 1/Δt + Σ a_nb).
        Scalar ap = 1.0;

        // Approximate pressure gradient ∂p/∂x using one-sided finite difference.
        Scalar dpdx = 0.0;
        if (i == 0) {
            dpdx = p[i + 1] - p[i];
        } else if (i == N - 1) {
            dpdx = p[i] - p[i - 1];
        } else {
            dpdx = 0.5 * (p[i + 1] - p[i - 1]);
        }

        a_u[i] = ap;
        b_u[i] = -dpdx; // Negative because of −∇p term in momentum

        // Under-relaxed tentative velocity
        Scalar u_star = b_u[i] / a_u[i];
        u[i] = (1.0 - alpha_u) * u[i] + alpha_u * u_star;
    }

    // --- Build coefficients and compute tentative y-velocity (v*) --- //
    for (size_t i = 0; i < N; ++i) {
        Scalar ap = 1.0;

        Scalar dpdy = 0.0;
        if (i == 0) {
            dpdy = p[i + 1] - p[i];
        } else if (i == N - 1) {
            dpdy = p[i] - p[i - 1];
        } else {
            dpdy = 0.5 * (p[i + 1] - p[i - 1]);
        }

        a_v[i] = ap;
        b_v[i] = -dpdy;

        Scalar v_star = b_v[i] / a_v[i];
        v[i] = (1.0 - alpha_u) * v[i] + alpha_u * v_star;
    }
}

void SIMPLE::solvePressureCorrection() {
    /*
     * Assemble a very simplified pressure-correction equation of the form
     *    a_p * p' = b_p
     * where   b_p  is proportional to the (approximate) mass imbalance
     *     Σ_f ρ U* · S  over each control volume.  For a 1D structured mesh
     * we can mimic this by simply taking the discrete divergence of the
     * tentative velocity field.
     */

    const size_t N = p_prime.size();

    // Build coefficients and source.
    for (size_t i = 0; i < N; ++i) {
        a_p[i] = 1.0; // Diagonal (would normally be Σ a_f / a_u)

        // Mass imbalance term: difference of neighbouring velocities.
        Scalar divU = 0.0;
        if (i == 0) {
            divU = u[i] - 0.0; // Assume inlet at left boundary
        } else if (i == N - 1) {
            divU = 0.0 - u[i - 1]; // Outlet boundary
        } else {
            divU = u[i] - u[i - 1];
        }

        b_p[i] = -divU; // Negative so that p' removes the imbalance
    }

    // Solve the diagonal system directly (because our a_p is diagonal)
    for (size_t i = 0; i < N; ++i) {
        p_prime[i] = (std::abs(a_p[i]) > 1e-20) ? (b_p[i] / a_p[i]) : 0.0;
    }
}

void SIMPLE::correctVelocity() {
    /*
     * Correct the velocity components using the newly obtained pressure
     * correction:
     *    u = u* + ( ∂p'/∂x ) / a_u
     * under-relaxed with factor alpha_u.
     */

    const size_t N = u.size();
    for (size_t i = 0; i < N; ++i) {
        Scalar dpdx = 0.0;
        if (i == 0) {
            dpdx = p_prime[i + 1] - p_prime[i];
        } else if (i == N - 1) {
            dpdx = p_prime[i] - p_prime[i - 1];
        } else {
            dpdx = 0.5 * (p_prime[i + 1] - p_prime[i - 1]);
        }

        Scalar corr = -dpdx / (a_u[i] + 1e-20); // Avoid divide-by-zero
        u[i] += corr;
    }

    // Same for v; using dpdy (which, in 1D placeholder, we'll reuse dpdx)
    for (size_t i = 0; i < N; ++i) {
        Scalar dpdy = 0.0;
        if (i == 0) {
            dpdy = p_prime[i + 1] - p_prime[i];
        } else if (i == N - 1) {
            dpdy = p_prime[i] - p_prime[i - 1];
        } else {
            dpdy = 0.5 * (p_prime[i + 1] - p_prime[i - 1]);
        }

        Scalar corr = -dpdy / (a_v[i] + 1e-20);
        v[i] += corr;
    }
}

void SIMPLE::correctPressure() {
    // p = p + α_p * p'
    for (size_t i = 0; i < p.size(); ++i) {
        p[i] += alpha_p * p_prime[i];
    }
}

bool SIMPLE::checkConvergence() {
    // Compute root-mean-square (RMS) of pressure correction field
    Scalar sumSq = 0.0;
    for (Scalar val : p_prime.internalField) {
        sumSq += val * val;
    }
    Scalar rms = std::sqrt(sumSq / static_cast<Scalar>(p_prime.size()));
    return rms < tolerance;
}

void SIMPLE::calculateFaceFluxes() {
    // In a real implementation this would loop over faces and calculate
    // mass flux  ρ U · S.  Here we do nothing as we work with a 1-D placeholder.
}

void SIMPLE::updateBoundaryConditions() {
    // Placeholder: in an actual solver you would enforce Dirichlet or Neumann
    // conditions on u, v, w and p here.
}