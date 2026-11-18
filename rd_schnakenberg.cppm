/*!
 * A 2D Reaction Diffusion system deriving from sm::rd_base.
 *
 * Author: Seb James
 */
module;

#include <vector>
#include <array>
#include <sstream>

export module rd_schnakenberg;

import rd_base;
import sm.hdfdata;

export namespace sm
{
    /*!
     * Two component Schnakenberg Reaction Diffusion system
     */
    template <class F>
    class rd_schnakenberg : public sm::rd_base<F>
    {
    public:
        /*!
         * Reactant A
         */
        alignas(alignof(std::vector<F>))
        std::vector<F> A;

        /*!
         * Reactant B
         */
        alignas(alignof(std::vector<F>))
        std::vector<F> B;

        /*
         * Schnakenberg parameters
         * F = k1 - k2 A + k3 A^2 B
         * G = k4        - k3 A^2 B
         */
        alignas(F) F k1 = F{1};
        alignas(F) F k2 = F{1};
        alignas(F) F k3 = F{1};
        alignas(F) F k4 = F{1};

        /*
         * The diffusion parameters.
         */
        alignas(F) F D_A = F{0.1};
        alignas(F) F D_B = F{0.1};

        /*!
         * Simple constructor; no arguments. Simply call rd_base constructor.
         */
        rd_schnakenberg() : sm::rd_base<F>() {}

        /*!
         * Destructor
         */
        ~rd_schnakenberg() {} // No operation. rd_base destructor will free hexgrid.

        /*!
         * Perform memory allocations, vector resizes and so on.
         */
        void allocate()
        {
            // Always call allocate() from the base class first.
            sm::rd_base<F>::allocate();
            // Resize and zero-initialise the various containers. Note that the size of a
            // 'vector variable' is given by the number of hexes in the hex grid which is
            // a member of this class (via its parent, rd_base)
            this->resize_vector_variable (this->A);
            this->resize_vector_variable (this->B);
        }

        /*!
         * Initialise variables and parameters. Carry out one-time computations required
         * of the model.
         */
        void init()
        {
            // Initialise A, B with noise
            this->noiseify_vector_variable (this->A, F{0.5}, 1);
            this->noiseify_vector_variable (this->B, F{0.6}, 1);
        }

        /*!
         * Save the variables to HDF5.
         */
        void save()
        {
            std::stringstream fname;
            fname << this->logpath << "/dat_";
            fname.width(5);
            fname.fill('0');
            fname << this->stepCount << ".h5";
            sm::hdfdata data(fname.str(), std::ios::out | std::ios::trunc);
            std::stringstream path;
            // The A variables
            path << "/A";
            data.add_contained_vals (path.str().c_str(), this->A);
            // The B variable
            path.str("");
            path.clear();
            path << "/B";
            data.add_contained_vals (path.str().c_str(), this->B);
        }

        /*!
         * Schnakenberg computation for reagent A
         */
        void compute_dAdt (std::vector<F>& A_, std::vector<F>& dAdt)
        {
            std::vector<F> lapA(this->nhex, F{0});
            this->compute_laplace (A_, lapA);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                dAdt[h] = this->k1 - (this->k2 * A_[h])
                + (this->k3 * A_[h] * A_[h] * this->B[h]) + this->D_A * lapA[h];
            }
        }

        /*!
         * Schnakenberg computation for reagent B
         */
        void compute_dBdt (std::vector<F>& B_, std::vector<F>& dBdt)
        {
            std::vector<F> lapB(this->nhex, F{0});
            this->compute_laplace (B_, lapB);
#pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                // G = k4        - k3 A^2 B
                dBdt[h] = this->k4 - (this->k3 * this->A[h] * this->A[h] * B_[h]) + this->D_B * lapB[h];
            }
        }

        /*!
         * Simulate one timestep of the model
         */
        void step()
        {
            this->stepCount++;

            // 1. 4th order Runge-Kutta computation for A
            {
                // Atst: "A at a test point". Atst is a temporary estimate for A.
                std::vector<F> Atst(this->nhex, F{0});
                std::vector<F> dAdt(this->nhex, F{0});
                std::vector<F> K1(this->nhex, F{0});
                std::vector<F> K2(this->nhex, F{0});
                std::vector<F> K3(this->nhex, F{0});
                std::vector<F> K4(this->nhex, F{0});

                /*
                 * Stage 1
                 */
                this->compute_dAdt (this->A, dAdt);
#pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; ++h) {
                    K1[h] = dAdt[h] * this->dt;
                    Atst[h] = this->A[h] + K1[h] * F{0.5} ;
                }

                /*
                 * Stage 2
                 */
                this->compute_dAdt (Atst, dAdt);
#pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; ++h) {
                    K2[h] = dAdt[h] * this->dt;
                    Atst[h] = this->A[h] + K2[h] * F{0.5};
                }

                /*
                 * Stage 3
                 */
                this->compute_dAdt (Atst, dAdt);
#pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; ++h) {
                    K3[h] = dAdt[h] * this->dt;
                    Atst[h] = this->A[h] + K3[h];
                }

                /*
                 * Stage 4
                 */
                this->compute_dAdt (Atst, dAdt);
#pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; ++h) {
                    K4[h] = dAdt[h] * this->dt;
                }

                /*
                 * Final sum together. This could be incorporated in the for loop for
                 * Stage 4, but I've separated it out for pedagogy.
                 */
#pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; ++h) {
                    A[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(F)6.0);
                }
            }

            // 2. 4th order Runge-Kutta computation of B
            {
                // Btst: "B at a test point". Btst is a temporary estimate for B.
                std::vector<F> Btst(this->nhex, F{0});
                std::vector<F> dBdt(this->nhex, F{0});
                std::vector<F> K1(this->nhex, F{0});
                std::vector<F> K2(this->nhex, F{0});
                std::vector<F> K3(this->nhex, F{0});
                std::vector<F> K4(this->nhex, F{0});

                /*
                 * Stage 1
                 */
                this->compute_dBdt (this->B, dBdt);
#pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; ++h) {
                    K1[h] = dBdt[h] * this->dt;
                    Btst[h] = this->B[h] + K1[h] * F{0.5} ;
                }

                /*
                 * Stage 2
                 */
                this->compute_dBdt (Btst, dBdt);
#pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; ++h) {
                    K2[h] = dBdt[h] * this->dt;
                    Btst[h] = this->B[h] + K2[h] * F{0.5};
                }

                /*
                 * Stage 3
                 */
                this->compute_dBdt (Btst, dBdt);
#pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; ++h) {
                    K3[h] = dBdt[h] * this->dt;
                    Btst[h] = this->B[h] + K3[h];
                }

                /*
                 * Stage 4
                 */
                this->compute_dBdt (Btst, dBdt);
#pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; ++h) {
                    K4[h] = dBdt[h] * this->dt;
                }

                /*
                 * Final sum together. This could be incorporated in the for loop for
                 * Stage 4, but I've separated it out for pedagogy.
                 */
#pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; ++h) {
                    B[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(F)6.0);
                }
            }
        }

    }; // rd_schnakenberg

} // namespace
