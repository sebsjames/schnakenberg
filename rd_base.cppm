module;

#include <cstdint>
#include <memory>
#include <sstream>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>
#include <filesystem>

export module rd_base;

import sm.mathconst;
import sm.random;
import sm.hexgrid;
import sm.hdfdata;
import sm.hexgrid.hdf;
import mplot.readcurves;

export namespace sm
{
    /*!
     * Base class for RD systems
     */
    template <typename F>
    class rd_base
    {
    public:

        /*!
         * Hex to hex d for the grid. Make smaller to increase the number
         * of Hexes being computed.
         */
        alignas(float) float hextohex_d = 0.01f;

        /*!
         * The 'span' of the initial hex of hexes that is created, onto which the boundary is then drawn.
         */
        alignas(float) float hexspan = 4.0f;

        /*!
         * Holds the number of hexes in the populated hexgrid
         */
        alignas(F) std::uint32_t nhex = 0;

        /*!
         * Over what length scale should some values fall off to zero
         * towards the boundary? Used in a couple of different locations.
         */
        alignas(F) F boundaryFalloffDist = F{0.02}; // 0.02 default

    protected:
        /*!
         * Our choice of dt.
         */
        alignas(F) F dt = F{0.00001};

        /*!
         * Compute half and sixth dt in constructor.
         */
        alignas(F) F halfdt = F{0};
        alignas(F) F sixthdt = F{0};

        /*!
         * Hex to hex distance. Populate this from hg.d after hg has been
         * initialised.
         */
        alignas(F) F d = F{1};
        alignas(F) F v = F{1};

        /*!
         * Parameters that depend on d and v:
         */
        alignas(F) F oneoverd = F{1}/this->d;
        alignas(F) F oneover2d = F{1}/(this->d+this->d);
        alignas(F) F oneover3d = F{1}/(F{3}*this->d);
        alignas(F) F oneover3dd = F{1} / F{3}*this->d*this->d;
        alignas(F) F twoover3dd = F{2} / F{3}*this->d*this->d;

        alignas(F) F oneoverv = F{1}/this->v;
        alignas(F) F twov = this->v+this->v;
        alignas(F) F oneover2v = F{1}/this->twov;
        alignas(F) F oneover4v = F{1}/(this->twov+this->twov);

    public:

        /*!
         * Track the number of computational steps that we've carried
         * out. Only to show a message saying "100 steps done...", but
         * that's reason enough.
         */
        alignas(F) std::uint32_t stepCount = 0;

        /*!
         * ALIGNAS REGION ENDS.
         *
         * Below here, there's no need to worry about alignas keywords.
         */

        /*!
         * Hold on to the ReadCurves object, so that the additional contours are available.
         */
        mplot::ReadCurves r;

        /*!
         * The hexgrid "background" for the Reaction Diffusion system.
         */
        std::unique_ptr<sm::hexgrid<F>> hg;

        /*!
         * The logpath for this model. Used when saving data out.
         */
        std::string logpath = "logs";

        /*!
         * Setter which attempts to ensure the path exists.
         */
        void setLogpath (const std::string p) {
            this->logpath = p;
            // Ensure log directory exists
            std::filesystem::create_directories (this->logpath);
        }

        /*!
         * Make the svgpath something that can be set by client code. If empty, then
         * set an elliptical boundary
         */
        std::string svgpath = "./trial.svg";

        /*!
         * Ellipse a and b values. Used to specify a boundary only if svgpath is empty.
         */
        float ellipse_a = 1.0f;
        float ellipse_b = 1.0f;

        /*!
         * Simple constructor; no arguments.
         */
        rd_base()
        {
            this->halfdt = this->dt/F{2};
            this->sixthdt = this->dt/F{6};
        }

        /*!
         * Destructor has nothing to do
         */
        ~rd_base() {}

        /*!
         * Utility functions to resize/zero vector-vectors that hold N
         * different RD variables.
         */
        void resize_vector_vector (std::vector<std::vector<F> >& vv, std::uint32_t N)
        {
            vv.resize (N);
            for (std::uint32_t i=0; i<N; ++i) {
                vv[i].resize (this->nhex, F{0});
            }
        }
        void zero_vector_vector (std::vector<std::vector<F> >& vv, std::uint32_t N)
        {
            for (std::uint32_t i=0; i<N; ++i) {
                vv[i].assign (this->nhex, F{0});
            }
        }
        void set_vector_vector (std::vector<std::vector<F> >& vv, std::uint32_t N, F value)
        {
            for (std::uint32_t i=0; i<N; ++i) {
                vv[i].assign (this->nhex, value);
            }
        }
        void resize_vector_vector (std::vector<std::vector<F> >& vv, std::uint32_t N, std::uint32_t M)
        {
            vv.resize (N);
            for (std::uint32_t i=0; i<N; ++i) {
                vv[i].resize (M, F{0});
            }
        }
        void zero_vector_vector (std::vector<std::vector<F> >& vv, std::uint32_t N, std::uint32_t M)
        {
            for (std::uint32_t i=0; i<N; ++i) {
                vv[i].assign (M, F{0});
            }
        }
        void set_vector_vector (std::vector<std::vector<F> >& vv, std::uint32_t N, std::uint32_t M, F value)
        {
            for (std::uint32_t i=0; i<N; ++i) {
                vv[i].assign (M, value);
            }
        }

        /*!
         * Resize/zero M vectors of N vectors of nhex Fs.
         */
        void resize_vector_vector_vector (std::vector<std::vector<std::vector<F> > >& vvv, std::uint32_t N, std::uint32_t M)
        {
            vvv.resize (M);
            for (std::uint32_t m=0; m<M; ++m) {
                vvv[m].resize (N);
                for (std::uint32_t i=0; i<N; ++i) {
                    vvv[m][i].resize (this->nhex, F{0});
                }
            }
        }
        void zero_vector_vector_vector (std::vector<std::vector<std::vector<F> > >& vvv, std::uint32_t N, std::uint32_t M)
        {
            for (std::uint32_t m=0; m<M; ++m) {
                for (std::uint32_t i=0; i<N; ++i) {
                    vvv[m][i].assign (M, F{0});
                }
            }
        }

        /*!
         * Resize/zero a variable that'll be nhex elements long
         */
        void resize_vector_variable (std::vector<F>& v) { v.resize (this->nhex, F{0}); }
        void zero_vector_variable (std::vector<F>& v) { v.assign (this->nhex, F{0}); }

        /*!
         * Resize/zero a parameter that'll be N elements long
         */
        void resize_vector_param (std::vector<F>& p, std::uint32_t N) { p.resize (N, F{0}); }
        void zero_vector_param (std::vector<F>& p, std::uint32_t N) { p.assign (N, F{0}); }

        /*!
         * Resize/zero a vector of M vectors of parameters that'll each be N
         * elements long
         */
        void resize_vector_vector_param (std::vector<std::vector<F> >& vp, std::uint32_t N, std::uint32_t M)
        {
            vp.resize (M);
            for (std::uint32_t m = 0; m<M; ++m) {
                vp[m].resize (N, F{0});
            }
        }
        void zero_vector_vector_param (std::vector<std::vector<F> >& vp, std::uint32_t N, std::uint32_t M)
        {
            for (std::uint32_t m = 0; m<M; ++m) {
                vp[m].assign (N, F{0});
            }
        }

        /*!
         * Resize/zero a gradient field
         */
        void resize_gradient_field (std::array<std::vector<F>, 2>& g)
        {
            g[0].resize (this->nhex, F{0});
            g[1].resize (this->nhex, F{0});
        }
        void zero_gradient_field (std::array<std::vector<F>, 2>& g)
        {
            g[0].assign (this->nhex, F{0});
            g[1].assign (this->nhex, F{0});
        }

        /*!
         * Resize/zero a vector size N containing arrays of two vector<F>s
         * which are the x and y components of a (mathematical) vector
         * field.
         */
        void resize_vector_array_vector (std::vector<std::array<std::vector<F>, 2> >& vav, std::uint32_t N)
        {
            vav.resize (N);
            for (std::uint32_t n = 0; n<N; ++n) {
                this->resize_gradient_field (vav[n]);
            }
        }
        void zero_vector_array_vector (std::vector<std::array<std::vector<F>, 2> >& vav, std::uint32_t N)
        {
            for (std::uint32_t i = 0; i<N; ++i) {
                this->zero_gradient_field (vav[i]);
            }
        }

        /*!
         * Resize/zero a vector M containing vectors of size N containing arrays of
         * two vector<F>s which are the x and y components of a (mathematical)
         * vector field.
         */
        void resize_vector_vector_array_vector (std::vector<std::vector<std::array<std::vector<F>, 2> > >& vvav, std::uint32_t N, std::uint32_t M)
        {
            vvav.resize (M);
            for (std::uint32_t j = 0; j<M; ++j) {
                vvav[j].resize (N);
                for (std::uint32_t n = 0; n<N; ++n) {
                    this->resize_gradient_field (vvav[j][n]);
                }
            }
        }
        void zero_vector_vector_array_vector (std::vector<std::vector<std::array<std::vector<F>, 2> > >& vvav, std::uint32_t N, std::uint32_t M)
        {
            for (std::uint32_t j = 0; j<M; ++j) {
                for (std::uint32_t i = 0; i<N; ++i) {
                    this->zero_gradient_field (vvav[j][i]);
                }
            }
        }

        /*!
         * Initialise a vector with noise, but with sigmoidal roll-off to
         * zero at the boundary.
         *
         * I apply a sigmoid to the boundary hexes, so that the noise
         * drops away towards the edge of the domain.
         */
        void noiseify_vector_variable (std::vector<F>& v, F offset, F gain)
        {
            sm::rand_uniform<F> rng;
            for (auto h : this->hg->hexen) {
                // boundarySigmoid. Jumps sharply (100, larger is
                // sharper) over length scale 0.05 to 1. So if
                // distance from boundary > 0.05, noise has normal
                // value. Close to boundary, noise is less.
                v[h.vi] = rng.get() * gain + offset;
                if (h.dist_to_boundary > -0.5) { // It's possible that dist_to_boundary is set to -1.0
                    F bSig = F{1} / ( F{1} + std::exp (-F{100}*(h.dist_to_boundary - this->boundaryFalloffDist)) );
                    v[h.vi] = v[h.vi] * bSig;
                }
            }
        }

        /*!
         * Perform memory allocations, vector resizes and so on.
         */
        virtual void allocate()
        {
            // Create a hexgrid. 3 is the 'x span' which determines how
            // many hexes are initially created. 0 is the z co-ordinate for the hexgrid.
            this->hg = std::make_unique<sm::hexgrid<F>>(this->hextohex_d, this->hexspan, 0);

            // Either set a boundary using the svgpath, or set it as an ellipse
            if (this->svgpath != "") {
                // Read the curves which make a boundary
                this->r.init (this->svgpath);
                // Set the boundary in the hexgrid
                auto rpath = this->r.getCorticalPath(); // always becurvepath<float, 3>
                this->hg->set_boundary (rpath);         // needs bezcurvepath<F, 3>
            } else {
                this->hg->set_elliptical_boundary (this->ellipse_a, this->ellipse_b);
            }
            // Compute the distances from the boundary
            this->hg->compute_distance_to_boundary();
            // Vector size comes from number of Hexes in the hexgrid
            this->nhex = this->hg->num();
            // Spatial d comes from the hexgrid, too.
            this->set_d(this->hg->get_d());
            this->set_v(this->hg->get_v());
        }

        /*!
         * Initialise variables and parameters. Carry out one-time
         * computations required of the model.
         */
        virtual void init() = 0;

    protected:
        /*!
         * Require private setters for d and v as there are several other
         * members that have to be updated at the same time.
         */
        virtual void set_d (F d_)
        {
            this->d = d_;
            this->oneoverd = 1.0/this->d;
            this->oneover2d = 1.0/(2*this->d);
            this->oneover3d = 1.0/(3*this->d);
            this->oneover3dd = 1.0 / (3*this->d*this->d);
            this->twoover3dd = 2.0 / (3*this->d*this->d);
        }

        virtual void set_v (F v_)
        {
            this->v = v_;
            this->oneoverv = 1.0/this->v;
            this->twov = this->v+this->v;
            this->oneover2v = 1.0/this->twov;
            this->oneover4v = 1.0/(this->twov+this->twov);
        }

    public:
        /*!
         * Public getters for d and v
         */
        F get_d() const { return this->d; }
        F get_v() const { return this->v; }

        /*!
         * Public accessors for dt
         */
        void set_dt (F _dt) {
            this->dt = _dt;
            this->halfdt = this->dt/2.0;
            this->sixthdt = this->dt/6.0;
        }
        F get_dt() const { return this->dt; }

    public:

        /*!
         * HDF5 file saving/loading methods.
         */
        /*!
         * Save a data frame
         */
        virtual void save() {}

        /*!
         * Save position information
         */
        void savePositions()
        {
            std::stringstream fname;
            fname << this->logpath << "/positions.h5";
            sm::hdfdata data(fname.str(), std::ios::out | std::ios::trunc);
            data.add_val ("/area", this->hg->num() * this->hg->get_hex_area());
            this->saveHexPositions (data);
        }

        /*!
         * Save positions of the hexes - note using two vector<float>s
         * that have been populated with the positions from the hexgrid,
         * to fit in with the HDF API.
         */
        void saveHexPositions (sm::hdfdata& dat)
        {
            dat.add_contained_vals ("/x", this->hg->d_x);
            dat.add_contained_vals ("/y", this->hg->d_y);

            // Add the neighbour information too.
            std::vector<float> x_ne = this->hg->d_x;
            std::vector<float> y_ne = this->hg->d_y;
            std::uint32_t count = 0;
            for (std::int32_t i : this->hg->d_n0) {
                if (i >= 0) {
                    x_ne[count] = this->hg->d_x[i];
                    y_ne[count] = this->hg->d_y[i];
                }
                ++count;
            }
            dat.add_contained_vals ("/x_ne", x_ne);
            dat.add_contained_vals ("/y_ne", y_ne);

            std::vector<float> x_nne = this->hg->d_x;
            std::vector<float> y_nne = this->hg->d_y;
            count = 0;
            for (std::int32_t i : this->hg->d_n1) {
                if (i >= 0) {
                    x_nne[count] = this->hg->d_x[i];
                    y_nne[count] = this->hg->d_y[i];
                }
                ++count;
            }
            dat.add_contained_vals ("/x_nne", x_nne);
            dat.add_contained_vals ("/y_nne", y_nne);

            std::vector<float> x_nnw = this->hg->d_x;
            std::vector<float> y_nnw = this->hg->d_y;
            count = 0;
            for (std::int32_t i : this->hg->d_n2) {
                if (i >= 0) {
                    x_nnw[count] = this->hg->d_x[i];
                    y_nnw[count] = this->hg->d_y[i];
                }
                ++count;
            }
            dat.add_contained_vals ("/x_nnw", x_nnw);
            dat.add_contained_vals ("/y_nnw", y_nnw);

            std::vector<float> x_nw = this->hg->d_x;
            std::vector<float> y_nw = this->hg->d_y;
            count = 0;
            for (std::int32_t i : this->hg->d_n3) {
                if (i >= 0) {
                    x_nw[count] = this->hg->d_x[i];
                    y_nw[count] = this->hg->d_y[i];
                }
                ++count;
            }
            dat.add_contained_vals ("/x_nw", x_nw);
            dat.add_contained_vals ("/y_nw", y_nw);

            std::vector<float> x_nsw = this->hg->d_x;
            std::vector<float> y_nsw = this->hg->d_y;
            count = 0;
            for (std::int32_t i : this->hg->d_n4) {
                if (i >= 0) {
                    x_nsw[count] = this->hg->d_x[i];
                    y_nsw[count] = this->hg->d_y[i];
                }
                ++count;
            }
            dat.add_contained_vals ("/x_nsw", x_nsw);
            dat.add_contained_vals ("/y_nsw", y_nsw);

            std::vector<float> x_nse = this->hg->d_x;
            std::vector<float> y_nse = this->hg->d_y;
            count = 0;
            for (std::int32_t i : this->hg->d_n5) {
                if (i >= 0) {
                    x_nse[count] = this->hg->d_x[i];
                    y_nse[count] = this->hg->d_y[i];
                }
                ++count;
            }
            dat.add_contained_vals ("/x_nse", x_nse);
            dat.add_contained_vals ("/y_nse", y_nse);

            // And hex to hex distance:
            dat.add_val ("/d", this->d);
        }

        /*
         * Computation methods
         */

        /*!
         * Normalise the vector of Fs f.
         */
        void normalise (std::vector<F>& f)
        {
            F maxf = -1e7;
            F minf = +1e7;

            // Determines min and max
            for (auto val : f) {
                if (val>maxf) { maxf = val; }
                if (val<minf) { minf = val; }
            }
            F scalef = 1.0 /(maxf - minf);

            std::vector<std::vector<F> > norm_a;
            this->resize_vector_vector (norm_a);
            for (std::uint32_t fi = 0; fi < f.size(); ++fi) {
                f[fi] = fmin (fmax (((f[fi]) - minf) * scalef, F{0}), F{1});
            }
        }

        /*!
         * Do a single step through the model.
         */
        virtual void step() = 0;

        /*!
         * 2D spatial integration of the function f. Result placed in gradf.
         *
         * For each Hex, work out the gradient in x and y directions
         * using whatever neighbours can contribute to an estimate.
         */
        void spacegrad2D (std::vector<F>& f, std::array<std::vector<F>, 2>& gradf) {

            // Note - East is positive x; North is positive y.
#pragma omp parallel for schedule(static)
            for (std::uint32_t hi=0; hi<this->nhex; ++hi) {

                // Find x gradient
                if (this->hg->has_ne (hi) && this->hg->has_nw (hi)) {
                    gradf[0][hi] = (f[this->hg->ne (hi)] - f[this->hg->nw (hi)]) * oneover2d;
                } else if (this->hg->has_ne (hi)) {
                    gradf[0][hi] = (f[this->hg->ne (hi)] - f[hi]) * oneoverd;
                } else if (this->hg->has_nw (hi)) {
                    gradf[0][hi] = (f[hi] - f[this->hg->nw (hi)]) * oneoverd;
                } else {
                    // zero gradient in x direction as no neighbours in
                    // those directions? Or possibly use the average of
                    // the gradient between the nw,ne and sw,se neighbours
                    gradf[0][hi] = F{0};
                }

                // Find y gradient
                if (this->hg->has_nnw (hi) && this->hg->has_nne (hi) && this->hg->has_nsw (hi) && this->hg->has_nse (hi)) {
                    // Full complement. Compute the mean of the nse->nne and nsw->nnw gradients
                    gradf[1][hi] = ( (f[this->hg->nne (hi)] - f[this->hg->nse (hi)]) + (f[this->hg->nnw (hi)] - f[this->hg->nsw (hi)]) ) * oneover4v;
                } else if (this->hg->has_nnw (hi) && this->hg->has_nne (hi)) {
                    gradf[1][hi] = ( (f[this->hg->nne (hi)] + f[this->hg->nnw (hi)]) * F{0.5} - f[hi]) * oneoverv;
                } else if (this->hg->has_nsw (hi) && this->hg->has_nse (hi)) {
                    gradf[1][hi] = (f[hi] - (f[this->hg->nse (hi)] + f[this->hg->nsw (hi)]) * F{0.5}) * oneoverv;
                } else if (this->hg->has_nnw (hi) && this->hg->has_nsw (hi)) {
                    gradf[1][hi] = (f[this->hg->nnw (hi)] - f[this->hg->nsw (hi)]) * oneover2v;
                } else if (this->hg->has_nne (hi) && this->hg->has_nse (hi)) {
                    gradf[1][hi] = (f[this->hg->nne (hi)] - f[this->hg->nse (hi)]) * oneover2v;
                } else {
                    // Leave grady at 0
                    gradf[1][hi] = F{0};
                }
            }
        }

        /*!
         * Compute laplacian of scalar field F, with result placed in lapF.
         */
        virtual void compute_laplace (const std::vector<F>& _F, std::vector<F>& lapF) {

            F norm  = F{2} / (F{3.0} * this->d * this->d);

#pragma omp parallel for schedule(static)
            for (std::uint32_t hi=0; hi<this->nhex; ++hi) {

                // 1. The D Del^2 term

                // Compute the sum around the neighbours
                F thesum = F{-6} * _F[hi];
                if (this->hg->has_ne (hi)) {
                    thesum += _F[this->hg->ne (hi)];
                } else {
                    thesum += _F[hi]; // A ghost neighbour-east with same value as Hex_0
                }
                if (this->hg->has_nne (hi)) {
                    thesum += _F[this->hg->nne (hi)];
                } else {
                    thesum += _F[hi];
                }
                if (this->hg->has_nnw (hi)) {
                    thesum += _F[this->hg->nnw (hi)];
                } else {
                    thesum += _F[hi];
                }
                if (this->hg->has_nw (hi)) {
                    thesum += _F[this->hg->nw (hi)];
                } else {
                    thesum += _F[hi];
                }
                if (this->hg->has_nsw (hi)) {
                    thesum += _F[this->hg->nsw (hi)];
                } else {
                    thesum += _F[hi];
                }
                if (this->hg->has_nse (hi)) {
                    thesum += _F[this->hg->nse (hi)];
                } else {
                    thesum += _F[hi];
                }

                lapF[hi] = norm * thesum;
            }
        }

    }; // rd_base

} // namespace
