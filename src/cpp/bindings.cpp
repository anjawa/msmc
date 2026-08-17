#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/vector.h>
#include <map>
#include <stdexcept>

#include "mayer/pair_potential.h"
#include "mayer/hard_sphere.h"
#include "mayer/lennard_jones.h"
#include "mcmc/metropolis_hasting.h"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_core, m) {
    nb::class_<PairPotential>(m, "PairPotential")
        .def("fij", &PairPotential::fij, "rij"_a)
        .def("prob_fij", &PairPotential::prob_fij, "rij"_a)
        .def("compute_integrand", [](const PairPotential& self,
                nb::ndarray<const double, nb::ndim<2>, nb::c_contig> pos,
                int virial_no) {
            if (pos.shape(1) != 3)
                throw std::invalid_argument("pos must have shape (N, 3)");
            return self.compute_integrand(pos.data(), virial_no);
        }, "pos"_a, "virial_no"_a);

    nb::class_<HardSphere, PairPotential>(m, "HardSphere")
        .def(nb::init<double>(), "sigma"_a);

    nb::class_<LennardJones, PairPotential>(m, "LennardJones")
        .def(nb::init<double, double, double, double>(),
             "sigma"_a, "epsilon"_a, "k"_a, "T"_a)
        .def("uij", &LennardJones::uij, "rij"_a);

    nb::class_<Model>(m, "Model")
        .def(nb::init<std::shared_ptr<PairPotential>, std::map<int, double>>(),
             "potential"_a, "virial"_a)
        .def_ro("virial", &Model::virial);

    nb::class_<MetropolisHasting>(m, "MetropolisHasting")
        .def(nb::init<Model, Model>(), "target_model"_a, "ref_model"_a)
        .def("sample_virial", &MetropolisHasting::sample_virial,
             "virial_no"_a, "num_samples"_a, "warmup"_a, "proposal_sigma"_a, "seed"_a)
        .def("sample_virial_overlap", &MetropolisHasting::sample_virial_overlap,
             "virial_no"_a, "num_samples"_a, "warmup"_a, "seed"_a);
}
