/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef PROFILETYPES_HPP
#define PROFILETYPES_HPP

#include <tuple>

#include "Config.hpp"

#include <boost/any.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/variant.hpp>

#include "Anisotropic.hpp"
#include "TanhDiffuse.hpp"
#include "TanhMembrane.hpp"
#include "Uniform.hpp"
#include "Yukawa.hpp"

typedef boost::mpl::vector<Uniform, Yukawa, Anisotropic, TanhDiffuse, TanhMembrane>
profile_types;

typedef boost::make_variant_over<profile_types>::type Permittivity;

namespace profiles {
    class isUniform : public boost::static_visitor<bool>
    {
    public:
        bool operator()(const Uniform & /* arg */) const { return true; }
        bool operator()(const Yukawa & /* arg */) const { return false; }
        bool operator()(const Anisotropic & /* arg */) const { return false; }
        bool operator()(const TanhDiffuse & /* arg */) const { return false; }
        bool operator()(const TanhMembrane & /* arg */) const { return false; }
    };

    inline bool uniform(const Permittivity & arg) {
        return boost::apply_visitor(isUniform(), arg);
    }

    class epsilonValue : public boost::static_visitor<boost::any>
    {
    public:
        double operator()(const Uniform & arg) const { return arg.epsilon; }
        std::tuple<double, double> operator()(const Yukawa & arg) const {
            return std::make_tuple(arg.epsilon, arg.kappa);
        }
        bool operator()(const Anisotropic & /* arg */) const { return false; }
        bool operator()(const TanhDiffuse & /* arg */) const { return false; }
        bool operator()(const TanhMembrane & /* arg */) const { return false; }
    };

    inline double epsilon(const Permittivity & arg) {
        return boost::any_cast<double>(boost::apply_visitor(epsilonValue(), arg));
    }
} // namespace profiles

#endif // PROFILETYPES_HPP
