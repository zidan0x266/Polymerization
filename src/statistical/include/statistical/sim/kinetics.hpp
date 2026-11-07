// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Coordinate-free kinetic Monte Carlo (Gillespie's direct method) of radical (co)polymerization,
// as in the statistical models of Gao et al. (Macromolecules 2009) and Lyu et al. (Macromolecules
// 2018): every unit can meet every other one. With N units in the system, each pair of a radical
// and a monomer, cross-linker or pendant vinyl reacts at the rate k_p1/N, k_px/N or k_p2/N, each
// pair of radicals recombines at k_term/N, and a dormant initiator decomposes at k_init. An event
// is chosen with probability proportional to these rates; its radical, and then its partner, are
// drawn at random. Without cyclization, a partner from the radical's own molecule is drawn again.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <string>
#include <vector>

#include "poly/chem/progress.hpp"
#include "poly/chem/scheme.hpp"
#include "poly/random.hpp"
#include "statistical/model/network.hpp"

namespace statistical {

/// Outcome of one event.
enum class Event : std::uint8_t {
  Reaction,  ///< a reaction took place
  Idle,      ///< the drawn radical had no allowed partner: only the time advanced
  Exhausted  ///< no reaction is possible any more
};

class Kinetics {
 public:
  /// `cyclization = false` forbids reactions between units of the same molecule.
  Kinetics(Network& network, poly::chem::RadicalChainScheme::Features features, poly::chem::RateConstants rates, bool cyclization);

  /// Draw and carry out the next event, advancing the time.
  Event step(poly::Rng& rng, poly::chem::ReactionProgress& progress);

  [[nodiscard]] double time() const noexcept { return time_; }
  [[nodiscard]] std::size_t radicals() const noexcept { return radicals_.size(); }
  [[nodiscard]] std::size_t monomers() const noexcept { return monomers_.size(); }
  [[nodiscard]] std::size_t crosslinkers() const noexcept { return crosslinkers_.size(); }
  [[nodiscard]] std::size_t pendants() const noexcept { return pendants_.size(); }
  [[nodiscard]] std::size_t dormant() const noexcept { return dormant_.size(); }

  /// Check that the lists of reactive groups agree with the unit states.
  [[nodiscard]] std::expected<void, std::string> validate() const;

 private:
  Event propagate(std::vector<Unit>& targets, poly::Rng& rng, poly::chem::ReactionProgress& progress);
  Event propagate_pendant(poly::Rng& rng, poly::chem::ReactionProgress& progress);
  Event terminate(poly::Rng& rng, poly::chem::ReactionProgress& progress);
  Event initiate(poly::Rng& rng, poly::chem::ReactionProgress& progress);
  /// Whether `a` and `b` may react: distinct units, not bonded, and in different molecules
  /// when cyclization is off.
  [[nodiscard]] bool allowed(Unit a, Unit b) const noexcept;
  void react(Unit a, Unit b, const poly::chem::Reaction& r, poly::chem::ReactionProgress& progress);

  Network* network_;
  poly::chem::RadicalChainScheme::Features features_;
  poly::chem::RateConstants rates_;
  bool cyclization_;
  double volume_;
  double time_{0.0};
  // Reactive groups. A unit appears once per radical in radicals_ and once per pendant vinyl in
  // pendants_; uniform draws from these lists are uniform over radicals and vinyl groups.
  std::vector<Unit> radicals_;
  std::vector<Unit> monomers_;      ///< unreacted monovinyl monomers
  std::vector<Unit> crosslinkers_;  ///< unreacted divinyl cross-linkers
  std::vector<Unit> pendants_;
  std::vector<Unit> dormant_;       ///< dormant initiators
};

}  // namespace statistical
