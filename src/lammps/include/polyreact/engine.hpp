// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The rule engine of fix poly/react: given a bead and a candidate partner (or a bead alone), it
// decides which reaction of a reaction file, if any, takes place. It knows nothing about LAMMPS;
// the fix finds the candidates and carries out the outcome. This header and reaction_file.hpp are
// C++17, since the fix includes them and is compiled with the standard of the LAMMPS it is
// loaded into.
//
// A reaction is possible for a pair if the types and functionalities match the equation in either
// order (the bead that samples the partner first), the distance lies in [min_cutoff, cutoff), and
// the bond conditions hold: a bond forms only between beads not yet bonded, and between beads of
// the same molecule or residue only if the reaction allows it; a dissociation needs a bond. Each
// possible reaction is tried with its rate; when several are possible for the same pair, the
// competition rule of the file chooses one:
//
//   independent  each is tried on its own; one of those that pass is taken at random (ChemLab)
//   partition    one draw divides [0, 1) by the rates, rescaled if they sum to more than one
//   legacy       two reactions a (listed first) and b, as in the original 3DBFM code: a takes
//                place with p = r_a, or r_a / (r_a + r_b) if r_a + r_b >= 1; otherwise b with
//                probability (1 - p) min(1, r_b + p); with more than two, as partition
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "polyreact/reaction_file.hpp"

namespace polyreact {

struct Bead {
  int type{0};
  int functionality{0};
  bool operator==(const Bead& o) const noexcept { return type == o.type && functionality == o.functionality; }
};

/// What the caller knows about a candidate pair.
struct PairContext {
  double distance{0.0};
  bool bonded{false};
  bool same_molecule{false};
  bool same_residue{false};
};

/// A reaction of the file with its types as indices.
struct Rule {
  std::size_t reaction{0};  ///< index into ReactionFile::reactions
  Kind kind{Kind::Association};
  std::array<int, 2> type{};
  std::array<int, 2> min{};
  std::array<int, 2> max{};
  std::array<int, 2> new_type{};
  std::array<int, 2> delta{};
  double rate{0.0};
  double cutoff{0.0};
  double min_cutoff{0.0};
  bool intramolecular{false};
  bool intraresidual{false};
  bool virtual_bond{false};
  std::size_t group{0};  ///< index into ReactionFile::groups
};

struct Outcome {
  std::size_t reaction{0};  ///< index into ReactionFile::reactions
  Kind kind{Kind::Association};
  bool swapped{false};  ///< the partner is the first bead of the equation
  bool creates_bond{false};
  bool breaks_bond{false};
  Bead center;   ///< the center bead after the reaction
  Bead partner;  ///< the partner after the reaction (unimolecular: unused)
};

class Engine {
 public:
  /// Most reactions that may compete for one pair of types.
  static constexpr std::size_t kMaxCandidates = 32;

  /// `types` gives the index (0 to num_types - 1) of every bead type named in the file, e.g. the
  /// LAMMPS atom type minus one; other indices are types that never react.
  Engine(const ReactionFile& file, const std::map<std::string, int, std::less<>>& types, int num_types);
  /// Index the types in the order in which they first appear in the file.
  explicit Engine(const ReactionFile& file);

  /// Decide the reaction between `center`, the bead that looks for a partner, and `partner`.
  template <class Rng>
  [[nodiscard]] std::optional<Outcome> pair(Bead center, Bead partner, const PairContext& context, Rng& rng) const;
  /// Decide a unimolecular reaction of `bead`.
  template <class Rng>
  [[nodiscard]] std::optional<Outcome> single(Bead bead, Rng& rng) const;

  [[nodiscard]] int num_types() const noexcept { return num_types_; }
  [[nodiscard]] const std::vector<std::string>& type_names() const noexcept { return names_; }
  /// Index of a type name, or -1.
  [[nodiscard]] int type_index(std::string_view name) const noexcept;
  [[nodiscard]] const std::vector<Rule>& rules() const noexcept { return rules_; }
  [[nodiscard]] Competition competition() const noexcept { return competition_; }
  /// Largest cutoff of the pair reactions (the radius in which partners are searched).
  [[nodiscard]] double max_cutoff() const noexcept { return max_cutoff_; }
  /// Whether beads of this type can take part in a pair reaction.
  [[nodiscard]] bool reactive(int type) const noexcept { return in_range(type) && reactive_[static_cast<std::size_t>(type)]; }
  /// Whether some pair reaction has beads of types a and b (in either order); a cheap test before
  /// the context of a pair is worked out.
  [[nodiscard]] bool may_react(int a, int b) const noexcept {
    return in_range(a) && in_range(b) &&
           !pair_table_[static_cast<std::size_t>(a) * static_cast<std::size_t>(num_types_) + static_cast<std::size_t>(b)].empty();
  }
  [[nodiscard]] bool has_unimolecular() const noexcept { return has_single_; }

 private:
  struct Candidate {
    std::uint32_t rule;
    bool swapped;
  };

  void build(const ReactionFile& file, const std::map<std::string, int, std::less<>>& types);
  [[nodiscard]] bool in_range(int type) const noexcept { return type >= 0 && type < num_types_; }
  [[nodiscard]] static bool matches(const Rule& r, int slot, Bead b) noexcept {
    return b.type == r.type[slot] && b.functionality >= r.min[slot] && b.functionality < r.max[slot];
  }
  [[nodiscard]] static Bead apply(const Rule& r, int slot, Bead b) noexcept { return {r.new_type[slot], b.functionality + r.delta[slot]}; }
  [[nodiscard]] static bool context_allows(const Rule& r, const PairContext& c) noexcept;
  template <class Rng>
  [[nodiscard]] std::optional<std::size_t> choose(const Candidate* candidates, std::size_t n, Rng& rng) const;

  int num_types_{0};
  std::vector<std::string> names_;
  std::vector<Rule> rules_;
  std::vector<std::vector<std::uint32_t>> pair_table_;    ///< [a * num_types + b]: pair rules for types a, b
  std::vector<std::vector<std::uint32_t>> single_table_;  ///< [a]: unimolecular rules
  std::vector<std::uint8_t> reactive_;
  Competition competition_{Competition::Independent};
  double max_cutoff_{0.0};
  bool has_single_{false};
};

// ---- implementation of the templates ------------------------------------------------------------

template <class Rng>
std::optional<std::size_t> Engine::choose(const Candidate* candidates, std::size_t n, Rng& rng) const {
  if (n == 0) return std::nullopt;
  const auto rate = [&](std::size_t i) { return rules_[candidates[i].rule].rate; };
  if (competition_ == Competition::Independent) {
    std::array<std::size_t, kMaxCandidates> passed{};
    std::size_t m = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (rng.uniform() < rate(i)) passed[m++] = i;
    }
    if (m == 0) return std::nullopt;
    if (m == 1) return passed[0];
    return passed[std::min(m - 1, static_cast<std::size_t>(rng.uniform() * static_cast<double>(m)))];
  }
  if (competition_ == Competition::Legacy && n <= 2) {
    if (n == 1) {
      if (rng.uniform() < rate(0)) return 0;
      return std::nullopt;
    }
    const double a = rate(0), b = rate(1);
    const double p = (a + b >= 1.0) ? a / (a + b) : a;
    if (rng.uniform() < p) return 0;
    const double r2 = rng.uniform();
    const double r3 = rng.uniform();
    if (r2 >= p && r3 < b + p) return 1;
    return std::nullopt;
  }
  double sum = 0.0;
  for (std::size_t i = 0; i < n; ++i) sum += rate(i);
  const double scale = sum > 1.0 ? 1.0 / sum : 1.0;
  const double u = rng.uniform();
  double cumulative = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    cumulative += rate(i) * scale;
    if (u < cumulative) return i;
  }
  return std::nullopt;
}

template <class Rng>
std::optional<Outcome> Engine::pair(Bead center, Bead partner, const PairContext& context, Rng& rng) const {
  if (!in_range(center.type) || !in_range(partner.type)) return std::nullopt;
  const auto& list = pair_table_[static_cast<std::size_t>(center.type) * static_cast<std::size_t>(num_types_) +
                                 static_cast<std::size_t>(partner.type)];
  std::array<Candidate, kMaxCandidates> candidates{};
  std::size_t n = 0;
  for (const std::uint32_t index : list) {
    const Rule& r = rules_[index];
    if (!context_allows(r, context)) continue;
    if (matches(r, 0, center) && matches(r, 1, partner)) {
      candidates[n++] = {index, false};
    } else if (matches(r, 0, partner) && matches(r, 1, center)) {
      candidates[n++] = {index, true};
    }
  }
  const auto chosen = choose(candidates.data(), n, rng);
  if (!chosen) return std::nullopt;
  const Candidate c = candidates[*chosen];
  const Rule& r = rules_[c.rule];
  Outcome out;
  out.reaction = r.reaction;
  out.kind = r.kind;
  out.swapped = c.swapped;
  out.creates_bond = r.kind == Kind::Association && !r.virtual_bond;
  out.breaks_bond = r.kind == Kind::Dissociation;
  out.center = apply(r, c.swapped ? 1 : 0, center);
  out.partner = apply(r, c.swapped ? 0 : 1, partner);
  return out;
}

template <class Rng>
std::optional<Outcome> Engine::single(Bead bead, Rng& rng) const {
  if (!in_range(bead.type)) return std::nullopt;
  std::array<Candidate, kMaxCandidates> candidates{};
  std::size_t n = 0;
  for (const std::uint32_t index : single_table_[static_cast<std::size_t>(bead.type)]) {
    if (matches(rules_[index], 0, bead)) candidates[n++] = {index, false};
  }
  const auto chosen = choose(candidates.data(), n, rng);
  if (!chosen) return std::nullopt;
  const Rule& r = rules_[candidates[*chosen].rule];
  Outcome out;
  out.reaction = r.reaction;
  out.kind = Kind::Unimolecular;
  out.center = apply(r, 0, bead);
  return out;
}

}  // namespace polyreact
