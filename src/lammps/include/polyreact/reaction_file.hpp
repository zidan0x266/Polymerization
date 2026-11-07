// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Reaction files of fix poly/react. The format is that of ChemLab (J. Krajniak et al., J. Comput.
// Chem. 2018, 39, 1764; Z. Zhang et al., Adv. Theory Simul. 2019, 2, 1800102): an INI file with a
// [general] section, [group_<name>] sections for the bonds that reactions create, [ext_<name>]
// sections for extensions, and [reaction_<name>] sections, each holding one equation such as
//
//   FA(3,4) + MA(1,2) -> FA(1):DA(2)
//
// A bead is a type with an integer functionality; the left side gives a window [min, max) for
// the functionality and the right side the new type and the change of functionality. Additions to
// ChemLab: unimolecular reactions (A(0,1) -> B(1)), types without a window (any functionality, no
// change), bond_type/angle_type in groups and the competition rule in [general].
#pragma once

#include <climits>
#include <filesystem>
#include <istream>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace polyreact {

enum class Kind {
  Association,   ///< A + B -> A':B'     a bond forms between two beads
  Dissociation,  ///< A:B -> A' + B'     the bond between two bonded beads breaks
  Exchange,      ///< A:B + C -> A':C' + B'
  Unimolecular   ///< A -> A'            a bead changes on its own
};

/// A bead on the left side of an equation: functionality in [min, max).
struct Reactant {
  std::string type;
  int min{INT_MIN};
  int max{INT_MAX};
  [[nodiscard]] bool accepts(int functionality) const noexcept { return functionality >= min && functionality < max; }
  bool operator==(const Reactant& o) const noexcept { return type == o.type && min == o.min && max == o.max; }
};

/// A bead on the right side: its new type and the change of its functionality.
struct Product {
  std::string type;
  int delta{0};
  bool operator==(const Product& o) const noexcept { return type == o.type && delta == o.delta; }
};

/// A parsed equation; products are in the order of the reactants (for an exchange A, B, C).
struct Equation {
  Kind kind{Kind::Association};
  std::vector<Reactant> reactants;
  std::vector<Product> products;
};

/// How one of several reactions possible for the same pair is chosen.
enum class Competition {
  Independent,  ///< each is tried with its rate; one of those that pass is taken at random (ChemLab)
  Partition,    ///< one draw divides [0, 1) by the rates, rescaled if they sum to more than one
  Legacy        ///< the rule of the original 3DBFM code for two reactions (see engine.hpp)
};

struct Reaction {
  std::string name;  ///< section name without "reaction_"
  std::string equation_text;
  Equation equation;
  double rate{0.0};  ///< probability per check
  std::optional<double> diss_rate;  ///< dissociation: probability per check, if not given by rate
  double cutoff{0.0};
  double min_cutoff{0.0};
  bool intramolecular{false};  ///< allow a bond between beads of the same molecule
  bool intraresidual{false};   ///< allow a bond between beads of the same residue (initial molecule)
  bool virtual_bond{false};    ///< change the beads without creating a bond
  bool active{true};
  std::string group;
  std::vector<std::string> exclude_extensions;
  int line{0};  ///< line of the section header
};

struct Group {
  std::string name;
  std::optional<int> bond_type;   ///< LAMMPS bond type of the new bonds
  std::optional<int> angle_type;  ///< LAMMPS angle type of the angles the new bonds close
  std::vector<std::string> extensions;
  std::map<std::string, std::string> options;  ///< further ChemLab keys (potential, ...), kept but unused
};

struct Extension {
  std::string name;
  std::string type;
  std::map<std::string, std::string> options;
};

struct General {
  long interval{1};
  bool nearest{false};
  long max_per_interval{-1};
  Competition competition{Competition::Independent};
};

struct ReactionFile {
  General general;
  std::vector<Group> groups;
  std::vector<Reaction> reactions;
  std::vector<Extension> extensions;
  std::vector<std::string> notes;  ///< accepted ChemLab keys that fix poly/react does not use

  [[nodiscard]] const Group* find_group(std::string_view name) const noexcept;
  /// Bead types in order of first appearance in the reactions.
  [[nodiscard]] std::vector<std::string> type_names() const;
};

/// Parse one equation; throws std::invalid_argument with a readable message.
[[nodiscard]] Equation parse_equation(std::string_view text);

/// Parse a reaction file; throws std::runtime_error with "<source>:<line>: <message>".
[[nodiscard]] ReactionFile parse_reaction_file(std::istream& in, std::string_view source = "reaction file");
[[nodiscard]] ReactionFile read_reaction_file(const std::filesystem::path& path);

[[nodiscard]] std::string_view to_string(Competition c) noexcept;
[[nodiscard]] std::optional<Competition> parse_competition(std::string_view s) noexcept;

}  // namespace polyreact
