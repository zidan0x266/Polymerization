// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Chemical state of a lattice unit. The legacy code enumerated ten "montype" integers; they are
// really the tuple (species, number of radicals, number of unreacted vinyl groups):
//
//   legacy 0: monovinyl, r=0, v=1   (unreacted monomer)      legacy 5: divinyl, r=0, v=1 (pendant)
//   legacy 1: monovinyl, r=1, v=0   (propagating radical)    legacy 6: divinyl, r=1, v=0
//   legacy 2: monovinyl, r=0, v=0   (saturated / dead)       legacy 7: divinyl, r=2, v=0
//   legacy 3: divinyl,   r=0, v=2   (unreacted cross-linker) legacy 8: divinyl, r=0, v=0 (fully reacted)
//   legacy 4: divinyl,   r=1, v=1                            legacy 9: dormant initiator
//
// Keeping the tuple makes the reaction rules generic (see reaction/scheme.hpp) while the legacy
// code is still available for output compatibility.
#pragma once

#include <cstdint>
#include <string_view>

namespace bfm {

enum class Species : std::uint8_t { Monovinyl = 0, Divinyl = 1, Initiator = 2 };

[[nodiscard]] constexpr int max_vinyls(Species s) noexcept {
  switch (s) {
    case Species::Monovinyl: return 1;
    case Species::Divinyl: return 2;
    case Species::Initiator: return 0;
  }
  return 0;
}

struct UnitState {
  Species species{Species::Monovinyl};
  std::uint8_t radicals{0};
  std::uint8_t vinyls{0};

  constexpr bool operator==(const UnitState&) const noexcept = default;

  [[nodiscard]] constexpr bool has_radical() const noexcept { return radicals > 0; }
  [[nodiscard]] constexpr bool has_vinyl() const noexcept { return vinyls > 0; }
  /// A monomer that has not reacted yet: its vinyl groups are all "first" vinyls.
  [[nodiscard]] constexpr bool is_unreacted_monomer() const noexcept {
    return species != Species::Initiator && radicals == 0 && vinyls == max_vinyls(species);
  }
  /// A vinyl group left on a unit that has already reacted once (a cross-linker pendant).
  [[nodiscard]] constexpr bool has_pendant_vinyl() const noexcept {
    return has_vinyl() && !is_unreacted_monomer();
  }

  static constexpr UnitState monovinyl_monomer() noexcept { return {Species::Monovinyl, 0, 1}; }
  static constexpr UnitState divinyl_monomer() noexcept { return {Species::Divinyl, 0, 2}; }
  static constexpr UnitState dormant_initiator() noexcept { return {Species::Initiator, 0, 0}; }
  /// An initiator fragment that already carries its radical (ATRP: all chains are live
  /// from the start; FRP: the product of an initiation event).
  static constexpr UnitState primary_radical() noexcept { return {Species::Initiator, 1, 0}; }
};

/// The legacy "montype" code of a unit (0..9), used for structural-unit statistics and for the
/// residue names in PDB/PSF output. `bonds` distinguishes a dormant initiator (9) from an
/// initiator fragment that has been consumed into a chain (legacy reported it as 2).
[[nodiscard]] constexpr int legacy_type_code(UnitState s, int bonds) noexcept {
  switch (s.species) {
    case Species::Monovinyl:
      if (s.radicals == 0 && s.vinyls == 1) return 0;
      if (s.radicals == 1 && s.vinyls == 0) return 1;
      return 2;
    case Species::Divinyl:
      if (s.radicals == 0 && s.vinyls == 2) return 3;
      if (s.radicals == 1 && s.vinyls == 1) return 4;
      if (s.radicals == 0 && s.vinyls == 1) return 5;
      if (s.radicals == 1 && s.vinyls == 0) return 6;
      if (s.radicals == 2 && s.vinyls == 0) return 7;
      return 8;
    case Species::Initiator:
      if (s.radicals > 0) return 1;
      return bonds == 0 ? 9 : 2;
  }
  return -1;
}

inline constexpr int kNumLegacyTypeCodes = 10;

[[nodiscard]] constexpr std::string_view legacy_residue_name(int code) noexcept {
  constexpr std::string_view names[kNumLegacyTypeCodes] = {"MON", "ONE", "TWO", "THR", "FOU",
                                                            "FIV", "SIX", "SEV", "EIG", "INI"};
  return (code >= 0 && code < kNumLegacyTypeCodes) ? names[code] : "UNK";
}

/// The legacy `monreactive` column: 1 initiator / primary radical, 2 monovinyl, 4 divinyl.
[[nodiscard]] constexpr int legacy_species_code(Species s) noexcept {
  switch (s) {
    case Species::Initiator: return 1;
    case Species::Monovinyl: return 2;
    case Species::Divinyl: return 4;
  }
  return 0;
}

}  // namespace bfm
