// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2026 Zidan Zhang <zhangzidan@gmail.com>

#include "polyreact/presets.hpp"

#include <cctype>
#include <format>
#include <iterator>
#include <stdexcept>

namespace polyreact {

namespace {

struct Propagation {
  std::string_view label, attacker, target;
};
struct Termination {
  std::string_view label, a, b;
};

// Reaction categories of the original code, in order.
constexpr std::array<Propagation, 20> kPropagations{{
    {"primary_radical+monomer", "PRI", "MON"},
    {"chain_radical+monomer", "ONE", "MON"},
    {"crosslinker_radical(4)+monomer", "FOU", "MON"},
    {"crosslinker_radical(6)+monomer", "SIX", "MON"},
    {"crosslinker_diradical(7)+monomer", "SEV", "MON"},
    {"primary_radical+crosslinker", "PRI", "THR"},
    {"chain_radical+crosslinker", "ONE", "THR"},
    {"primary_radical+pendant(4)", "PRI", "FOU"},
    {"chain_radical+pendant(4)", "ONE", "FOU"},
    {"primary_radical+pendant(5)", "PRI", "FIV"},
    {"chain_radical+pendant(5)", "ONE", "FIV"},
    {"crosslinker_radical(4)+crosslinker", "FOU", "THR"},
    {"crosslinker_radical(6)+crosslinker", "SIX", "THR"},
    {"crosslinker_diradical(7)+crosslinker", "SEV", "THR"},
    {"crosslinker_radical(4)+pendant(4)", "FOU", "FOU"},
    {"crosslinker_radical(4)+pendant(5)", "FOU", "FIV"},
    {"crosslinker_radical(6)+pendant(4)", "SIX", "FOU"},
    {"crosslinker_diradical(7)+pendant(4)", "SEV", "FOU"},
    {"crosslinker_radical(6)+pendant(5)", "SIX", "FIV"},
    {"crosslinker_diradical(7)+pendant(5)", "SEV", "FIV"},
}};

constexpr std::array<Termination, 15> kTerminations{{
    {"term:primary+primary", "PRI", "PRI"},
    {"term:primary+chain", "PRI", "ONE"},
    {"term:chain+chain", "ONE", "ONE"},
    {"term:primary+radical(4)", "PRI", "FOU"},
    {"term:chain+radical(4)", "ONE", "FOU"},
    {"term:primary+radical(6)", "PRI", "SIX"},
    {"term:chain+radical(6)", "ONE", "SIX"},
    {"term:primary+diradical(7)", "PRI", "SEV"},
    {"term:chain+diradical(7)", "ONE", "SEV"},
    {"term:radical(4)+radical(4)", "FOU", "FOU"},
    {"term:radical(4)+radical(6)", "FOU", "SIX"},
    {"term:radical(4)+diradical(7)", "FOU", "SEV"},
    {"term:radical(6)+radical(6)", "SIX", "SIX"},
    {"term:radical(6)+diradical(7)", "SIX", "SEV"},
    {"term:diradical(7)+diradical(7)", "SEV", "SEV"},
}};

/// The unit after it gives up one radical.
std::string_view without_radical(std::string_view type) {
  if (type == "PRI" || type == "ONE") return "TWO";
  if (type == "FOU") return "FIV";
  if (type == "SIX") return "EIG";
  if (type == "SEV") return "SIX";
  throw std::logic_error("not a radical unit");
}

/// The unit after a radical adds across one of its vinyl groups.
std::string_view attacked(std::string_view type) {
  if (type == "MON") return "ONE";
  if (type == "THR") return "FOU";
  if (type == "FOU") return "SEV";
  if (type == "FIV") return "SIX";
  throw std::logic_error("not a unit with a vinyl group");
}

}  // namespace

std::optional<Mechanism> parse_mechanism(std::string_view s) noexcept {
  std::string v(s);
  for (char& c : v) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  if (v == "frp") return Mechanism::Frp;
  if (v == "atrp") return Mechanism::Atrp;
  return std::nullopt;
}

std::string radical_reaction_file(Mechanism mechanism, const RadicalRates& k, const PresetOptions& o) {
  for (const double rate : {k.k_init, k.k_p1, k.k_p2, k.k_term, k.k_px.value_or(0.0)}) {
    if (!(rate >= 0 && rate <= 1)) throw std::invalid_argument("rate constants are probabilities per check in [0, 1]");
  }
  if (o.interval < 1 || !(o.cutoff > 0)) throw std::invalid_argument("interval and cutoff must be positive");
  const bool frp = mechanism == Mechanism::Frp;
  const double k_px = k.k_px.value_or(k.k_p1);
  std::string text;
  auto out = std::back_inserter(text);
  std::format_to(out,
                 "# {}\n"
                 "# of monovinyl monomers and divinyl cross-linkers, with the chemistry of the 3DBFM module.\n"
                 "# Bead types: MON monomer, THR cross-linker, PRI initiator fragment with its radical, ONE\n"
                 "# monovinyl radical, TWO saturated unit, FOU/FIV/SIX/SEV/EIG partially or fully reacted\n"
                 "# cross-linker, INI dormant initiator. Rates are probabilities per check:\n"
                 "# k_init={} k_p1={} k_px={} k_p2={} k_term={}\n\n"
                 "[general]\n"
                 "interval = {}\n"
                 "competition = {}\n\n"
                 "[group_radical]\n"
                 "bond_type = {}\n",
                 frp ? "Free radical polymerization (FRP)" : "Atom transfer radical polymerization (ATRP)", frp ? k.k_init : 0.0, k.k_p1, k_px,
                 k.k_p2, frp ? k.k_term : 0.0, o.interval, to_string(o.competition), o.bond_type);
  if (o.angle_type) std::format_to(out, "angle_type = {}\n", *o.angle_type);
  const auto reaction = [&](std::string_view label, std::string_view equation, double rate, bool pair) {
    std::format_to(out, "\n[reaction_{}]\nreaction = {}\nrate = {}\n", label, equation, rate);
    if (pair) std::format_to(out, "cutoff = {}\nintramolecular = {}\nintraresidual = true\n", o.cutoff, o.cyclization);
    std::format_to(out, "group = radical\n");
  };
  for (const Propagation& p : kPropagations) {
    const double rate = p.target == "MON" ? k.k_p1 : (p.target == "THR" ? k_px : k.k_p2);
    reaction(p.label, std::format("{} + {} -> {}:{}", p.attacker, p.target, without_radical(p.attacker), attacked(p.target)), rate, true);
  }
  if (frp) {
    for (const Termination& t : kTerminations) {
      reaction(t.label, std::format("{} + {} -> {}:{}", t.a, t.b, without_radical(t.a), without_radical(t.b)), k.k_term, true);
    }
    reaction("initiation", "INI -> PRI", k.k_init, false);
  }
  return text;
}

}  // namespace polyreact
