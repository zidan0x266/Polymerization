// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2026 Zidan Zhang <zhangzidan@gmail.com>
//
// polyreact_preset: write the reaction file of fix poly/react for radical copolymerization with
// the chemistry of the 3DBFM module.
//
//   polyreact_preset frp|atrp [k_init=..] [k_p1=..] [k_px=..] [k_p2=..] [k_term=..] [interval=..]
//                    [cutoff=..] [bond_type=..] [angle_type=..] [competition=..] [cyclization=..]

#include <charconv>
#include <exception>
#include <iostream>
#include <ostream>  // std::println overload for std::ostream
#include <print>
#include <stdexcept>
#include <string>
#include <string_view>

#include "polyreact/presets.hpp"

namespace {

constexpr std::string_view kUsage =
    "usage: polyreact_preset frp|atrp [key=value ...]\n\n"
    "Writes the reaction file of fix poly/react for radical copolymerization of monomers and\n"
    "cross-linkers (chemistry of the 3dbfm module) to standard output. Keys: k_init, k_p1, k_px,\n"
    "k_p2, k_term (probabilities per check), interval, cutoff, bond_type, angle_type,\n"
    "competition (legacy, partition, independent), cyclization (true or false).\n";

double number(std::string_view key, std::string_view value) {
  double x = 0;
  const auto [ptr, ec] = std::from_chars(value.data(), value.data() + value.size(), x);
  if (ec != std::errc{} || ptr != value.data() + value.size()) throw std::invalid_argument(std::string(key) + ": not a number");
  return x;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2 || std::string_view(argv[1]) == "-h" || std::string_view(argv[1]) == "--help") {
      std::print("{}", kUsage);
      return argc < 2 ? 1 : 0;
    }
    const auto mechanism = polyreact::parse_mechanism(argv[1]);
    if (!mechanism) throw std::invalid_argument("the mechanism must be frp or atrp");
    polyreact::RadicalRates k;
    polyreact::PresetOptions o;
    for (int i = 2; i < argc; ++i) {
      const std::string_view arg(argv[i]);
      const auto eq = arg.find('=');
      if (eq == std::string_view::npos) throw std::invalid_argument(std::string(arg) + ": expected key=value");
      const std::string_view key = arg.substr(0, eq), value = arg.substr(eq + 1);
      if (key == "k_init") k.k_init = number(key, value);
      else if (key == "k_p1") k.k_p1 = number(key, value);
      else if (key == "k_px") k.k_px = number(key, value);
      else if (key == "k_p2") k.k_p2 = number(key, value);
      else if (key == "k_term") k.k_term = number(key, value);
      else if (key == "interval") o.interval = static_cast<long>(number(key, value));
      else if (key == "cutoff") o.cutoff = number(key, value);
      else if (key == "bond_type") o.bond_type = static_cast<int>(number(key, value));
      else if (key == "angle_type") o.angle_type = static_cast<int>(number(key, value));
      else if (key == "competition") {
        const auto c = polyreact::parse_competition(value);
        if (!c) throw std::invalid_argument("competition must be legacy, partition or independent");
        o.competition = *c;
      } else if (key == "cyclization") {
        if (value != "true" && value != "false") throw std::invalid_argument("cyclization must be true or false");
        o.cyclization = value == "true";
      } else {
        throw std::invalid_argument(std::string("unknown key ") + std::string(key));
      }
    }
    std::print("{}", polyreact::radical_reaction_file(*mechanism, k, o));
    return 0;
  } catch (const std::exception& e) {
    std::println(std::cerr, "polyreact_preset: {}", e.what());
    return 1;
  }
}
