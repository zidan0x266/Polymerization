// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2026 Zidan Zhang <zhangzidan@gmail.com>

#include "polyreact/reaction_file.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <format>
#include <fstream>
#include <set>
#include <stdexcept>
#include <utility>

namespace polyreact {

namespace {

std::string_view trim(std::string_view s) {
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.remove_suffix(1);
  return s;
}

std::string lower(std::string_view s) {
  std::string out(s);
  for (char& c : out) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return out;
}

std::vector<std::string_view> split(std::string_view s, char sep) {
  std::vector<std::string_view> parts;
  std::size_t start = 0;
  for (std::size_t i = 0; i <= s.size(); ++i) {
    if (i == s.size() || s[i] == sep) {
      parts.push_back(trim(s.substr(start, i - start)));
      start = i + 1;
    }
  }
  return parts;
}

bool valid_name(std::string_view s) {
  if (s.empty() || !(std::isalpha(static_cast<unsigned char>(s[0])) || s[0] == '_')) return false;
  return std::ranges::all_of(s, [](char c) { return std::isalnum(static_cast<unsigned char>(c)) || c == '_'; });
}

int parse_int(std::string_view s, std::string_view what) {
  s = trim(s);
  int value = 0;
  const auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
  if (ec != std::errc{} || ptr != s.data() + s.size()) throw std::invalid_argument(std::format("{} '{}' is not an integer", what, s));
  return value;
}

/// "NAME" or "NAME(a)" or "NAME(a, b)": the name and the numbers in parentheses.
std::pair<std::string, std::vector<int>> parse_bead(std::string_view token) {
  token = trim(token);
  const auto open = token.find('(');
  const std::string_view name = trim(token.substr(0, open));
  if (!valid_name(name)) throw std::invalid_argument(std::format("'{}' is not a bead type", token));
  std::vector<int> numbers;
  if (open != std::string_view::npos) {
    if (token.back() != ')') throw std::invalid_argument(std::format("'{}': missing ')'", token));
    for (const std::string_view n : split(token.substr(open + 1, token.size() - open - 2), ',')) numbers.push_back(parse_int(n, "functionality"));
  }
  return {std::string(name), numbers};
}

Reactant parse_reactant(std::string_view token) {
  auto [name, numbers] = parse_bead(token);
  if (numbers.empty()) return {std::move(name), INT_MIN, INT_MAX};
  if (numbers.size() != 2) throw std::invalid_argument(std::format("'{}': a reactant takes a window (min, max)", trim(token)));
  if (numbers[0] >= numbers[1]) throw std::invalid_argument(std::format("'{}': the window [min, max) is empty", trim(token)));
  return {std::move(name), numbers[0], numbers[1]};
}

Product parse_product(std::string_view token) {
  auto [name, numbers] = parse_bead(token);
  if (numbers.size() > 1) throw std::invalid_argument(std::format("'{}': a product takes one change of functionality", trim(token)));
  return {std::move(name), numbers.empty() ? 0 : numbers[0]};
}

bool parse_bool(std::string_view key, std::string_view value) {
  const std::string v = lower(trim(value));
  if (v == "1" || v == "true" || v == "yes" || v == "on") return true;
  if (v == "0" || v == "false" || v == "no" || v == "off") return false;
  throw std::invalid_argument(std::format("{}: '{}' is not a boolean", key, value));
}

double parse_double(std::string_view key, std::string_view value) {
  value = trim(value);
  double x = 0;
  const auto [ptr, ec] = std::from_chars(value.data(), value.data() + value.size(), x);
  if (ec != std::errc{} || ptr != value.data() + value.size()) throw std::invalid_argument(std::format("{}: '{}' is not a number", key, value));
  return x;
}

long parse_long(std::string_view key, std::string_view value) {
  value = trim(value);
  long x = 0;
  const auto [ptr, ec] = std::from_chars(value.data(), value.data() + value.size(), x);
  if (ec != std::errc{} || ptr != value.data() + value.size()) throw std::invalid_argument(std::format("{}: '{}' is not an integer", key, value));
  return x;
}

std::vector<std::string> parse_list(std::string_view value) {
  std::vector<std::string> out;
  for (const std::string_view item : split(value, ',')) {
    if (!item.empty()) out.emplace_back(item);
  }
  return out;
}

}  // namespace

Equation parse_equation(std::string_view text) {
  const auto arrow = text.find("->");
  if (arrow == std::string_view::npos || text.find("->", arrow + 2) != std::string_view::npos) {
    throw std::invalid_argument(std::format("'{}': an equation needs exactly one '->'", trim(text)));
  }
  const std::string_view left = trim(text.substr(0, arrow)), right = trim(text.substr(arrow + 2));
  const auto plus = split(left, '+');
  Equation eq;
  const auto expect = [&](bool ok, std::string_view form) {
    if (!ok) throw std::invalid_argument(std::format("'{}' does not have the form {}", trim(text), form));
  };
  if (plus.size() == 2 && plus[0].contains(':')) {
    eq.kind = Kind::Exchange;
    const auto ab = split(plus[0], ':');
    const auto rp = split(right, '+');
    expect(ab.size() == 2 && !plus[1].contains(':') && rp.size() == 2, "A:B + C -> A':C' + B'");
    const auto ac = split(rp[0], ':');
    expect(ac.size() == 2 && !rp[1].contains(':'), "A:B + C -> A':C' + B'");
    eq.reactants = {parse_reactant(ab[0]), parse_reactant(ab[1]), parse_reactant(plus[1])};
    eq.products = {parse_product(ac[0]), parse_product(rp[1]), parse_product(ac[1])};
  } else if (plus.size() == 2) {
    eq.kind = Kind::Association;
    const auto rp = split(right, ':');
    expect(!plus[1].contains(':') && rp.size() == 2 && !right.contains('+'), "A + B -> A':B'");
    eq.reactants = {parse_reactant(plus[0]), parse_reactant(plus[1])};
    eq.products = {parse_product(rp[0]), parse_product(rp[1])};
  } else if (plus.size() == 1 && left.contains(':')) {
    eq.kind = Kind::Dissociation;
    const auto ab = split(left, ':');
    const auto rp = split(right, '+');
    expect(ab.size() == 2 && rp.size() == 2 && !right.contains(':'), "A:B -> A' + B'");
    eq.reactants = {parse_reactant(ab[0]), parse_reactant(ab[1])};
    eq.products = {parse_product(rp[0]), parse_product(rp[1])};
  } else if (plus.size() == 1) {
    eq.kind = Kind::Unimolecular;
    expect(!right.contains('+') && !right.contains(':'), "A -> A'");
    eq.reactants = {parse_reactant(left)};
    eq.products = {parse_product(right)};
  } else {
    expect(false, "A + B -> A':B', A:B -> A' + B', A:B + C -> A':C' + B' or A -> A'");
  }
  return eq;
}

std::string_view to_string(Competition c) noexcept {
  switch (c) {
    case Competition::Independent: return "independent";
    case Competition::Partition: return "partition";
    case Competition::Legacy: return "legacy";
  }
  return "?";
}

std::optional<Competition> parse_competition(std::string_view s) noexcept {
  const std::string v = lower(trim(s));
  if (v == "independent") return Competition::Independent;
  if (v == "partition") return Competition::Partition;
  if (v == "legacy") return Competition::Legacy;
  return std::nullopt;
}

const Group* ReactionFile::find_group(std::string_view name) const noexcept {
  const auto it = std::ranges::find(groups, name, &Group::name);
  return it == groups.end() ? nullptr : &*it;
}

std::vector<std::string> ReactionFile::type_names() const {
  std::vector<std::string> names;
  const auto add = [&](const std::string& n) {
    if (std::ranges::find(names, n) == names.end()) names.push_back(n);
  };
  for (const Reaction& r : reactions) {
    for (const Reactant& x : r.equation.reactants) add(x.type);
    for (const Product& x : r.equation.products) add(x.type);
  }
  return names;
}

ReactionFile parse_reaction_file(std::istream& in, std::string_view source) {
  ReactionFile file;
  enum class Section { None, General, Group, Extension, Reaction } section = Section::None;
  std::set<std::string> seen;
  int line_number = 0;
  std::string line;
  const auto fail = [&](std::string_view message, int at) {
    throw std::runtime_error(std::format("{}:{}: {}", source, at, message));
  };
  const auto note = [&](std::string_view key, std::string_view where) {
    file.notes.push_back(std::format("line {}: '{}' in {} is accepted for ChemLab compatibility but not used", line_number, key, where));
  };

  while (std::getline(in, line)) {
    ++line_number;
    std::string_view s = trim(line);
    if (s.empty() || s.front() == '#' || s.front() == ';') continue;
    for (std::size_t i = 1; i < s.size(); ++i) {  // inline comment: " #..."
      if (s[i] == '#' && std::isspace(static_cast<unsigned char>(s[i - 1]))) {
        s = trim(s.substr(0, i));
        break;
      }
    }
    if (s.front() == '[') {
      if (s.back() != ']') fail("a section header must end with ']'", line_number);
      const std::string name(trim(s.substr(1, s.size() - 2)));
      if (!seen.insert(name).second) fail(std::format("section [{}] appears twice", name), line_number);
      if (name == "general") {
        section = Section::General;
      } else if (name.starts_with("group_")) {
        section = Section::Group;
        file.groups.emplace_back().name = name.substr(6);
      } else if (name.starts_with("ext_")) {
        section = Section::Extension;
        file.extensions.emplace_back().name = name.substr(4);
      } else if (name.starts_with("reaction_")) {
        section = Section::Reaction;
        Reaction& r = file.reactions.emplace_back();
        r.name = name.substr(9);
        r.line = line_number;
      } else {
        fail(std::format("unknown section [{}] (expected general, group_*, ext_* or reaction_*)", name), line_number);
      }
      continue;
    }
    // key: value or key = value, split at the first ':' or '=' (as Python's ConfigParser)
    const auto sep = s.find_first_of(":=");
    if (sep == std::string_view::npos || sep == 0) fail(std::format("expected 'key: value' or 'key = value', got '{}'", s), line_number);
    const std::string key = lower(trim(s.substr(0, sep)));
    const std::string_view value = trim(s.substr(sep + 1));
    try {
      switch (section) {
        case Section::None: fail("a key before the first section", line_number); break;
        case Section::General: {
          General& g = file.general;
          if (key == "interval") {
            g.interval = parse_long(key, value);
          } else if (key == "nearest") {
            g.nearest = parse_bool(key, value);
          } else if (key == "max_per_interval") {
            g.max_per_interval = parse_long(key, value);
          } else if (key == "competition") {
            const auto c = parse_competition(value);
            if (!c) fail("competition must be independent, partition or legacy", line_number);
            g.competition = *c;
          } else if (key == "pair_distances_filename" || key == "bond_limit") {
            note(key, "[general]");
          } else {
            fail(std::format("unknown key '{}' in [general]", key), line_number);
          }
          break;
        }
        case Section::Group: {
          Group& g = file.groups.back();
          if (key == "bond_type") {
            g.bond_type = static_cast<int>(parse_long(key, value));
          } else if (key == "angle_type") {
            g.angle_type = static_cast<int>(parse_long(key, value));
          } else if (key == "extensions") {
            g.extensions = parse_list(value);
          } else if (key == "potential" || key == "potential_options" || key == "connectivity_map" || key == "eq_length" ||
                     key == "final_type" || key == "alpha") {
            g.options[key] = std::string(value);
            note(key, std::format("[group_{}]", g.name));
          } else {
            fail(std::format("unknown key '{}' in [group_{}]", key, g.name), line_number);
          }
          break;
        }
        case Section::Extension: {
          Extension& e = file.extensions.back();
          if (key == "ext_type") {
            e.type = std::string(value);
          } else {
            e.options[key] = std::string(value);
          }
          break;
        }
        case Section::Reaction: {
          Reaction& r = file.reactions.back();
          if (key == "reaction") {
            r.equation_text = std::string(value);
            r.equation = parse_equation(value);
          } else if (key == "rate") {
            r.rate = parse_double(key, value);
          } else if (key == "diss_rate") {
            r.diss_rate = parse_double(key, value);
          } else if (key == "cutoff") {
            r.cutoff = parse_double(key, value);
          } else if (key == "min_cutoff") {
            r.min_cutoff = parse_double(key, value);
          } else if (key == "intramolecular") {
            r.intramolecular = parse_bool(key, value);
          } else if (key == "intraresidual") {
            r.intraresidual = parse_bool(key, value);
          } else if (key == "virtual") {
            r.virtual_bond = parse_bool(key, value);
          } else if (key == "active") {
            r.active = parse_bool(key, value);
          } else if (key == "group") {
            r.group = std::string(value);
          } else if (key == "exclude_extensions") {
            r.exclude_extensions = parse_list(value);
          } else if (key == "alpha" || key == "sigma" || key == "eq_distance") {
            note(key, std::format("[reaction_{}]", r.name));
          } else {
            fail(std::format("unknown key '{}' in [reaction_{}]", key, r.name), line_number);
          }
          break;
        }
      }
    } catch (const std::invalid_argument& e) {
      fail(e.what(), line_number);
    }
  }

  // Cross-checks once everything is read, so the order of the sections does not matter.
  for (const Extension& e : file.extensions) {
    if (e.type.empty()) fail(std::format("[ext_{}] has no ext_type", e.name), 0);
  }
  for (const Group& g : file.groups) {
    for (const std::string& x : g.extensions) {
      if (std::ranges::find(file.extensions, x, &Extension::name) == file.extensions.end()) {
        fail(std::format("[group_{}] uses the undefined extension '{}'", g.name, x), 0);
      }
    }
  }
  for (const Reaction& r : file.reactions) {
    const auto at = r.line;
    if (r.equation_text.empty()) fail(std::format("[reaction_{}] has no reaction", r.name), at);
    if (!(r.rate >= 0 && r.rate <= 1)) fail(std::format("[reaction_{}]: rate is a probability per check in [0, 1]", r.name), at);
    if (r.diss_rate && !(*r.diss_rate >= 0 && *r.diss_rate <= 1)) fail(std::format("[reaction_{}]: diss_rate must be in [0, 1]", r.name), at);
    if (r.diss_rate && r.equation.kind != Kind::Dissociation) fail(std::format("[reaction_{}]: diss_rate is only for dissociation", r.name), at);
    if (r.equation.kind != Kind::Unimolecular) {
      if (!(r.cutoff > 0)) fail(std::format("[reaction_{}] needs a positive cutoff", r.name), at);
      if (!(r.min_cutoff >= 0 && r.min_cutoff < r.cutoff)) fail(std::format("[reaction_{}]: min_cutoff must be in [0, cutoff)", r.name), at);
    }
    if (r.group.empty()) fail(std::format("[reaction_{}] has no group", r.name), at);
    if (!file.find_group(r.group)) fail(std::format("[reaction_{}] refers to the undefined group '{}'", r.name, r.group), at);
  }
  if (file.general.interval < 1) fail("interval must be at least 1", 0);
  return file;
}

ReactionFile read_reaction_file(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error(std::format("cannot open the reaction file {}", path.string()));
  return parse_reaction_file(in, path.string());
}

}  // namespace polyreact
