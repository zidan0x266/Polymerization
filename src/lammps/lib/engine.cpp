// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2026 Zidan Zhang <zhangzidan@gmail.com>

#include "polyreact/engine.hpp"

#include <algorithm>
#include <format>
#include <stdexcept>

namespace polyreact {

Engine::Engine(const ReactionFile& file, const std::map<std::string, int, std::less<>>& types, int num_types) : num_types_(num_types) {
  build(file, types);
}

Engine::Engine(const ReactionFile& file) {
  std::map<std::string, int, std::less<>> types;
  for (const std::string& name : file.type_names()) types.emplace(name, static_cast<int>(types.size()));
  num_types_ = static_cast<int>(types.size());
  build(file, types);
}

int Engine::type_index(std::string_view name) const noexcept {
  const auto it = std::ranges::find(names_, name);
  return it == names_.end() ? -1 : static_cast<int>(it - names_.begin());
}

bool Engine::context_allows(const Rule& r, const PairContext& c) noexcept {
  if (c.distance < r.min_cutoff || c.distance >= r.cutoff) return false;
  if (r.kind == Kind::Dissociation) return c.bonded;
  if (c.bonded && !r.virtual_bond) return false;  // no second bond between the same beads
  if (c.same_molecule && !r.intramolecular) return false;
  if (c.same_residue && !r.intraresidual) return false;
  return true;
}

void Engine::build(const ReactionFile& file, const std::map<std::string, int, std::less<>>& types) {
  if (num_types_ < 1) throw std::invalid_argument("the reaction file needs at least one bead type");
  names_.assign(static_cast<std::size_t>(num_types_), std::string());
  for (const auto& [name, index] : types) {
    if (index < 0 || index >= num_types_) throw std::invalid_argument(std::format("type '{}' has the index {} out of range", name, index));
    names_[static_cast<std::size_t>(index)] = name;
  }
  const auto index_of = [&](const std::string& name, const Reaction& r) {
    const auto it = types.find(name);
    if (it == types.end()) throw std::invalid_argument(std::format("[reaction_{}]: unknown bead type '{}'", r.name, name));
    return it->second;
  };
  competition_ = file.general.competition;
  pair_table_.assign(static_cast<std::size_t>(num_types_) * static_cast<std::size_t>(num_types_), {});
  single_table_.assign(static_cast<std::size_t>(num_types_), {});
  reactive_.assign(static_cast<std::size_t>(num_types_), 0);

  for (std::size_t i = 0; i < file.reactions.size(); ++i) {
    const Reaction& reaction = file.reactions[i];
    if (!reaction.active) continue;
    const Equation& eq = reaction.equation;
    if (eq.kind == Kind::Exchange) {
      throw std::invalid_argument(std::format("[reaction_{}]: exchange reactions are not supported yet", reaction.name));
    }
    Rule r;
    r.reaction = i;
    r.kind = eq.kind;
    r.rate = eq.kind == Kind::Dissociation ? reaction.diss_rate.value_or(reaction.rate) : reaction.rate;
    r.cutoff = reaction.cutoff;
    r.min_cutoff = reaction.min_cutoff;
    r.intramolecular = reaction.intramolecular;
    r.intraresidual = reaction.intraresidual;
    r.virtual_bond = reaction.virtual_bond;
    const auto group = std::ranges::find(file.groups, reaction.group, &Group::name);
    r.group = static_cast<std::size_t>(group - file.groups.begin());
    for (std::size_t s = 0; s < eq.reactants.size(); ++s) {
      r.type[s] = index_of(eq.reactants[s].type, reaction);
      r.min[s] = eq.reactants[s].min;
      r.max[s] = eq.reactants[s].max;
      r.new_type[s] = index_of(eq.products[s].type, reaction);
      r.delta[s] = eq.products[s].delta;
    }
    const auto index = static_cast<std::uint32_t>(rules_.size());
    rules_.push_back(r);
    if (eq.kind == Kind::Unimolecular) {
      single_table_[static_cast<std::size_t>(r.type[0])].push_back(index);
      has_single_ = true;
      continue;
    }
    const auto a = static_cast<std::size_t>(r.type[0]), b = static_cast<std::size_t>(r.type[1]);
    const auto n = static_cast<std::size_t>(num_types_);
    pair_table_[a * n + b].push_back(index);
    if (a != b) pair_table_[b * n + a].push_back(index);
    reactive_[a] = reactive_[b] = 1;
    max_cutoff_ = std::max(max_cutoff_, r.cutoff);
  }
  for (const auto& list : pair_table_) {
    if (list.size() > kMaxCandidates) throw std::invalid_argument("more than 32 reactions for one pair of bead types");
  }
  for (const auto& list : single_table_) {
    if (list.size() > kMaxCandidates) throw std::invalid_argument("more than 32 unimolecular reactions for one bead type");
  }
}

}  // namespace polyreact
