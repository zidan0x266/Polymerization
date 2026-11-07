// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Configuration-file reader shared by all modules. Files hold `key = value` lines with optional
// [section] headers and `#` or `;` comments; `key=value` arguments on the command line override
// them. Each module's Config supplies a setter that parses one key and returns an error message
// instead of throwing, so that errors can be reported with their file and line.
#pragma once

#include <algorithm>
#include <cctype>
#include <expected>
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace poly::ini {

using Error = std::string;
template <class T>
using Parsed = std::expected<T, Error>;

[[nodiscard]] inline std::string trim(std::string_view s) {
  const auto notspace = [](unsigned char c) { return !std::isspace(c); };
  const auto b = std::find_if(s.begin(), s.end(), notspace);
  const auto e = std::find_if(s.rbegin(), s.rend(), notspace).base();
  return (b < e) ? std::string(b, e) : std::string();
}

[[nodiscard]] inline std::string lower(std::string s) {
  std::ranges::transform(s, s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

/// Lowercase key without its optional section prefix ("Run.Seed" -> "seed").
[[nodiscard]] inline std::string normalize_key(std::string_view raw) {
  std::string key = lower(trim(raw));
  if (const auto dot = key.rfind('.'); dot != std::string::npos) key = key.substr(dot + 1);
  return key;
}

[[nodiscard]] inline Parsed<bool> parse_bool(std::string_view key, std::string_view value) {
  const std::string v = lower(trim(value));
  if (v == "true" || v == "yes" || v == "on" || v == "1") return true;
  if (v == "false" || v == "no" || v == "off" || v == "0") return false;
  return std::unexpected(std::format("config key '{}': expected a boolean, got '{}'", key, value));
}

template <class T>
[[nodiscard]] Parsed<T> parse_number(std::string_view key, std::string_view value) {
  std::istringstream in{trim(value)};
  T out{};
  in >> out;
  if (in.fail()) return std::unexpected(std::format("config key '{}': cannot parse '{}'", key, value));
  std::string rest;
  in >> rest;
  if (!rest.empty()) return std::unexpected(std::format("config key '{}': trailing text '{}'", key, rest));
  return out;
}

/// Comma-separated list of numbers; empty entries are skipped.
[[nodiscard]] inline Parsed<std::vector<double>> parse_double_list(std::string_view key, std::string_view value) {
  std::vector<double> out;
  std::string token;
  std::istringstream in{trim(value)};
  while (std::getline(in, token, ',')) {
    const std::string t = trim(token);
    if (t.empty()) continue;
    auto x = parse_number<double>(key, t);
    if (!x) return std::unexpected(std::move(x.error()));
    out.push_back(*x);
  }
  return out;
}

/// Assign a parsed value or propagate its error.
template <class T, class U>
std::expected<void, Error> assign(T& target, Parsed<U>&& parsed) {
  if (!parsed) return std::unexpected(std::move(parsed.error()));
  target = static_cast<T>(*parsed);
  return {};
}

/// Read a configuration file, calling `set(key, value)` for every entry. Unknown keys and
/// malformed lines throw std::invalid_argument prefixed with file and line.
template <class Setter>
void read_file(const std::filesystem::path& path, Setter&& set) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error(std::format("cannot open configuration file {}", path.string()));
  std::string line;
  int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    if (const auto hash = line.find('#'); hash != std::string::npos) line.erase(hash);
    if (const auto semi = line.find(';'); semi != std::string::npos) line.erase(semi);
    const std::string t = trim(line);
    if (t.empty()) continue;
    if (t.front() == '[' && t.back() == ']') continue;  // sections are decorative
    const auto eq = t.find('=');
    if (eq == std::string::npos) {
      throw std::invalid_argument(std::format("{}:{}: expected 'key = value'", path.string(), lineno));
    }
    if (auto r = set(std::string_view(t).substr(0, eq), std::string_view(t).substr(eq + 1)); !r) {
      throw std::invalid_argument(std::format("{}:{}: {}", path.string(), lineno, r.error()));
    }
  }
}

/// Apply `key=value` command-line overrides (a leading `--` on the key is ignored).
template <class Setter>
void apply_overrides(const std::vector<std::string>& overrides, Setter&& set) {
  for (const std::string& o : overrides) {
    const auto eq = o.find('=');
    if (eq == std::string::npos) throw std::invalid_argument(std::format("override '{}' is not of the form key=value", o));
    std::string_view key = std::string_view(o).substr(0, eq);
    if (key.starts_with("--")) key.remove_prefix(2);
    if (auto r = set(key, std::string_view(o).substr(eq + 1)); !r) throw std::invalid_argument(r.error());
  }
}

}  // namespace poly::ini
