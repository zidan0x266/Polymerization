// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Output files that report failure: every writer of the modules opens its file through here.
#pragma once

#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <stdexcept>

namespace poly::io {

/// Open `path` for writing in `mode`; throws std::runtime_error if it cannot be opened.
[[nodiscard]] inline std::ofstream open_output(const std::filesystem::path& path, std::ios::openmode mode = std::ios::out) {
  std::ofstream out(path, mode);
  if (!out) throw std::runtime_error(std::format("cannot open {} for writing", path.string()));
  return out;
}

/// Throw std::runtime_error if writing `out`, the stream of `path`, failed.
inline void check_written(const std::ofstream& out, const std::filesystem::path& path) {
  if (!out) throw std::runtime_error(std::format("error while writing {}", path.string()));
}

}  // namespace poly::io
