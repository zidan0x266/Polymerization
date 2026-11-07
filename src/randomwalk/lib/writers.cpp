// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "rw/io/writers.hpp"

#include <charconv>
#include <concepts>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rw {

namespace {

/// Floor division, for image flags of unwrapped coordinates.
int floor_div(int a, int b) noexcept { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); }

/// Text output through a fixed-size buffer, with integers formatted by std::to_chars.
class Buffered {
 public:
  explicit Buffered(std::ostream& out) : out_(out) { buf_.reserve(kChunk + 256); }
  ~Buffered() { flush(); }
  Buffered(const Buffered&) = delete;
  Buffered& operator=(const Buffered&) = delete;

  Buffered& operator<<(std::string_view s) {
    buf_.append(s);
    if (buf_.size() >= kChunk) flush();
    return *this;
  }
  Buffered& operator<<(char c) {
    buf_.push_back(c);
    return *this;
  }
  Buffered& operator<<(std::integral auto v) {
    char tmp[24];
    const auto end = std::to_chars(tmp, tmp + sizeof tmp, v).ptr;
    return *this << std::string_view(tmp, static_cast<std::size_t>(end - tmp));
  }
  void flush() {
    out_.write(buf_.data(), static_cast<std::streamsize>(buf_.size()));
    buf_.clear();
  }

 private:
  static constexpr std::size_t kChunk = std::size_t{1} << 20;
  std::ostream& out_;
  std::string buf_;
};

/// Positions are lattice points, so every coordinate is one of 2n values: format each once.
std::vector<std::string> coordinate_table(const Packing& packing) {
  std::vector<std::string> table(static_cast<std::size_t>(packing.extent()));
  for (std::size_t i = 0; i < table.size(); ++i) table[i] = std::format("{:.6f}", packing.unit() * static_cast<double>(i));
  return table;
}

std::ofstream open_output(const std::filesystem::path& path) {
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error(std::format("cannot write {}", path.string()));
  return out;
}

void check_written(const std::ofstream& out, const std::filesystem::path& path) {
  if (!out) throw std::runtime_error(std::format("error while writing {}", path.string()));
}

}  // namespace

void write_lammps_data(std::ostream& stream, const Packing& packing, double mass, std::string_view title) {
  const auto coord = coordinate_table(packing);
  const int extent = packing.extent();
  Buffered out(stream);

  out << title << "\n\n";
  out << packing.size() << " atoms\n" << kNumBeadTypes << " atom types\n";
  out << packing.n_bonds() << " bonds\n1 bond types\n\n";
  const std::string box = std::format("{:.6f}", packing.box_length());
  for (const char* axis : {"x", "y", "z"}) out << "0 " << box << ' ' << axis << "lo " << axis << "hi\n";

  out << "\nAtom Type Labels\n\n";
  for (int t = 0; t < kNumBeadTypes; ++t) out << t + 1 << ' ' << kTypeLabels[t] << '\n';
  out << "\nMasses\n\n";
  for (int t = 0; t < kNumBeadTypes; ++t) out << t + 1 << ' ' << std::format("{}", mass) << '\n';

  out << "\nAtoms # molecular\n\n";
  for (std::size_t i = 0; i < packing.size(); ++i) {
    const Vec3i p = packing.beads[i];
    const Vec3i image{floor_div(p.x, extent), floor_div(p.y, extent), floor_div(p.z, extent)};
    out << i + 1 << ' ' << packing.molecule(i) << ' ' << static_cast<int>(packing.type(i)) << ' ';
    out << coord[p.x - image.x * extent] << ' ' << coord[p.y - image.y * extent] << ' ' << coord[p.z - image.z * extent];
    out << ' ' << image.x << ' ' << image.y << ' ' << image.z << '\n';
  }

  out << "\nBonds\n\n";
  std::size_t bond = 0;
  for (std::size_t c = 0; c < packing.n_chains; ++c) {
    const std::size_t first = c * packing.chain_length + 1;  // one-based atom ID
    for (std::size_t v = 0; v + 1 < packing.chain_length; ++v) {
      out << ++bond << " 1 " << first + v << ' ' << first + v + 1 << '\n';
    }
  }
}

void write_lammps_data(const std::filesystem::path& path, const Packing& packing, double mass, std::string_view title) {
  auto out = open_output(path);
  write_lammps_data(out, packing, mass, title);
  check_written(out, path);
}

void write_xyz(std::ostream& stream, const Packing& packing) {
  const auto coord = coordinate_table(packing);
  const int extent = packing.extent();
  const auto wrap = [extent](int c) { return ((c % extent) + extent) % extent; };
  Buffered out(stream);
  const std::string box = std::format("{:.6f}", packing.box_length());
  out << packing.size() << '\n';
  out << "Lattice=\"" << box << " 0 0 0 " << box << " 0 0 0 " << box
      << "\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n";
  for (std::size_t i = 0; i < packing.size(); ++i) {
    const Vec3i p = packing.beads[i];
    out << static_cast<int>(packing.type(i)) << ' ' << coord[wrap(p.x)] << ' ' << coord[wrap(p.y)] << ' ' << coord[wrap(p.z)]
        << '\n';
  }
}

void write_xyz(const std::filesystem::path& path, const Packing& packing) {
  auto out = open_output(path);
  write_xyz(out, packing);
  check_written(out, path);
}

}  // namespace rw
