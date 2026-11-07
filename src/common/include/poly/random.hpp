// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Random number source shared by all modules: xoshiro256++ (Blackman & Vigna, 2019), seeded
// through splitmix64. It passes the standard statistical test batteries, is several times
// faster than std::mt19937_64, and every 64-bit output can be split into two independent
// 32-bit halves.
#pragma once

#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <limits>
#include <random>
#include <utility>

namespace poly {

/// splitmix64 step: expands one 64-bit seed into the generator state.
[[nodiscard]] constexpr std::uint64_t splitmix64(std::uint64_t& x) noexcept {
  std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

/// xoshiro256++, a UniformRandomBitGenerator usable with <random> distributions.
class Xoshiro256pp {
 public:
  using result_type = std::uint64_t;

  explicit constexpr Xoshiro256pp(std::uint64_t seed) noexcept {
    for (auto& word : s_) word = splitmix64(seed);
  }

  static constexpr result_type min() noexcept { return 0; }
  static constexpr result_type max() noexcept { return std::numeric_limits<result_type>::max(); }

  constexpr result_type operator()() noexcept {
    const std::uint64_t result = std::rotl(s_[0] + s_[3], 23) + s_[0];
    const std::uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = std::rotl(s_[3], 45);
    return result;
  }

 private:
  std::array<std::uint64_t, 4> s_{};
};

class Rng {
 public:
  explicit Rng(std::uint64_t seed) : seed_(seed), engine_(seed) {}

  /// A fresh seed from std::random_device mixed with the wall clock.
  [[nodiscard]] static std::uint64_t entropy_seed() {
    std::random_device rd;
    std::uint64_t seed = (static_cast<std::uint64_t>(rd()) << 32) ^ static_cast<std::uint64_t>(rd());
    const auto now = static_cast<std::uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    seed ^= now * 0x9E3779B97F4A7C15ULL;
    if (seed == 0) seed = 0x5DEECE66DULL;
    return seed;
  }

  [[nodiscard]] std::uint64_t seed() const noexcept { return seed_; }

  /// Uniform double in [0, 1) with 53 random bits.
  [[nodiscard]] double uniform01() noexcept { return static_cast<double>(engine_() >> 11) * 0x1.0p-53; }

  /// Uniform integer in [0, n), n >= 1. Lemire's nearly-divisionless method.
  [[nodiscard]] std::uint32_t below(std::uint32_t n) noexcept {
    return bounded(static_cast<std::uint32_t>(engine_() >> 32), n);
  }

  /// Two independent uniform integers in [0, n1) and [0, n2) from a single 64-bit draw.
  [[nodiscard]] std::pair<std::uint32_t, std::uint32_t> below2(std::uint32_t n1, std::uint32_t n2) noexcept {
    const std::uint64_t r = engine_();
    return {bounded(static_cast<std::uint32_t>(r >> 32), n1), bounded(static_cast<std::uint32_t>(r), n2)};
  }

  [[nodiscard]] bool bernoulli(double p) noexcept { return uniform01() < p; }

  [[nodiscard]] Xoshiro256pp& engine() noexcept { return engine_; }

 private:
  /// Map 32 random bits to [0, n) without bias; redraws on the rare rejection.
  [[nodiscard]] std::uint32_t bounded(std::uint32_t x, std::uint32_t n) noexcept {
    std::uint64_t m = static_cast<std::uint64_t>(x) * n;
    auto low = static_cast<std::uint32_t>(m);
    if (low < n) {
      const std::uint32_t threshold = (0u - n) % n;
      while (low < threshold) {
        m = static_cast<std::uint64_t>(static_cast<std::uint32_t>(engine_() >> 32)) * n;
        low = static_cast<std::uint32_t>(m);
      }
    }
    return static_cast<std::uint32_t>(m >> 32);
  }

  std::uint64_t seed_;
  Xoshiro256pp engine_;
};

}  // namespace poly
