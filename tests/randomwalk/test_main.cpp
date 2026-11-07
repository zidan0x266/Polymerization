// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Self-contained unit tests (no framework dependency). The most important one checks the
// probabilities with which a walk picks its next site against an independent reimplementation
// of the original program's neighbor scan.

#include <cmath>
#include <cstdio>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "rw/io/config.hpp"
#include "rw/io/writers.hpp"
#include "rw/lattice/fcc.hpp"
#include "rw/model/packing.hpp"
#include "rw/walk/packer.hpp"

namespace {

int g_failures = 0;
int g_checks = 0;

#define CHECK(cond)                                                                                  \
  do {                                                                                               \
    ++g_checks;                                                                                      \
    if (!(cond)) {                                                                                   \
      ++g_failures;                                                                                  \
      std::fprintf(stderr, "CHECK failed: %s  (%s:%d)\n", #cond, __FILE__, __LINE__);                \
    }                                                                                                \
  } while (0)

using namespace rw;

bool adjacent(Vec3i a, Vec3i b, int extent) {
  const auto d = [extent](int u) { return ((u % extent) + extent + extent / 2) % extent - extent / 2; };
  const Vec3i v{d(a.x - b.x), d(a.y - b.y), d(a.z - b.z)};
  return v.norm2() == 2;
}

// ------------------------------------------------------------------------------------------------
void test_lattice() {
  for (int n : {2, 3, 5}) {
    FccLattice lattice(n);
    CHECK(lattice.num_sites() == static_cast<std::size_t>(4 * n * n * n));
    std::set<std::tuple<int, int, int>> points;
    for (Site s = 0; s < lattice.num_sites(); ++s) {
      const Vec3i p = lattice.position(s);
      CHECK((p.x + p.y + p.z) % 2 == 0);
      CHECK(p.x >= 0 && p.y >= 0 && p.z >= 0 && p.x < 2 * n && p.y < 2 * n && p.z < 2 * n);
      CHECK(lattice.site(p) == s);
      points.insert({p.x, p.y, p.z});
      std::set<Site> neighbors;
      for (int k = 0; k < FccLattice::kCoordination; ++k) neighbors.insert(lattice.neighbor(p, k));
      CHECK(neighbors.size() == 12 && !neighbors.contains(s));
    }
    CHECK(points.size() == lattice.num_sites());
    CHECK(lattice.wrap({-1, 2 * n, 3}) == (Vec3i{2 * n - 1, 0, 3}));
  }
  for (const Vec3i& d : FccLattice::kNeighbors) CHECK(d.norm2() == 2);

  // neighbor counts follow occupy/vacate exactly
  FccLattice lattice(3);
  std::vector<Site> filled;
  for (Site s = 0; s < lattice.num_sites(); s += 7) {
    lattice.occupy(s);
    filled.push_back(s);
  }
  for (std::size_t i = 0; i < filled.size(); i += 3) lattice.vacate(filled[i]);
  std::size_t count = 0;
  bool counts_ok = true;
  for (Site s = 0; s < lattice.num_sites(); ++s) {
    int occupied_neighbors = 0;
    for (Site t = 0; t < lattice.num_sites(); ++t) {
      if (lattice.occupied(t) && adjacent(lattice.position(s), lattice.position(t), lattice.extent())) ++occupied_neighbors;
    }
    counts_ok = counts_ok && occupied_neighbors == lattice.occupied_neighbors(s);
    count += lattice.occupied(s);
  }
  CHECK(counts_ok);
  CHECK(count == lattice.occupied_count());
}

// ------------------------------------------------------------------------------------------------
void test_angles() {
  const AngleSet all{60, 90, 120, 180};
  for (const Vec3i& p : FccLattice::kNeighbors) {
    for (const Vec3i& q : FccLattice::kNeighbors) {
      // angle between the bond back to the previous bead (-p) and the bond forward (q)
      const double angle = std::acos(-p.dot(q) / 2.0) * 180.0 / 3.14159265358979323846;
      const int degrees = static_cast<int>(std::lround(angle));
      CHECK(all.allows_steps(p, q) == (degrees != 0));
      for (int a : AngleSet::kValues) CHECK(AngleSet{a}.allows_steps(p, q) == (degrees == a));
    }
  }
  AngleSet set;
  CHECK(set.empty());
  CHECK(!set.allow(45) && set.empty());
  CHECK(set.allow(90) && set.allows(90) && !set.allows(120));
  CHECK((WalkRules{}.angles == AngleSet{60, 90, 120}));
}

// ------------------------------------------------------------------------------------------------
// The original program's loops: cubic cell offsets x, y, z in -1, 0, 1, then basis sites
// (0,0,0), (1/2,1/2,0), (1/2,0,1/2), (0,1/2,1/2); nearest neighbors kept in that order.
std::vector<Vec3i> original_scan(Vec3i head) {
  const Vec3i cell{head.x / 2 * 2, head.y / 2 * 2, head.z / 2 * 2};  // doubled cell origin
  const Vec3i basis[4] = {{0, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}};
  std::vector<Vec3i> order;
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      for (int k = -1; k <= 1; ++k) {
        for (const Vec3i& b : basis) {
          const Vec3i site = cell + Vec3i{2 * i, 2 * j, 2 * k} + b;
          if ((site - head).norm2() == 2) order.push_back(site - head);
        }
      }
    }
  }
  return order;
}

void test_scan_order() {
  for (const Vec3i head : {Vec3i{4, 4, 4}, Vec3i{5, 5, 4}, Vec3i{5, 4, 5}, Vec3i{4, 5, 5}}) {
    const auto expected = original_scan(head);
    CHECK(expected.size() == 12);
    const auto& order = legacy_scan_order(head);
    bool same = true;
    for (int i = 0; i < 12; ++i) same = same && FccLattice::kNeighbors[order[i]] == expected[i];
    CHECK(same);
  }
}

// ------------------------------------------------------------------------------------------------
// Selection probabilities of the third bead of a chain on an empty lattice, for every weighting,
// against probabilities computed directly from the rules. The first two beads are both
// occupied and both have served as heads.
std::map<int, double> expected_third_step(SiteWeighting w, AngleSet angles, Vec3i b0, Vec3i b1, int extent) {
  const Vec3i previous = b1 - b0;
  // empty neighbors in scan order (the order only matters for legacy)
  const std::vector<Vec3i> scan = w == SiteWeighting::Legacy
                                      ? original_scan(b1)
                                      : std::vector<Vec3i>(FccLattice::kNeighbors.begin(), FccLattice::kNeighbors.end());
  std::vector<Vec3i> free_sites;
  for (const Vec3i& d : scan) {
    if (d != -previous) free_sites.push_back(d);
  }
  const auto weight = [&](Vec3i d) {
    if (w == SiteWeighting::Uniform) return 1.0;
    const Vec3i s = b1 + d;
    return 12.0 - adjacent(s, b0, extent) - adjacent(s, b1, extent);
  };
  std::vector<Vec3i> allowed;
  for (const Vec3i& d : free_sites) {
    if (angles.allows_steps(previous, d)) allowed.push_back(d);
  }
  std::map<int, double> p;
  double total = 0;
  for (std::size_t i = 0; i < allowed.size(); ++i) {
    const double wt = (w == SiteWeighting::Legacy) ? weight(free_sites[i]) : weight(allowed[i]);
    const Vec3i d = allowed[i];
    p[(d.x + 1) * 9 + (d.y + 1) * 3 + (d.z + 1)] += wt;
    total += wt;
  }
  for (auto& [key, value] : p) value /= total;
  return p;
}

void test_selection_probabilities() {
  struct Case {
    SiteWeighting weighting;
    AngleSet angles;
  };
  for (const Case c : {Case{SiteWeighting::Legacy, {60, 90, 120}}, Case{SiteWeighting::Legacy, {90, 120}},
                       Case{SiteWeighting::FreeNeighbors, {60, 90, 120}}, Case{SiteWeighting::Uniform, {90, 120, 180}}}) {
    poly::Rng rng(2026);
    // counts[(head parity, previous step)][next step]
    std::map<std::pair<int, int>, std::map<int, double>> counts;
    std::map<std::pair<int, int>, std::pair<Vec3i, Vec3i>> geometry;
    const auto key = [](Vec3i d) { return (d.x + 1) * 9 + (d.y + 1) * 3 + (d.z + 1); };
    constexpr int kSamples = 300000;
    for (int i = 0; i < kSamples; ++i) {
      Packer packer(4, WalkRules{c.angles, c.weighting, 1}, rng);
      std::vector<Vec3i> beads(3);
      if (!packer.grow_chain(beads)) continue;
      const int parity = (beads[1].x & 1) | ((beads[1].y & 1) << 1);
      const std::pair<int, int> group{parity, key(beads[1] - beads[0])};
      counts[group][key(beads[2] - beads[1])] += 1;
      geometry.try_emplace(group, beads[0], beads[1]);
    }
    // chi-square over all (group, next step) cells
    double chi2 = 0;
    int dof = 0;
    for (const auto& [group, observed] : counts) {
      double n = 0;
      for (const auto& [k, v] : observed) n += v;
      const auto [b0, b1] = geometry[group];
      const auto expected = expected_third_step(c.weighting, c.angles, b0, b1, 8);
      for (const auto& [k, v] : observed) CHECK(expected.contains(k));  // only allowed steps occur
      for (const auto& [k, p] : expected) {
        const double o = observed.contains(k) ? observed.at(k) : 0.0;
        chi2 += (o - n * p) * (o - n * p) / (n * p);
      }
      dof += static_cast<int>(expected.size()) - 1;
    }
    CHECK(dof > 100);
    CHECK(chi2 / dof < 1.35);  // about 5 standard deviations for dof > 300
    if (chi2 / dof >= 1.35) {
      std::fprintf(stderr, "  chi2/dof = %.3f (dof %d) for weighting %s\n", chi2 / dof, dof,
                   std::string(to_string(c.weighting)).c_str());
    }
  }
}

// ------------------------------------------------------------------------------------------------
void test_box_size() {
  CHECK(cells_for_occupancy(18 * 5680, 0.05) == 80);  // the paper system: 113.137 box edge
  for (std::size_t n : {100u, 400u, 4000u, 99328u, 102240u, 106280u}) {
    for (double occ : {0.01, 0.05, 0.1, 0.3}) {
      const int cells = cells_for_occupancy(n, occ);
      CHECK(static_cast<double>(n) / (4.0 * cells * cells * cells) < occ);
      CHECK(cells == 1 || static_cast<double>(n) / (4.0 * (cells - 1) * (cells - 1) * (cells - 1)) >= occ);
    }
  }
  CHECK(cells_for_occupancy(400, 0.1) == 11);  // exact cube: the occupancy must be strictly below
}

// ------------------------------------------------------------------------------------------------
void test_packing() {
  for (SiteWeighting w : {SiteWeighting::Legacy, SiteWeighting::FreeNeighbors, SiteWeighting::Uniform}) {
    for (AngleSet angles : {AngleSet{60, 90, 120}, AngleSet{90}, AngleSet{180}, AngleSet{60, 90, 120, 180}}) {
      poly::Rng rng(7);
      Packer packer(10, WalkRules{angles, w, 1000}, rng);
      const Packing p = packer.pack(60, 12, 30, 0.97);
      CHECK(p.size() == 60 * 12 + 30 && p.n_bonds() == 60 * 11);
      CHECK(validate(p, angles).has_value());
      CHECK(packer.lattice().occupied_count() == p.size());
    }
  }

  poly::Rng a(11), b(11), c(12);
  const Packing pa = Packer(8, WalkRules{}, a).pack(20, 8, 10, 1.0);
  const Packing pb = Packer(8, WalkRules{}, b).pack(20, 8, 10, 1.0);
  const Packing pc = Packer(8, WalkRules{}, c).pack(20, 8, 10, 1.0);
  CHECK(pa.beads == pb.beads);
  CHECK(pa.beads != pc.beads);

  CHECK(pa.type(0) == BeadType::End && pa.type(1) == BeadType::Interior && pa.type(7) == BeadType::End);
  CHECK(pa.type(8) == BeadType::End && pa.type(160) == BeadType::Crosslinker);
  CHECK(pa.molecule(0) == 1 && pa.molecule(7) == 1 && pa.molecule(8) == 2 && pa.molecule(160) == 21 && pa.molecule(169) == 30);

  // the check catches broken packings
  const AngleSet angles = WalkRules{}.angles;
  Packing overlap = pa;
  overlap.beads[165] = overlap.beads[3];
  CHECK(!validate(overlap, angles).has_value());
  Packing stretched = pa;
  for (std::size_t v = 4; v < 8; ++v) stretched.beads[v] = stretched.beads[v] + Vec3i{4, 0, 0};
  CHECK(!validate(stretched, angles).has_value());
  CHECK(!validate(pa, AngleSet{180}).has_value());

  // a crowded lattice ends in a clear error instead of an endless loop
  poly::Rng rng(3);
  Packer crowded(2, WalkRules{AngleSet{180}, SiteWeighting::Uniform, 5}, rng);
  bool threw = false;
  try {
    (void)crowded.pack(2, 12, 0, 1.0);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK(threw);
}

// ------------------------------------------------------------------------------------------------
void test_data_file() {
  poly::Rng rng(5);
  const Packing p = Packer(4, WalkRules{}, rng).pack(6, 10, 3, 1.0);
  std::ostringstream data;
  write_lammps_data(data, p, 1.0, "test");
  const std::string text = data.str();
  std::istringstream in(text);
  std::string line;
  std::getline(in, line);
  CHECK(line == "test");
  std::map<std::string, std::size_t> header;
  std::vector<std::vector<double>> atoms;
  std::vector<std::vector<long>> bonds;
  std::string section;
  double box = 0;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream words(line);
    if (line.ends_with("xhi")) words >> box >> box;
    if (line.ends_with(" atoms") || line.ends_with(" bonds") || line.ends_with(" atom types") || line.ends_with(" bond types")) {
      std::size_t n = 0;
      words >> n;
      header[line.substr(line.find(' ') + 1)] = n;
      continue;
    }
    if (std::isalpha(static_cast<unsigned char>(line[0]))) {
      section = line;
      continue;
    }
    if (section == "Atoms # molecular") {
      std::vector<double> v(9);
      for (double& x : v) words >> x;
      atoms.push_back(v);
    } else if (section == "Bonds") {
      std::vector<long> v(4);
      for (long& x : v) words >> x;
      bonds.push_back(v);
    }
  }
  CHECK(header["atoms"] == 63 && header["bonds"] == 54 && header["atom types"] == 4 && header["bond types"] == 1);
  CHECK(atoms.size() == 63 && bonds.size() == 54);
  CHECK(std::fabs(box - p.box_length()) < 1e-6);
  bool positions_ok = true;
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const auto& a = atoms[i];
    positions_ok = positions_ok && a[0] == static_cast<double>(i + 1) && a[1] == static_cast<double>(p.molecule(i)) &&
                   a[2] == static_cast<double>(p.type(i));
    const double unwrapped[3] = {a[3] + a[6] * box, a[4] + a[7] * box, a[5] + a[8] * box};
    const Vec3i q = p.beads[i];
    positions_ok = positions_ok && std::fabs(unwrapped[0] - q.x * p.unit()) < 1e-5 &&
                   std::fabs(unwrapped[1] - q.y * p.unit()) < 1e-5 && std::fabs(unwrapped[2] - q.z * p.unit()) < 1e-5;
    positions_ok = positions_ok && a[3] >= 0 && a[3] < box && a[4] >= 0 && a[4] < box && a[5] >= 0 && a[5] < box;
  }
  CHECK(positions_ok);
  bool bonds_ok = true;
  for (const auto& b : bonds) {
    const auto i = static_cast<std::size_t>(b[2] - 1), j = static_cast<std::size_t>(b[3] - 1);
    bonds_ok = bonds_ok && b[1] == 1 && j == i + 1 && p.molecule(i) == p.molecule(j) &&
               (p.beads[j] - p.beads[i]).norm2() == 2;
  }
  CHECK(bonds_ok);
  CHECK(text.contains("Atom Type Labels\n\n1 bead\n2 end\n3 xlink\n4 xlink_full\n"));
}

// ------------------------------------------------------------------------------------------------
void test_config() {
  Config cfg;
  CHECK(cfg.box_cells() == 80);
  CHECK(cfg.set("bond_angles", " 90 ,120").has_value() && (cfg.walk.angles == AngleSet{90, 120}));
  CHECK(!cfg.set("bond_angles", "45").has_value());
  CHECK(!cfg.set("bond_angles", "").has_value());
  CHECK(cfg.set("walk.site_weighting", "Uniform").has_value() && cfg.walk.weighting == SiteWeighting::Uniform);
  CHECK(!cfg.set("site_weighting", "random").has_value());
  CHECK(!cfg.set("chain_lenght", "5").has_value());
  CHECK(cfg.set("cells", "12").has_value() && cfg.box_cells() == 12);
  cfg.apply_overrides({"chain_length=8", "--n_chains=10", "n_crosslinkers=5", "seed=99"});
  cfg.validate();

  // to_ini round trip
  const auto dir = std::filesystem::temp_directory_path() / "rw_test_config";
  std::filesystem::create_directories(dir);
  {
    std::ofstream out(dir / "c.ini");
    out << cfg.to_ini();
  }
  const Config back = Config::from_file(dir / "c.ini");
  CHECK(back.to_ini() == cfg.to_ini());
  std::filesystem::remove_all(dir);

  bool threw = false;
  try {
    Config bad;
    bad.chain_length = 1;
    bad.validate();
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);
  threw = false;
  try {
    Config full;
    full.cells = 5;  // 500 sites for 102240 beads
    full.validate();
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);
}

}  // namespace

int main() {
  test_lattice();
  test_angles();
  test_scan_order();
  test_selection_probabilities();
  test_box_size();
  test_packing();
  test_data_file();
  test_config();
  std::cout << g_checks - g_failures << " / " << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
