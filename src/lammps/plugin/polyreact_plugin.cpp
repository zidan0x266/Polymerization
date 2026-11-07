// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Registers fix poly/react with LAMMPS: plugin load polyreactplugin.so

#include "lammpsplugin.h"
#include "version.h"

#include "fix_poly_react.h"

using namespace LAMMPS_NS;

static Fix *create_fix(LAMMPS *lmp, int argc, char **argv)
{
  return new FixPolyReact(lmp, argc, argv);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  auto register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;
  plugin.style = "fix";
  plugin.name = "poly/react";
  plugin.info = "fix poly/react: polymerization reactions from a reaction file";
  plugin.author = "Zidan Zhang (zhangzidan@gmail.com)";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &create_fix;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
