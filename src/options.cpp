#include "euler.hpp"

#include <cstdlib>
#include <iostream>
#include <limits>

namespace euler {
Options ParseOptions(int argc, char **argv) {
  Options o;
  bool assembly_selected = false;
  bool requested_matrix_free = false;
  for (int i = 1; i < argc; i++) {
    std::string s = argv[i];
    auto value = [&]() {
      require(i + 1 < argc, "Missing option value");
      return std::string(argv[++i]);
    };
    if (s == "--hdg-first") {
      o.hdg_first = true;
    } else if (s == "--nx") {
      o.nx = std::stoi(value());
    } else if (s == "--ny") {
      o.ny = std::stoi(value());
    } else if (s == "--k") {
      o.k = std::stoi(value());
    } else if (s == "--cfl") {
      o.cfl = std::stod(value());
      require(o.cfl > 0, "--cfl must be positive");
    } else if (s == "--steps") {
      o.steps = std::stoi(value());
      require(o.steps > 0, "--steps must be positive");
    } else if (s == "--newton-steps") {
      o.newton_steps = std::stoi(value());
      require(o.newton_steps > 0, "--newton-steps must be positive");
    } else if (s == "--check-update") {
      o.check_update = true;
    } else if (s == "--timing-only") {
      o.timing_only = true;
    } else if (s == "--ilu-ordering") {
      o.ilu_ordering = value();
    } else if (s == "--dt") {
      o.dt = std::stod(value());
    } else if (s == "--T") {
      o.end = std::stod(value());
    } else if (s == "--alpha-deg") {
      o.alpha = std::stod(value());
    } else if (s == "--safety") {
      o.safety = std::stod(value());
    } else if (s == "--volume-order") {
      o.volume_order = std::stoi(value());
    } else if (s == "--assembly") {
      o.assembly = value();
      assembly_selected = true;
    } else if (s == "--form") {
      o.form = value();
    } else if (s == "--solver") {
      o.solver = value();
    } else if (s == "--precond") {
      o.precond = value();
    } else if (s == "--output") {
      o.output = value();
    } else if (s == "--matrix-free") {
      requested_matrix_free = true;
    } else if (s == "--no-save-solution") {
      o.save_solution = false;
    } else if (s == "--linear-tol") {
      o.linear_tol = std::stod(value());
    } else if (s == "--newton-tol") {
      o.tol = std::stod(value());
    } else if (s == "--restart") {
      o.restart = std::stoi(value());
    } else if (s == "--max-krylov") {
      o.max_krylov = std::stoi(value());
    } else if (s == "--max-newton") {
      o.max_newton = std::stoi(value());
    } else if (s == "--verify") {
      o.verify = true;
    } else if (s == "--help") {
      std::cout <<
          "--nx N --ny N --k 0..8 --alpha-deg A (default 45)\n"
          "--dt DT | --cfl C (CFL wins; h=min(hx,hy), speed includes safety)\n"
          "--T T (default .2) | --steps N (steps wins; no T clipping)\n"
          "--max-newton N --newton-tol E | --newton-steps N (fixed updates/step)\n"
          "--form dg|hdg|both --solver gmres|direct|dense\n"
          "--precond block-jacobi|block-ilu0 --ilu-ordering none|mdf\n"
          "--assembly full|element|partial --matrix-free --safety 1.1\n"
          "--linear-tol E --restart N --max-krylov N --output PREFIX\n"
          "--check-update --timing-only --no-save-solution --hdg-first --verify\n";
      std::exit(EXIT_SUCCESS);
    } else {
      throw std::runtime_error("Unknown option: " + s);
    }
  }
  require(o.nx >= 2 && o.ny >= 2 && o.k >= 0 && o.k <= 8 && o.dt > 0 && o.end >= 0 &&
              o.end < 1 / pi && o.safety >= 1 && o.tol > 0 && o.linear_tol > 0,
          "Invalid mesh, degree, time or safety");
  long double cells = static_cast<long double>(o.nx) * o.ny, q = o.k + 1.L;
  require(4 * cells * q * q <= std::numeric_limits<int>::max() &&
              8 * cells * q <= std::numeric_limits<int>::max(),
          "Mesh/degree exceeds MFEM int DOF capacity");
  require(std::isfinite(o.alpha) && std::isfinite(o.dt) && std::isfinite(o.end) &&
              std::isfinite(o.safety) && std::isfinite(o.tol) &&
              std::isfinite(o.linear_tol),
          "Nonfinite option");
  require(o.form == "dg" || o.form == "hdg" || o.form == "both", "Invalid form");
  require(o.solver == "gmres" || o.solver == "direct" || o.solver == "dense",
          "Invalid solver");
  require(o.restart > 0 && o.max_krylov > 0 && o.max_newton > 0,
          "Iteration limits must be positive");
  require((o.cfl == -1 || (std::isfinite(o.cfl) && o.cfl > 0)) &&
              (o.steps == -1 || o.steps > 0) &&
              (o.newton_steps == -1 || o.newton_steps > 0),
          "CFL, steps and fixed Newton updates must be positive when supplied");
  require(o.alpha > 0 && o.alpha < 90, "Require 0 < alpha-deg < 90");
  const double angle = o.alpha*pi/180;
  require(std::isfinite(2/std::cos(angle)) && std::isfinite(2/std::sin(angle)),
          "Angle produces nonfinite domain lengths");
  require(o.ilu_ordering == "none" || o.ilu_ordering == "mdf", "Invalid ILU ordering");
  // Keep the requested method collocated. Other quadratures are different
  // discretizations and should not silently enter the timing comparison.
  require(o.volume_order == -1 || o.volume_order == 2 * o.k + 1,
          "Volume rule must be (k+1)^2 GL points (order 2k+1)");
  o.volume_order = -1;
  if (o.solver == "gmres") {
    require(o.precond == "block-jacobi" || o.precond == "block-ilu0",
            "Invalid preconditioner");
    if (requested_matrix_free) {
      o.assembly = "partial";
    }
    o.matrix_free = o.assembly != "full";
  } else {
    require(!requested_matrix_free, "Direct LU requires an assembled matrix");
    if (!assembly_selected) {
      o.assembly = "full";
    }
    require(o.assembly == "full", "Direct LU requires --assembly full");
    o.matrix_free = false;
    o.precond = "none";
  }
  require(o.assembly == "full" || o.assembly == "element" || o.assembly == "partial",
          "Invalid assembly");
  if (o.timing_only) { o.save_solution = false; }
  return o;
}
} // namespace euler
