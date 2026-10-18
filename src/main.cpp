#include "euler.hpp"
#include <filesystem>
#include <iostream>
int main(int argc, char **argv) {
  using namespace euler;
  try {
    Options o = ParseOptions(argc, argv);
    auto parent = std::filesystem::path(o.output).parent_path();
    if (!parent.empty()) {
      std::filesystem::create_directories(parent);
    }
    if (o.verify) {
      Verify(o);
      return 0;
    }
    Problem p(o);
    RunResult dg, hdg;
    if (o.form == "both" && o.hdg_first) {
      hdg = Run(p, "hdg", o.output + "_hdg");
      dg = Run(p, "dg", o.output + "_dg");
    } else {
      if (o.form != "hdg") {
        dg = Run(p, "dg", o.output + "_dg");
      }
      if (o.form != "dg") {
        hdg = Run(p, "hdg", o.output + "_hdg");
      }
    }
    if (o.form == "both") {
      Vector d = dg.u;
      d -= hdg.u;
      double diff = d.Normlinf();
      std::ofstream f(o.output + "_comparison.csv");
      f << "max_nodal_difference,rho_L2_difference,dg_over_hdg_mean_newton\n"
        << std::setprecision(17) << diff << ',' << p.Difference(dg.u, hdg.u) << ','
        << (dg.updates && hdg.updates && hdg.newton_seconds > 0
               ? (dg.newton_seconds / dg.updates) / (hdg.newton_seconds / hdg.updates) : 0.) << '\n';
      std::cout << "DG/HDG max nodal difference=" << diff << '\n';
      require(diff < 1e-7, "DG/HDG final solutions disagree");
    }
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << '\n';
    return 1;
  }
}
