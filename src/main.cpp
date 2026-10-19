#include "euler.hpp"
#include <filesystem>
#include <iostream>
int main(int argc, char **argv) {
  using namespace euler;
  try {
    Options o = ParseOptions(argc, argv);
    const EulerGas gas(3.); // Gamma belongs to the selected test problem.
    auto parent = std::filesystem::path(o.output).parent_path();
    if (!parent.empty()) {
      std::filesystem::create_directories(parent);
    }
    if (o.verify) {
      Verify(o, gas);
      return 0;
    }
    Problem p(o, gas);
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
      const auto l2 = p.Differences(dg.u, hdg.u);
      double combined_l2 = 0.;
      for (double value : l2) { combined_l2 += value*value; }
      combined_l2 = std::sqrt(combined_l2);
      std::ofstream f(o.output + "_comparison.csv");
      f << "max_nodal_difference,rho_L2_difference,rhou_L2_difference,"
           "rhov_L2_difference,energy_L2_difference,solution_L2_difference,"
           "dg_over_hdg_mean_newton\n"
        << std::setprecision(17) << diff << ',' << l2[0] << ',' << l2[1] << ','
        << l2[2] << ',' << l2[3] << ',' << combined_l2 << ','
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
