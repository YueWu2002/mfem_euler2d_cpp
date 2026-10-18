#include "euler.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
namespace euler {
namespace {
// With exact local recovery, epsilon_DG = -B D^{-1} epsilon_S.
// D is diagonal in the collocated trace basis. Avoid constructing a second
// global DG Jacobian solely to validate an HDG solve during timed runs.
void RecoveredResidual(const Problem &p, const Hybrid &h, Condensed &s,
                       const Vector &dl, Vector &result) {
  auto &eps = s.eps_scratch; auto &scaled = s.scaled_scratch;
  s.op->Mult(dl,eps); eps -= s.rhs;
  for (int f = 0; f < p.nf; ++f) {
    for (int i = 0; i < p.t; ++i) { scaled[f*p.t+i] = eps[f*p.t+i]/h.D[f](i,i); }
  }
  auto &local = s.local_scratch; auto &cell = s.cell_scratch;
  for (int e = 0; e < p.n; ++e) {
    for (int i = 0; i < 4*p.t; ++i) { local[i] = scaled[h.ids[e][i]]; }
    h.B[e].Mult(local,cell);
    for (int i = 0; i < p.m; ++i) { result[e*p.m+i] = -cell[i]; }
  }
}
}
RunResult Run(Problem &p, const std::string &form, const std::string &prefix) {
  const auto &o = p.opt;
  RunResult out;
  out.u = p.Initial();
  require(p.Admissible(out.u), "Initial projection not positive; refine mesh");
  const auto initial = p.Integral(out.u);
  std::ofstream log;
  if (!o.timing_only) {
    log.open(prefix + "_iterations.csv");
    require(bool(log), "Cannot open iteration output");
    log << "form,matrix_free,step,newton,time,dt,cfl,residual_physical_rms,damping,"
           "krylov_iterations,linear_relative_residual,dg_linear_relative_residual,"
           "residual_s,validation_s,jacobian_s,condensation_s,assembly_s,setup_s,"
           "solve_s,recovery_s,line_search_s,iteration_s,update_check_relative\n"
        << std::setprecision(17);
  }
  // Storage allocation is outside the time loop; every numerical refresh and
  // factorization is still performed and timed at the current Newton state.
  Vector residual(p.nu), rhs(p.nu), er(p.nu), du(p.nu), trial(p.nu), trial_residual(p.nu);
  Hybrid hybrid;
  Condensed condensed;
  if (form == "hdg") { p.PrepareHybrid(hybrid,condensed); }
  std::unique_ptr<ElementOperator> dg_jacobian;
  if (form == "dg" && o.assembly == "partial") { dg_jacobian = p.PartialDG(out.u,false); }
  const int global_size = form == "dg" ? p.nu : p.nl;
  LinearSolver linear_solver(global_size, form == "dg" ? p.m : p.t, o,
                             form == "dg" ? p.cell_metric : p.trace_metric);
  LinearResult lr, attempt_result;
  lr.x.SetSize(global_size); attempt_result.x.SetSize(global_size);
  Options linear = o;
  double time = 0, begin = clock_s(), dt_sum = 0., max_cfl = 0.;
  int step = 0;
  while (o.steps > 0 ? step < o.steps : time < o.end - 1e-14) {
    // Freeze once from the preceding time level; CFL selection reuses exactly
    // these speeds and does not refresh them inside Newton or line search.
    p.Freeze(out.u, o.dt);
    double rate = p.WaveRate();
    double dt = o.cfl > 0 ? o.cfl / rate : o.dt;
    if (o.steps < 0) { dt = std::min(dt, o.end - time); }
    require(std::isfinite(dt) && dt > 0 && time + dt < 1 / pi,
            "Step leaves the smooth simple-wave interval T < 1/pi");
    p.dt = dt;
    const double physical_scale = 1 / (dt * std::sqrt(p.lx * p.ly));
    const double old_norm = p.CellNorm(p.old);
    double r0 = -1.;
    bool converged = false, residual_ready = false;
    const int limit = o.newton_steps > 0 ? o.newton_steps : o.max_newton;
    for (int it = 0; it <= limit; ++it) {
      double istart = clock_s();
      if (!residual_ready) { p.Residual(out.u,residual); }
      const double residual_norm = p.CellNorm(residual);
      const double rms = physical_scale * residual_norm;
      if (r0 < 0) { r0 = rms; }
      if (o.newton_steps < 0 && rms <= o.tol * std::max(1., r0)) {
        converged = true;
        break;
      }
      if (it == limit) { break; }
      double residual_time = clock_s() - istart;
      double astart = clock_s(), jt = 0, ct = 0, rt = 0, vt = 0;
      double dgrel = 0, update_check = -1., diagnostic_time = 0.;
      lr.count = 0; lr.residual = lr.assembly = lr.setup = lr.solve = 0.;
      // A common mass norm criterion, including a roundoff allowance for small
      // Newton right-hand sides. All stored cell equations have the same scale.
      const double allowance = 128 * std::numeric_limits<double>::epsilon() *
          (p.CellNorm(out.u) + old_norm + residual_norm);
      const double target = std::max(10 * o.linear_tol * residual_norm, allowance);
      if (form == "dg") {
        if (o.assembly == "partial") { dg_jacobian->Refresh(out.u); }
        else { dg_jacobian = p.DG(out.u); }
        jt = clock_s() - astart;
        rhs = residual; rhs *= -1;
        linear_solver.Solve(*dg_jacobian,rhs,o,lr);
        du = lr.x;
        astart = clock_s();
        dg_jacobian->Mult(du, er);
        er += residual;
        vt = clock_s() - astart;
      } else {
        // Reconstruct Lambda=Phi(U) at each Newton state. On this manifold
        // r_H=r_DG and g=0 algebraically, so avoid redundant flux/residual passes.
        // Only S enters the global solver, and only local backsolves recover du.
        auto &h = hybrid; auto &s = condensed;
        p.HDG(out.u,h);
        jt = clock_s() - astart;
        astart = clock_s();
        p.Condense(h,residual,s);
        ct = clock_s() - astart;
        linear.linear_tol = o.linear_tol;
        bool valid = false;
        for (int attempt = 0; attempt < 4; ++attempt) {
          linear_solver.Solve(*s.op,s.rhs,linear,attempt_result,attempt == 0);
          lr.count += attempt_result.count;
          lr.assembly += attempt_result.assembly;
          lr.setup += attempt_result.setup;
          lr.solve += attempt_result.solve;
          lr.residual = attempt_result.residual;
          lr.x = attempt_result.x;
          astart = clock_s();
          p.Recover(s,lr.x,du);
          rt += clock_s() - astart;
          astart = clock_s();
          RecoveredResidual(p,h,s,lr.x,er);
          valid = p.CellNorm(er) <= target;
          vt += clock_s() - astart;
          if (valid) { break; }
          require(o.solver == "gmres", "Direct Schur solve fails DG residual criterion");
          linear.linear_tol *= .01;
        }
        require(valid, "Recovered DG residual criterion failed after Schur refinement");
      }
      dgrel = residual_norm > 0 ? p.CellNorm(er) / residual_norm : 0.;
      require(p.CellNorm(er) <= target, "DG Newton linear residual too large");
      if (o.check_update && form == "hdg") {
        astart = clock_s();
        auto dg = p.PartialDG(out.u);
        Vector independent(p.nu), rhs = residual;
        dg->Mult(du, independent);
        independent += residual;
        require(p.CellNorm(independent) <= std::max(4*target, 1e-12),
                "Independent recovered DG equation check failed");
        rhs *= -1;
        auto reference = Solve(*dg, rhs, p.m, o, p.cell_metric);
        Vector difference = du;
        difference -= reference.x;
        update_check = p.CellNorm(difference) / std::max(p.CellNorm(reference.x), 1e-30);
        require(p.CellNorm(difference) <=
                    std::max(1e-8 * p.CellNorm(reference.x), 4*allowance),
                "DG/HDG one-update equivalence failed");
        diagnostic_time = clock_s() - astart;
      }
      double ls = clock_s(), damping = 1.;
      bool accept = false;
      for (int back = 0; back < 30; ++back) {
        trial = out.u;
        trial.Add(damping, du);
        if (p.Admissible(trial)) {
          p.Residual(trial,trial_residual);
          if (p.CellNorm(trial_residual) <=
                (1 - 1e-4*damping)*residual_norm + allowance) {
            accept = true;
            break;
          }
        }
        damping *= .5;
      }
      require(accept, "Newton line search failed");
      out.u = trial;
      // Reuse the accepted forward evaluation at the next Newton state.
      // Its computation is already included in this update's line-search time.
      residual = trial_residual;
      residual_ready = true;
      double lst = clock_s() - ls;
      // Stop timing before any CSV formatting/writes. Optional independent
      // reference solves are diagnostics and excluded from Newton timings.
      double iteration_time = clock_s() - istart - diagnostic_time;
      ++out.updates;
      out.newton_seconds += iteration_time;
      if (log.is_open()) {
        log << form << ',' << o.matrix_free << ',' << step << ',' << it << ','
            << time + dt << ',' << dt << ',' << rate*dt << ',' << rms << ','
            << damping << ',' << lr.count << ',' << lr.residual << ',' << dgrel
            << ',' << residual_time << ',' << vt << ',' << jt << ',' << ct << ','
            << lr.assembly << ',' << lr.setup << ',' << lr.solve << ',' << rt << ','
            << lst << ',' << iteration_time << ',' << update_check << '\n';
      }
    }
    require(converged || o.newton_steps > 0, "Newton failed to converge");
    time += dt;
    dt_sum += dt;
    max_cfl = std::max(max_cfl, rate*dt);
    ++step;
  }
  out.seconds = clock_s() - begin;
  out.final_time = time;
  out.time_steps = step;
  out.errors = p.Errors(out.u, time);
  const auto final = p.Integral(out.u);
  for (int c = 0; c < 4; ++c) {
    out.drift = std::max(out.drift, std::abs(final[c] - initial[c]));
  }
  if (o.save_solution) { p.Save(out.u, prefix); }
  std::ofstream summary(prefix + "_summary.csv");
  require(bool(summary), "Cannot open summary output");
  summary << "form,nx,ny,k,gamma,alpha_deg,T,dt,cfl,actual_max_cfl,mean_dt,"
             "volume_dofs,trace_dofs,matrix_free,assembly,solver,preconditioner,"
             "ilu_ordering,requested_steps,fixed_newton_steps,time_steps,"
             "newton_updates,newton_s,mean_newton_s,time_loop_s,rho_L2,rhou_L2,"
             "rhov_L2,energy_L2,max_conservation_drift,check_update\n"
          << std::setprecision(17);
  const double mean = out.updates ? out.newton_seconds / out.updates : 0.;
  summary << form << ',' << o.nx << ',' << o.ny << ',' << o.k << ",3," << o.alpha << ','
          << time << ',' << o.dt << ',' << o.cfl << ',' << max_cfl << ','
          << (step ? dt_sum/step : 0.) << ',' << p.nu << ',' << p.nl << ','
          << o.matrix_free << ',' << o.assembly << ',' << o.solver << ',' << o.precond
          << ',' << o.ilu_ordering << ',' << o.steps << ',' << o.newton_steps << ','
          << step << ',' << out.updates << ',' << out.newton_seconds << ',' << mean
          << ',' << out.seconds;
  for (double x : out.errors) { summary << ',' << x; }
  summary << ',' << out.drift << ',' << o.check_update << '\n';
  std::cout << form << " k=" << o.k << " " << o.nx << "x" << o.ny
            << " updates=" << out.updates << " mean_newton_s=" << mean
            << " rho_L2=" << out.errors[0] << " drift=" << out.drift << '\n';
  return out;
}
} // namespace euler
