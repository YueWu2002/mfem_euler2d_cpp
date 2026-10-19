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
                       const Vector &schur_residual, Vector &result) {
  auto &scaled = s.scaled_scratch;
  for (int f = 0; f < p.nf; ++f) {
#pragma omp simd
    for (int i = 0; i < p.t; ++i) {
      scaled[f*p.t+i] = schur_residual[f*p.t+i]/h.D[f*p.t+i];
    }
  }
  auto &local = s.local_scratch; auto &cell = s.cell_scratch;
  for (int e = 0; e < p.n; ++e) {
    for (int i = 0; i < 4*p.t; ++i) { local[i] = scaled[h.ids[e][i]]; }
    h.B[e].Mult(local,cell);
#pragma omp simd
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
    log << "form,matrix_free,step,stage,newton,time,stage_time,dt,cfl,"
           "residual_physical_rms,damping,"
           "krylov_iterations,linear_relative_residual,dg_linear_relative_residual,"
           "residual_s,validation_s,jacobian_s,condensation_s,assembly_s,setup_s,"
           "solve_s,recovery_s,line_search_s,iteration_s,update_check_relative\n"
        << std::setprecision(17);
  }
  // Storage allocation is outside the time loop; every numerical refresh and
  // factorization is still performed and timed at the current Newton state.
  Vector residual(p.nu), rhs(p.nu), er(p.nu), du(p.nu), trial(p.nu),
      trial_residual(p.nu), step_state(p.nu), stage_base(p.nu);
  std::array<Vector,4> slope;
  for (auto &v : slope) { v.SetSize(p.nu); }
  Hybrid hybrid;
  Condensed condensed;
  if (form == "hdg") { p.PrepareHybrid(hybrid,condensed); }
  std::unique_ptr<ElementOperator> dg_jacobian;
  if (form == "dg" && o.assembly == "partial") { dg_jacobian = p.PartialDG(out.u,false); }
  const int global_size = form == "dg" ? p.nu : p.nl;
  LinearSolver linear_solver(global_size, form == "dg" ? p.m : p.t, o,
                             form == "dg" ? p.cell_metric : p.trace_metric);
  if(o.solver=="gmres" && o.precond=="block-ilu0") {
    linear_solver.Prepare(p.SolverGraph(form=="hdg"),o);
  }
  LinearResult lr, attempt_result;
  lr.x.SetSize(global_size); lr.equation_residual.SetSize(global_size);
  attempt_result.x.SetSize(global_size);
  attempt_result.equation_residual.SetSize(global_size);
  Options linear = o;
  double time = 0, begin = clock_s(), dt_sum = 0., max_cfl = 0.;
  int step = 0;
  // SSP-DIRK3(4,3,3), applied as a purely implicit RK method.  The common
  // diagonal keeps every stage Jacobian in the existing I-a_ii*dt*L' form.
  constexpr double alpha = 0.24169426078821;
  constexpr double beta = alpha/4.;
  constexpr double eta = 0.12915286960590;
  const double ssp_a[4][4] = {
      {alpha, 0., 0., 0.},
      {-alpha, alpha, 0., 0.},
      {0., 1.-alpha, alpha, 0.},
      {beta, eta, .5-beta-eta-alpha, alpha}};
  const double ssp_b[4] = {0., 1./6., 1./6., 2./3.};
  const double ssp_c[4] = {alpha, 0., 1., .5};
  bool final_step = time >= o.end;
  while (o.steps > 0 ? step < o.steps : !final_step) {
    // Beginning-of-step bounds select dt and also serve the first stage.
    step_state = out.u;
    // Initial, accepted Newton, and reconstructed DIRK states are checked at
    // their creation sites, so this time-level state is already admissible.
    p.Freeze(step_state, o.dt, false);
    double rate = p.WaveRate();
    double dt = o.cfl > 0 ? o.cfl / rate : o.dt;
    if (o.steps < 0) {
      final_step = time + dt >= o.end;
      if (final_step) { dt = o.end - time; }
    }
    require(std::isfinite(dt) && dt > 0 && time + dt < 4 / ((p.gas.Gamma()+1)*pi),
            "Step leaves the smooth simple-wave interval T < 4/((gamma+1)*pi)");
    const int stage_count = o.integrator == "be" ? 1 : 4;
    for (int stage = 0; stage < stage_count; ++stage) {
      const double diagonal = o.integrator == "be" ? 1. : ssp_a[stage][stage];
      stage_base = step_state;
      if (o.integrator != "be") {
        for (int j = 0; j < stage; ++j) {
          stage_base.Add(dt*ssp_a[stage][j],slope[j]);
        }
      }
      // Freeze at the stage's initial guess (previous stage solution for
      // stage > 0). Keep these bounds for every Newton and line-search state.
      // Stage zero reuses the beginning-of-step estimate used for CFL.
      if (stage > 0) { p.Freeze(out.u, diagonal*dt, false); }
      p.dt = diagonal*dt;
      p.old = stage_base;
      const double physical_scale = 1 / (p.dt * std::sqrt(p.lx * p.ly));
      const double old_norm = p.CellNorm(stage_base);
      double r0 = -1.;
      bool converged = false, residual_ready = false;
      const int limit = o.newton_steps > 0 ? o.newton_steps : o.max_newton;
      for (int it = 0; it <= limit; ++it) {
        double istart = clock_s();
        if (!residual_ready) {
          p.Residual(out.u,residual);
        }
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
        double dgrel = 0, er_norm = 0, update_check = -1., diagnostic_time = 0.;
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
          // Solve already formed J*du-rhs = J*du+residual for its mandatory
          // true-residual check. Reuse it instead of applying J a second time.
          astart = clock_s();
          er_norm = p.CellNorm(lr.equation_residual);
          vt = clock_s() - astart;
        } else {
          // Reconstruct Lambda=Phi(U) at each Newton state. On this manifold
          // r_H=r_DG and g=0 algebraically, so avoid redundant flux/residual passes.
          // Only S enters the global solver, and only local backsolves recover du.
          auto &h = hybrid; auto &s = condensed;
          // Streaming combines Jacobian assembly and exact local condensation.
          // Both are included in condensation_s and hence in iteration_s.
          p.StreamCondense(out.u,residual,h,s);
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
            astart = clock_s();
            p.Recover(s,attempt_result.x,du);
            rt += clock_s() - astart;
            astart = clock_s();
            RecoveredResidual(p,h,s,attempt_result.equation_residual,er);
            er_norm = p.CellNorm(er);
            valid = er_norm <= target;
            vt += clock_s() - astart;
            if (valid) { break; }
            require(o.solver == "gmres", "Direct Schur solve fails DG residual criterion");
            linear.linear_tol *= .01;
          }
          require(valid, "Recovered DG residual criterion failed after Schur refinement");
        }
        dgrel = residual_norm > 0 ? er_norm / residual_norm : 0.;
        require(er_norm <= target, "DG Newton linear residual too large");
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
          const double stage_time = time +
              (o.integrator == "be" ? 1. : ssp_c[stage])*dt;
          log << form << ',' << o.matrix_free << ',' << step << ',' << stage << ','
              << it << ',' << time + dt << ',' << stage_time << ',' << dt << ','
              << rate*dt << ',' << rms << ','
              << damping << ',' << lr.count << ',' << lr.residual << ',' << dgrel
              << ',' << residual_time << ',' << vt << ',' << jt << ',' << ct << ','
              << lr.assembly << ',' << lr.setup << ',' << lr.solve << ',' << rt << ','
              << lst << ',' << iteration_time << ',' << update_check << '\n';
        }
      }
      require(converged || o.newton_steps > 0, "Newton failed to converge");
      slope[stage] = out.u;
      slope[stage] -= stage_base;
      slope[stage] *= 1./p.dt;
    }
    if (o.integrator != "be") {
      out.u = step_state;
      for (int stage = 0; stage < 4; ++stage) {
        out.u.Add(dt*ssp_b[stage],slope[stage]);
      }
      require(p.Admissible(out.u), "SSP-DIRK3 step produced an inadmissible state");
    }
    time = o.steps < 0 && final_step ? o.end : time + dt;
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
  summary << "form,nx,ny,k,gamma,integrator,stages,T,dt,cfl,"
             "actual_max_cfl,mean_dt,"
             "volume_dofs,trace_dofs,matrix_free,assembly,solver,preconditioner,"
             "ilu_ordering,requested_steps,fixed_newton_steps,time_steps,"
             "newton_updates,newton_s,mean_newton_s,time_loop_s,rho_L2,rhou_L2,"
             "rhov_L2,energy_L2,max_conservation_drift,check_update\n"
          << std::setprecision(17);
  const double mean = out.updates ? out.newton_seconds / out.updates : 0.;
  summary << form << ',' << o.nx << ',' << o.ny << ',' << o.k << ','
          << p.gas.Gamma() << ','
          << o.integrator << ',' << (o.integrator == "be" ? 1 : 4) << ','
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
