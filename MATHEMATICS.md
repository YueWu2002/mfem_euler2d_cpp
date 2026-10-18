# Discrete equations and metrics

This implementation follows the frozen-HLL construction and paired Gram metric
in `../../GPT/note.pdf`, sections 2.4, 3.1 and 3.2.1. The two sibling C++ projects
are implementation references; they do not define a different quadrature here.

## Rotated gamma=3 simple wave

Let `d=(cos(alpha),sin(alpha))`, `z=d dot (x,y)`, and solve

\[
 \xi+t\sin(\pi\xi)=z \pmod 2,\qquad j=\sin(\pi\xi).
\]

Then `v=1+j/2`, `rho=(1-j/2)/sqrt(3)`, `p=rho^3`,
`U=(rho,rho*v*d_x,rho*v*d_y,rho^3/2+rho*v^2/2)`.
The first shock time is `1/pi`. Each domain translation changes z by 2, so
periodicity is compatible with every supported angle, not only pi/4.
The supported box has positive lengths and requires `0<alpha<pi/2`.

## Target DG equation

Use Qk GL nodal basis functions and exactly `(k+1)^2` tensor GL volume points;
use `k+1` shared GL points on each face. All four components share the scalar
basis. With physical cell Gram matrix `M_K=area*M_ref`, backward Euler is

\[
 r_K^w=M_K(U_K-U_K^n)-\Delta t\,(F(U),\nabla\phi)_K
       +\Delta t\,\langle\widehat F_{K},\phi\rangle_{\partial K}=0.
\]

The code stores `r_K=M_K^{-1}r_K^w`. This is a row scaling of the same DG
method; it must be accounted for when choosing residual norms. The GL
quadrature mass is diagonal and exact on products of Qk basis functions.
The nonlinear flux integrals are the specified collocated quadrature, not an
exact nonlinear integral or an overintegrated variant. Initial projection and
final error evaluation use a separate accurate rule outside the time loop.

At a face with fixed positive-axis normal, define `a=-SL>=0`, `b=SR>=0`.
With square-root-density Roe averages,

\[
 a=\eta\max(0,c_L-v_{n,L},\widetilde c-\widetilde v_n),\quad
 b=\eta\max(0,c_R+v_{n,R},\widetilde c+\widetilde v_n),\qquad\eta=1.1.
\]

Here `c=sqrt(3p/rho)`, and Roe sound speed comes from the averaged total
enthalpy and full velocity. These are exactly the reference implementations'
Einfeldt bounds. Freeze a,b from the preceding time level throughout Newton,
GMRES and line search. The physical flux is still nonlinear in the current U.
The oriented HLL flux is

\[
 \widehat F=\frac{bF_n(U_L)+aF_n(U_R)+ab(U_L-U_R)}{a+b}.
\]

Zero a or b gives the corresponding one-sided physical flux, without a
separate differentiability branch. Require `a+b>0`.

## Equivalent nodal trace construction

The trace space is Pk on each edge, with Lagrange basis at exactly the shared
GL points, so its evaluation matrix is the identity. Set

\[
 \Phi(U_L,U_R)=\frac{aU_L+bU_R+F_n(U_L)-F_n(U_R)}{a+b},
\]
\[
 H_L=F_n(U_L)+a(U_L-\Lambda),\qquad
 H_R=-F_n(U_R)+b(U_R-\Lambda).
\]

These are outward fluxes; the right outward normal is the negative of the
fixed face normal. The stored skeleton row for each node/component is

\[
 g_{f,q}=w_q[(a+b)\Lambda_q-aU_{L,q}-bU_{R,q}-F_n(U_L)+F_n(U_R)].
\]

Physical face length `ell_f` has been divided out of the skeleton equation.
Positive weights give `g=0 <=> Lambda=Phi` at every point; substitution returns
the chosen DG HLL flux exactly. No projection of a nonlinear face flux onto
fewer moments occurs.

Let `r_H` be the mass-inverted cell residual. With frozen coefficients,
`A=dr_H/dU`, `B=dr_H/dLambda`, `C=dg/dU`, `D=dg/dLambda`, where A is cell-local,
D is face-local and diagonal (`w_q(a+b)` at each component). Then

\[
 R_{DG}=r_H-BD^{-1}g,\qquad J_{DG}=A-BD^{-1}C.
\]

Thus the exact cell Newton correction is independent of the current trace,
including off-constraint states. Production reconstructs the trace at each
Newton state: on this manifold `r_H=R_DG`, `g=0`. It uses these identities
without recomputing redundant trace residuals. This is the note's optional
trace reconstruction, not nonlinear elimination of cell states.

HDG advancement **always eliminates cell increments and globally solves S**:

\[
 X_K=A_K^{-1}B_K,\quad y_K=A_K^{-1}r_{H,K},\qquad
 S=D-CA^{-1}B,\quad b_S=-g+CA^{-1}r_H,
\]
\[
 S\delta\Lambda=b_S,\qquad
 \delta U_K=-y_K-X_KR_K\delta\Lambda.
\]

Factor each A_K once by pivoted LU and solve for X_K and y_K with that factor;
never explicitly form A_K^{-1}. Partial action evaluates `D*v-C*(X*v)` without
forming a global S or even the full cell-local C*X. Element/full modes materialize
these products. Global direct modes also receive S only, never the augmented
HDG matrix. DG and HDG share the same DG-based positivity/Armijo acceptance.

The periodic grid has 4*Nx*Ny*(k+1)^2 cell DOFs and 8*Nx*Ny*(k+1) trace DOFs,
so trace/cell DOF ratio is `2/(k+1)`. This alone is not a speedup estimate.

## Dual norm and basis-independent right GMRES

For an **unscaled weak residual** e, the required norm is

\[
 \|e\|_{H^{-1}}=\sqrt{e^T H^{-1}e},
\]

not `sqrt(e^T H e)`. Let H be the fixed cell or skeleton L2 Gram matrix and
`H=L L^T`. The note's normalized system and Jacobi PC are

\[
 \widetilde K=L^{-1}KL^{-T},\quad \widetilde b=L^{-1}b,
 \quad\widetilde P=\operatorname{blockdiag}(\widetilde K_{ii}).
\]

Right-preconditioned GMRES minimizes the original dual residual norm.
Equivalently, operate directly in weak coordinates with inverse-H-weighted
Arnoldi products and the physical block PC. The code uses this equivalent
weighted Arnoldi implementation, with reusable MFEM vectors/matrices.

**Stored-equation weights differ because the equations have been row scaled:**

| System | Stored residual e_s | Weak residual e_w | Arnoldi weight W, giving `e_s^T W e_s=e_w^T H^{-1}e_w` |
|---|---|---|---|
| DG | `M^{-1} e_w` | `M e_s` | `M=area*M_ref` |
| HDG Schur | `ell_f^{-1} e_w` | `ell_f e_s` | `ell_f/w_q` per nodal component |

The HDG Gram is `H_trace=diag(ell_f*w_q)`, counting each global face once.
A physically scaled cell residual is `r_w/dt`; therefore nonlinear convergence
uses `sqrt(r_stored^T M r_stored)/dt/sqrt(domain area)`. The common relative
DG linear-residual criterion is unchanged by division by the same dt.

Using inverse M again on the **already mass-inverted** DG residual would apply
the wrong metric. Dropping ell_f in the skeleton weight is likewise wrong
when vertical and horizontal face lengths differ.

Block Jacobi uses exact diagonal blocks J_KK for DG and **S_ff**, including
local elimination contributions, for HDG; it does not use D_ff alone.
Fixed-order block ILU(0) is also covariant under local block basis changes in
exact arithmetic when the full block sparsity pattern is fixed. MDF can change
the block ordering under basis changes because it depends on numerical block
norms; use `--ilu-ordering none` for invariance comparisons.

For a nonorthogonal local basis transformation T, use `K'=T^T K T`,
`b'=T^T b`, `H'=T^T H T`. H' is generally dense and **must not be replaced by
its diagonal**. The optional verification explicitly Cholesky-whitens this
full transformed Gram, compares physical solutions and restarted-GMRES counts,
and checks both PCs. Production is restricted to the stated GL nodal basis,
where the Gram is diagonal. Basis invariance does not assert identical finite
Krylov iterates between the distinct DG and Schur systems.

## Schur stopping and update validation

With exact local recovery and `epsilon_S=S*deltaLambda-b_S`,

\[
 \epsilon_{DG}=J_{DG}\delta U+R_{DG}=-BD^{-1}\epsilon_S.
\]

Production checks both the true weighted Schur residual and this recovered
cell residual in the common cell metric. If necessary it tightens the Schur
relative tolerance by 100, at most three further continuations from the current correction, reusing
the same assembled matrix and PC factors. All retry costs are
included in the Newton timer. There is a floating-point absolute allowance
proportional to machine epsilon and cell state norms.

`--check-update` independently forms a DG action and solves the DG Newton
equation at the same U, then compares all conservative components of deltaU.
`--verify` also checks analytic Jacobians by finite differences, partial/full
actions and off-manifold traces. These are supplied tests, not results already
obtained. Positivity checks at evaluation points are not a proof of polynomial
positivity everywhere or unconditional Newton convergence.
