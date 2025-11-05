# Flavoured Leptogenesis: symbolic + numeric framework
import numpy as np
import matplotlib.pyplot as plt
from scipy import special, integrate
import sympy as sp
sp.init_printing(use_unicode=True)
# MODE: 'symbolic' أو 'numeric'
# MODE = 'numeric'

Mpl_red = 2.435e18       # reduced Planck mass [GeV]
gstar = 106.75
gstar_s = 106.75
vEW = 246.0              # Higgs VEV [GeV] (standard)
gN = 2.0                 # degrees of freedom for Majorana N

z = sp.symbols('z', positive=True)        # z = M1/T (dimensionless)
M1_s, M2_s, M3_s = sp.symbols('M1 M2 M3', positive=True)
# Yukawa matrix symbols Y_{alpha i} (alpha = e,mu,tau indexed 0,1,2)
Y = sp.Matrix([[sp.symbols(f'Y_e{i+1}') for i in range(3)],
               [sp.symbols(f'Y_mu{i+1}') for i in range(3)],
               [sp.symbols(f'Y_tau{i+1}') for i in range(3)]])
# flavoured CP asymmetries eps_{i,alpha}
eps_sym = sp.Matrix([[sp.symbols(f'eps_{i+1}_{a}') for a in range(3)] for i in range(3)])

# helper symbolic functions (leave Bessel as symbolic functions)
K1 = sp.Function('K1')  # K1(x)
K2 = sp.Function('K2')  # K2(x)

# equilibrium yield (general) for Ni (symbolic expression)
# Y_N^eq(z; Mi) = (45 g_N / (4 pi^4 g*_s)) * ( (Mi/M1)^2 * z^2 * K2((Mi/M1)*z) )
pref_Yeq = (45 * gN) / (4 * sp.pi**4 * gstar_s)
r1, r2, r3 = sp.symbols('r1 r2 r3', positive=True)  # r_i = M_i / M1
Yeq_sym_i = lambda ri: pref_Yeq * ri**2 * z**2 * sp.Function('K2')(ri*z)

# Decay parameter D_i (symbolic form) — using argument x_i = r_i * z
# D_i(z) = (Gamma_i / H(T=Mi/z)) * z * K1(x_i)/K2(x_i)
# We'll leave Gamma_i symbolic as Gamma_i (depends on Yukawas)
Gamma1_s, Gamma2_s, Gamma3_s = sp.symbols('Gamma1 Gamma2 Gamma3', positive=True)
H_sym = lambda T: sp.sqrt(8*sp.pi**3 * gstar / 90) * T**2 / Mpl_red  # symbolic form (T remains symbol)
# Represent H at T = Mi/z as H_i(z) ~ sqrt(...) * (Mi/z)^2 / Mpl -> we'll show pattern
# For clarity keep D_i symbolic:
ri = sp.symbols('ri', positive=True)
D_sym = sp.Function('D')  # represent D_i symbolically

# Flavoured Boltzmann equations (symbolic form)
# For each i (heavy neutrino): dY_Ni/dz = - D_i(z) * (Y_Ni - Y_Ni^eq)
# For each flavour alpha: dY_Delta_alpha/dz = - sum_i eps_{i,alpha} D_i (Y_Ni - Y_Ni^eq) - W_alpha(z) Y_Delta_alpha
Y_N1_sym, Y_N2_sym, Y_N3_sym = sp.symbols('Y_N1 Y_N2 Y_N3')
YDL1, YDL2, YDL3 = sp.symbols('Y_DL_e Y_DL_mu Y_DL_tau')

Y_N1 = sp.Function('Y_N1')(z)
Y_N2 = sp.Function('Y_N2')(z)
Y_N3 = sp.Function('Y_N3')(z)

dY_N1 = -sp.Function('D1')(z) * (Y_N1 - Yeq_sym_i(r1))
dY_N2 = -sp.Function('D2')(z) * (Y_N2 - Yeq_sym_i(r2))
dY_N3 = -sp.Function('D3')(z) * (Y_N3 - Yeq_sym_i(r3))

# flavoured dYDL (symbolic, show summation compactly)
# dY_DL_alpha = - sum_i eps_{i,alpha} D_i (Y_Ni - Y_Ni^eq) - W_alpha(z) * Y_DL_alpha
dYDL_alpha = lambda alpha: -sum(
    eps_sym[i, alpha] * sp.Function(f'D{i+1}')(z) *
    (sp.Function(f'Y_N{i+1}')(z) - Yeq_sym_i(sp.symbols(f'r{i+1}')))
    for i in range(3)
) - sp.Function(f'W_{["e","mu","tau"][alpha]}')(z) * sp.Function(f'Y_DL_{["e","mu","tau"][alpha]}')(z)
# If MODE == 'symbolic' -> print and return clean symbolic forms
if MODE == 'symbolic':
    print("\n=== SYMBOLIC REPRESENTATION (summary) ===\n")
    print("Equilibrium yield (symbolic) for Ni:")
    print("Y_N^eq_i(z) = pref * r_i^2 * z^2 * K2(r_i z)   where pref = 45 g_N /(4 pi^4 g*_s)")
    sp.pprint(pref_Yeq)
    print("\nExample: Y_N^eq (for r_i) =")
    sp.pprint(Yeq_sym_i(ri))
    print("\nSymbolic Boltzmann (heavy neutrinos):")
    print("dY_N1/dz = - D1(z) * (Y_N1(z) - Y_N1^eq(z))")
    sp.pprint(dY_N1)
    print("\nSymbolic flavoured Boltzmann (generic alpha):")
    sp.pprint(dYDL_alpha(0))   # print for alpha=0 (electron flavour)
    print("\nNote: D_i(z) left as symbolic function (depends on Gamma_i and H(T=Mi/z)).")
    print("You can substitute explicit Gamma_i expressions (e.g., Gamma_i = (Y^†Y)_{ii} Mi/(8 pi))")
    print("\nSphaleron conversion (symbolic): Y_B = c_s * sum_alpha Y_DL_alpha")
    c_s = sp.Rational(-28,79)
    sp.pprint(c_s)

    print("\n\n=== END SYMBOLIC SUMMARY ===\n")
    # End of symbolic mode
    raise SystemExit("Symbolic mode: finished printing symbolic framework. Set MODE='numeric' to run numeric solver.")
# 4) NUMERIC MODE: substitue values and solve
# If here, MODE=='numeric'
# (Note: careful with epsilon magnitudes — we print a warning if any |eps|>1)
# Default numeric parameters (you can edit these benchmark values)
M1 = 1e10
M2 = 5e11
M3 = 1e13
M_list = np.array([M1, M2, M3])
r = M_list / M1

# Yukawa (example approximate values — editable)
Yukawas = np.array([1e-3, 5e-4, 2e-4])   # these are diagonal effective couplings used for Gamma estimate

# Decay widths (simple Yukawa formula Gamma = y^2 M /(8 pi))
Gamma = (Yukawas**2 * M_list) / (8.0 * np.pi)

# CP asymmetries (from text examples) -- CAUTION: some paper values can be >1; check physicality.
eps = np.array([-0.066, -1.012, 0.0005])   # default from provided benchmarks (user provided)
if np.any(np.abs(eps) > 1.0):
    print("WARNING: Some |eps| > 1. Ensure these are per-decay asymmetries or rescale appropriately for physicality.")

# H(T) and s(T) numeric functions
def H_of_T_num(T, gstar_local=gstar, Mpl=Mpl_red):
    return np.sqrt(8 * np.pi**3 * gstar_local / 90.0) * T**2 / Mpl

def s_of_T_num(T, gstar_s_local=gstar_s):
    return (2.0 * np.pi**2 / 45.0) * gstar_s_local * T**3

# Equilibrium yield numeric
def Y_N_eq_num(z, Mi, Mref=M1):
    T = Mref / z
    x = Mi / T  # = (Mi/Mref) * z
    K2 = special.kv(2, x)
    n_eq = gN / (2.0 * np.pi**2) * Mi**2 * T * K2
    s = s_of_T_num(T)
    return n_eq / s

# D_i(z) numeric (corrected: use H evaluated at T = Mi/z and Bessel args = x_i = (Mi/M1)*z)
def H_at_z(z, Mi):
    T = Mi / z
    return H_of_T_num(T)

def D_i_num(z, i):
    xi = r[i] * z
    K1 = special.kv(1, xi)
    K2 = special.kv(2, xi)
    # Avoid division by zero for tiny xi
    if K2 == 0:
        return 0.0
    return (Gamma[i] / H_at_z(z, M_list[i])) * z * K1 / K2

# Simple washout terms
# W_ID (inverse-decay) standard approximations
def W_ID_num(z, i):
    # using often-used form W_ID = 0.5 * D_i(z) * Y_N_eq(z) / Y_l_eq (we approximate Y_l_eq ~ constant ~ 3/4 * 45/(2*pi^4*gstar_s) ???)
    # For simplicity, use W_ID = 0.5 * D_i * (Y_N_eq / Y_ref) with Y_ref = 1e-3 typical scale (user can refine)
    Yref = 1e-3
    return 0.5 * D_i_num(z, i) * (Y_N_eq_num(z, M_list[i]) / Yref)

# ΔL=2 washout rough approximation (uses sum m_nu^2)
mnu_vals = np.array([0.0, 8.6e-3, 0.05])  # sample light neutrino masses [eV] (editable)
mbar2 = np.sum(mnu_vals**2)
# conversion eV^2 -> GeV^2 : 1 eV = 1e-9 GeV
mbar2_GeV2 = mbar2 * (1e-9)**2

def W_DeltaL2_num(z, i):
    Mi = M_list[i]
    # approximate pre-factor from earlier derivation
    return 4.49e-3 * (Mi * Mpl_red) / (gstar**1.5 * vEW**4) * (mbar2_GeV2) / (z**2)

# Boltzmann RHS numeric (3 Ni + 3 flavours tracked as vector of length 6: Y_N1,Y_N2,Y_N3, Y_DL_e,mu,tau)
def boltzmann_numeric(z, y):
    YN = y[0:3]
    YDL = y[3:6]
    dY = np.zeros(6)
    # Heavy neutrino equations
    for i in range(3):
        Di = D_i_num(z, i)
        Yeqi = Y_N_eq_num(z, M_list[i])
        dY[i] = - Di * (YN[i] - Yeqi)
    # flavour equations (flavour-dependent washout included: inverse decay + DeltaL2)
    for alpha in range(3):
        prod = 0.0
        wash = 0.0
        for i in range(3):
            Di = D_i_num(z, i)
            Yeqi = Y_N_eq_num(z, M_list[i])
            prod += eps[i] * Di * (YN[i] - Yeqi)   # NOTE: here eps is flavour-averaged; extend eps->matrix for true flavoured
            wash += W_ID_num(z, i) + W_DeltaL2_num(z, i)
        dY[3+alpha] = - prod - wash * YDL[alpha]
    return dY

# Initial conditions (thermal abundance for N's and zero lepton asymmetry)
z0 = 0.1
zmax = 200.0
z_vals_num = np.logspace(np.log10(z0), np.log10(zmax), 300)

Y0 = np.zeros(6)
# start N densities at equilibrium at z0
for i in range(3):
    Y0[i] = Y_N_eq_num(z0, M_list[i])
Y0[3:6] = 0.0

# Solve with stiff solver BDF for speed
sol = integrate.solve_ivp(boltzmann_numeric, [z_vals_num[0], z_vals_num[-1]], Y0,
                          t_eval=z_vals_num, method='BDF', rtol=1e-6, atol=1e-12, max_step=1.0)

if not sol.success:
    print("Warning: ODE solver failed:", sol.message)

# Extract results
Y_Ns = sol.y[0:3,:]
Y_DL_flavours = sol.y[3:6,:]
Y_DL_total = np.sum(Y_DL_flavours, axis=0)

# Sphaleron conversion: Y_B = c_s * Y_{B-L} (approx c_s = -28/79)
c_s = -28.0/79.0
Y_B = c_s * Y_DL_total

# Final values
Y_DL_final = Y_DL_total[-1]
Y_B_final = Y_B[-1]

print("\n=== Numeric run summary ===")
print(f"Final lepton asymmetry (sum flavours) Y_ΔL = {Y_DL_final:.3e}")
print(f"Converted baryon asymmetry Y_B = {Y_B_final:.3e}")
print("Note: check magnitude & sign; physical observed Y_B ~ +8.7e-11 (positive)")

# Plotting
plt.figure(figsize=(8,6))
for i in range(3):
    plt.loglog(z_vals_num, Y_Ns[i,:], label=f'Y_N{i+1}')
plt.loglog(z_vals_num, np.abs(Y_DL_total), linestyle='--', label='|Y_ΔL| (sum flavours)')
plt.xlabel('z = M1/T')
plt.ylabel('Yield')
plt.title('Numeric evolution (example benchmarks)')
plt.legend()
plt.grid(which='both', ls=':', lw=0.6)
plt.show()



