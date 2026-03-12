"""
Author: Francesco Fiorini, francesco.fiorini@phd.unipi.it
=============================================================================
Entanglement Purification – Monte-Carlo Simulation (Pauli Tracking)
Protocols : No-Purif, Ss-SpX, Ss-SpZ, Ds-Sp, Ss-Dp, Ds-Dp
Noise     : FIBER ONLY (Attenuation + Depolarizing over full distance d)
Gates     : IDEAL 
Method    : Monte Carlo (Repeated empirical shots tracking Pauli errors)
=============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# PARAMETERS
# ==============================================================================
DISTANCES_KM     = np.arange(0, 31, 2)   # Total fiber lengths [km]
ALPHA_DB_PER_KM  = 0.2                   # Attenuation [dB/km]
DEPOL_FIBER_RATE = 0.01                  # Depolarizing prob per 10 km
N_MONTE_CARLO    = 50000                 # MC shots per distance point
SEED             = 42

np.random.seed(SEED)

# ==============================================================================
# 1. CHANNEL MODEL
# ==============================================================================
def initial_fidelity(d_km: float) -> float:
    """ Computes Werner State fidelity passing through a noisy fiber channel. """
    # 1. Photon Loss (Transmittance)
    eta = 10 ** (-ALPHA_DB_PER_KM * d_km / 10)
    
    # 2. Depolarizing noise in the fiber (Exponential decay with distance)
    depol_survival = max(0.0, (1.0 - DEPOL_FIBER_RATE) ** (d_km / 10.0))
    
    # 3. Combined Fidelity mapping to Werner State
    F_eff = eta * depol_survival + (1 - eta) * 0.25  
    return float(np.clip(F_eff, 0.0, 1.0))

# ==============================================================================
# 2. PAULI TRACKING ENGINE (Ideal Gates)
# ==============================================================================
_mul_table = {
    ('I','I'):'I', ('I','X'):'X', ('I','Y'):'Y', ('I','Z'):'Z',
    ('X','I'):'X', ('X','X'):'I', ('X','Y'):'Z', ('X','Z'):'Y',
    ('Y','I'):'Y', ('Y','X'):'Z', ('Y','Y'):'I', ('Y','Z'):'X',
    ('Z','I'):'Z', ('Z','X'):'Y', ('Z','Y'):'X', ('Z','Z'):'I',
}
def pmul(a: str, b: str) -> str: return _mul_table[(a, b)]

_bell_from_pauli = {'I': 'Phi+', 'X': 'Psi+', 'Z': 'Phi-', 'Y': 'Psi-'}
def get_bell_state(pair: tuple) -> str: return _bell_from_pauli[pmul(*pair)]
def pair_fidelity(pair: tuple) -> float: return 1.0 if get_bell_state(pair) == 'Phi+' else 0.0  

def measure_z_parity(pair: tuple) -> bool: return get_bell_state(pair) in ('Phi+', 'Phi-')
def measure_x_parity(pair: tuple) -> bool: return get_bell_state(pair) in ('Phi+', 'Psi+')

def sample_pair(F: float) -> tuple:
    r = np.random.random()
    if   r < F:              return ('I', 'I')
    elif r < F + (1-F)/3:    return ('Z', 'I')
    elif r < F + 2*(1-F)/3:  return ('X', 'I')
    else:                    return ('Y', 'I')

def cnot_propagate(ctrl: tuple, tgt: tuple) -> tuple:
    eAc, eBc = ctrl
    eAt, eBt = tgt

    nAc, nAt = eAc, eAt
    if eAc in ('X', 'Y'): nAt = pmul(eAt, 'X')
    if eAt in ('Z', 'Y'): nAc = pmul(eAc, 'Z')

    nBc, nBt = eBc, eBt
    if eBc in ('X', 'Y'): nBt = pmul(eBt, 'X')
    if eBt in ('Z', 'Y'): nBc = pmul(eBc, 'Z')

    return (nAc, nBc), (nAt, nBt)

# ==============================================================================
# 3. PROTOCOLS
# ==============================================================================
def protocol_ss_spX(pairs: list):
    ctrl, tgt = pairs[0], pairs[1]
    ctrl, tgt = cnot_propagate(ctrl, tgt)
    if measure_z_parity(tgt): return ctrl
    return None

def protocol_ss_spZ(pairs: list):
    ctrl, tgt = pairs[0], pairs[1]
    # CNOT reversed for Z-error purification
    tgt, ctrl = cnot_propagate(tgt, ctrl)
    if measure_x_parity(tgt): return ctrl
    return None

def protocol_ds_sp(pairs: list):
    ctrl, anc1, anc2 = pairs[0], pairs[1], pairs[2]
    # Primary to ancilla 1 (X purification)
    ctrl, anc1 = cnot_propagate(ctrl, anc1)
    # Ancilla 2 to ancilla 1 (Double selection verification)
    anc2, anc1 = cnot_propagate(anc2, anc1)
    

    if measure_z_parity(anc1) and measure_x_parity(anc2): return ctrl
    return None

def protocol_ss_dp(pairs: list):
    ctrl, anc_x, anc_z = pairs[0], pairs[1], pairs[2]
    # Phase 1: X purification
    ctrl, anc_x = cnot_propagate(ctrl, anc_x)
    # Phase 2: Z purification (reversed CNOT)
    anc_z, ctrl = cnot_propagate(anc_z, ctrl)
    
    if measure_z_parity(anc_x) and measure_x_parity(anc_z): return ctrl
    return None

def protocol_ds_dp(pairs: list):
    ctrl, cd, ef, gh, ij = pairs[0], pairs[1], pairs[2], pairs[3], pairs[4]
    
    # --- X-error Double Selection Block ---
    ctrl, cd = cnot_propagate(ctrl, cd)     # X purification
    ef, cd   = cnot_propagate(ef, cd)       # Verification check on X
    
    # --- Z-error Double Selection Block ---
    gh, ctrl = cnot_propagate(gh, ctrl)     # Reversed: Z purification
    
    gh, ij   = cnot_propagate(gh, ij)       # Reversed verification check on Z


    if (measure_z_parity(cd) and measure_x_parity(ef) and
        measure_x_parity(gh) and measure_z_parity(ij)): 
        return ctrl
    return None

# ==============================================================================
# 4. SWEEP EXECUTION
# ==============================================================================
protocols = {
    'No Purif': (None,             1),
    'Ss-SpX':   (protocol_ss_spX,  2),
    'Ss-SpZ':   (protocol_ss_spZ,  2),
    'Ds-Sp':    (protocol_ds_sp,   3),
    'Ss-Dp':    (protocol_ss_dp,   3),
    'Ds-Dp':    (protocol_ds_dp,   5),
}

results = {name: {'fidelity': [], 'p_success': []} for name in protocols}

for d in DISTANCES_KM:
    F0 = initial_fidelity(d)
    for name, (func, n_pairs) in protocols.items():
        fids, n_success = [], 0
        for _ in range(N_MONTE_CARLO):
            raw_pairs = [sample_pair(F0) for _ in range(n_pairs)]
            if func is None:
                fids.append(pair_fidelity(raw_pairs[0]))
                n_success += 1
            else:
                out = func(raw_pairs)
                if out is not None:
                    fids.append(pair_fidelity(out))
                    n_success += 1
        results[name]['fidelity'].append(float(np.mean(fids)) if fids else 0.0)
        results[name]['p_success'].append(n_success / N_MONTE_CARLO)

# ==============================================================================
# 5. PLOTTING (Side-by-Side 14x6 format)
# ==============================================================================
colors  = ['#333333', '#1f77b4', '#9467bd', '#ff7f0e', '#2ca02c', '#d62728']
markers = ['s', '^', 'v', 'o', 'D', 'P']
lstyles = ['dotted', 'solid', 'solid', 'dashed', 'dashdot', (0, (5,2))]

plt.rcParams.update({'font.size': 14, 'axes.labelsize': 15, 'axes.titlesize': 16})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Fidelity
for i, name in enumerate(protocols):
    ax1.plot(DISTANCES_KM, results[name]['fidelity'], label=name, color=colors[i],
             linestyle=lstyles[i], marker=markers[i], markersize=7, linewidth=2.0)
    
ax1.axhline(0.5, color='red', linestyle='--', linewidth=1.5, label='Threshold Fidelity = 0.5')
ax1.set_xlabel('Fiber length (km)')
ax1.set_ylabel('Mean Output Fidelity')
ax1.set_ylim(0.35, 1.05)
ax1.grid(True, alpha=0.3)
ax1.legend(ncol=2, loc='upper right')

# Plot 2: Success Probability
for i, name in enumerate(protocols):
    if name == 'No Purif': continue
    ax2.plot(DISTANCES_KM, results[name]['p_success'], label=name, color=colors[i],
             linestyle=lstyles[i], marker=markers[i], markersize=7, linewidth=2.0)

ax2.set_xlabel('Fiber length (km)')
ax2.set_ylabel('Parity Success Probability')
ax2.set_ylim(0, 1.05)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper right')


ax1.text(-0.1, 1.05, 'a)', transform=ax1.transAxes, size=20, weight='bold')
ax2.text(-0.1, 1.05, 'b)', transform=ax2.transAxes, size=20, weight='bold')

plt.tight_layout()
plt.savefig('FidelityandProbSuccess.png', dpi=300, bbox_inches='tight')
plt.show()
