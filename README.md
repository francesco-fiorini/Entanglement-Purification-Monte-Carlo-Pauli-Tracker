# Entanglement Purification – Monte Carlo Pauli Tracker

Fast Monte-Carlo simulator for entanglement purification protocols over a lossy + depolarizing optical fiber.

## Table of Contents
- [What id does](#whatitdoes)
- [Protocols Included](#protocols)
- [Results](#results)
- [How to run](#howtorun)
- [Contributing](#contributing)
- [License](#license)

### What it does
- Models a realistic optical fiber (0–30 km) with attenuation (0.2 dB/km) + cumulative depolarizing noise (1 % per 10 km).
- Implements 5 standard recurrence purification protocols using **Pauli error tracking** (no density matrices → very fast).
- Runs 20 000 shots per distance point → accurate fidelity and success probability curves.

### Protocols included
| Protocol   | Pairs used | Purpose                     |
|------------|------------|-----------------------------|
| No Purif   | 1          | Raw link                    |
| Ss-SpX     | 2          | Single-selection bit-flip   |
| Ss-SpZ     | 2          | Single-selection phase-flip |
| Ds-Sp      | 3          | Double-selection single-purif |
| Ss-Dp      | 3          | Single-selection double-purif (X+Z) |
| Ds-Dp      | 5          | Deep double-selection double-purif |

### Results
Generates two clean 4:3 figures:
- `purification_fidelity_4x3.png` → Fidelity vs distance (with F=0.5 threshold)
- `purification_success_4x3.png` → Success probability vs distance

### How to run
```bash
pip install numpy matplotlib
python purification_sim.py

## Contributing
Contributions are welcome! Please contact francesco.fiorini@phd.unipi.it for suggested changes.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
