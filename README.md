# deterministic-slice-time-harness

**The capstone of the shadow-resolution arc: 0-loss slice-time extrapolation, or `Held`.** Pure Rust, zero deps, `json=0`, no JSON no Node.

The falsifiable, honest form of *"retrospectively and prospectively predict the next slices at 0 loss, independent of Schrodinger."* Stated as a **deterministic-simulation theorem, not physics.**

## The theorem (honest tags)
- **CANON/DESIGN** — In a deterministic metatagged simulated universe, if multi-cylinder / Q-PRISM / tomography reconstructs the full **state** of a bounded slice **and** the transition rule is known, then past and future slices are **0-loss computable** — they carry no new entropy relative to `(state, rule)`. Return exact bytes or `Held`.
- **MEASURED** (`cargo test`) — recovery half (enough shadows recover exactly, too few `Held`) + deterministic forward/backward extrapolation byte-identical + **new entropy DETECTED** (a truth that gained entropy ≠ the determined prediction).
- **UNVERIFIED** — live cross-fabric Hilbra run; real-world quantum prediction; any claim that genuinely-new entropy is predicted losslessly (Shannon/Fano/DPI forbid it).

## The 6 steps
1. Metatag particles in a bounded 2D torus grid (`SimState`).
2. Emit multiple shadow projections per slice (`project_shadows` — residues on K coprime cylinders = tomographic angles).
3. Reconstruct the slice via multi-cylinder / Q-PRISM (`reconstruct` — inverse Radon / CRT).
4. Apply the **deterministic, reversible** rule forward and backward (`SimState::at(±dt)` — constant-velocity torus drift).
5. Compare predicted slices byte-for-byte against generated truth (`run_deterministic`).
6. Return `Held` when state/rule/projection capacity is insufficient; **detect** (never manufacture) when the truth gained new entropy (`verify_or_detect_entropy`).

## Why "independent of Schrodinger"
A simulation has **no Heisenberg** (state is read exactly from memory) and **no Born collapse** (the rule is deterministic) — it is **Laplacian because it was never quantum**. That is the entire content; it is not a claim about real physics. For a real/stochastic universe only the *determined* component is 0-loss.

## The arc it closes
papers (shadow wall) → capstone (recovery, not inversion) → Path-1 (recall) → Path-2 (no store) → 3D Q-PRISM (any slice) → **deterministic slice-time (any *when*)**.

## License
MIT OR Apache-2.0.
