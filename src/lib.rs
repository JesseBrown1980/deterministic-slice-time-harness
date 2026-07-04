//! # deterministic-slice-time-harness — the capstone: slice-time extrapolation at 0 loss, or Held
//!
//! The falsifiable form of "retrospectively and prospectively predict the next slices at 0 loss,
//! independent of Schrodinger." Stated HONESTLY as a deterministic-simulation theorem, not physics:
//!
//! CANON/DESIGN: in a deterministic metatagged simulated universe, if multi-cylinder / Q-PRISM /
//!   tomography reconstructs the full STATE of a bounded slice AND the transition rule is known, then
//!   past and future slices are 0-loss COMPUTABLE - they carry no new entropy relative to (state,
//!   rule). The harness returns exact bytes or Held.
//! MEASURED (this crate, `cargo test`): the recovery half (enough shadows recover exactly, too few
//!   Hold) + the deterministic forward/backward extrapolation is byte-identical + NEW ENTROPY is
//!   DETECTED (the determined prediction != a truth that gained entropy).
//! UNVERIFIED: live cross-fabric Hilbra run; real-world quantum prediction; ANY claim that genuinely
//!   new entropy can be predicted losslessly (Shannon/Fano/DPI forbid it).
//!
//! Why "independent of Schrodinger": a simulation has no Heisenberg (state is read exactly from
//! memory) and no Born collapse (rules are deterministic) - it is Laplacian because it was never
//! quantum. That is the whole content; it is not a claim about real physics.
//!
//! Pure Rust, ZERO deps, HBI/HBP json=0, no JSON no Node.

// ---------------------------------------------------------------- multi-cylinder recovery (vendored)
pub const CYLINDERS: [u64; 4] = [33_554_467, 33_554_393, 33_554_213, 33_550_609]; // pairwise-coprime moduli
const BLOCK: usize = 6; // 48-bit blocks; any 2 of these ~2^25 moduli exceed 2^48

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Held {
    /// too few shadows to reconstruct the state (Shannon roof): product(moduli) < 2^(8*block)
    InsufficientShadows,
    /// the reconstructed bytes are not a valid serialized state
    BadState,
    /// the deterministic prediction diverged from truth: the truth gained entropy the rule can't fix
    PredictionDivergedNewEntropy,
}

fn blocks(d: &[u8]) -> Vec<u128> {
    d.chunks(BLOCK).map(|c| {
        let mut v = 0u128;
        for &b in c { v = (v << 8) | b as u128; }
        if c.len() < BLOCK { v <<= 8 * (BLOCK - c.len()) as u32; }
        v
    }).collect()
}
fn mod_inv(a: u128, m: u128) -> Option<u128> {
    fn e(a: i128, b: i128) -> (i128, i128, i128) { if a == 0 { (b, 0, 1) } else { let (g, x, y) = e(b % a, a); (g, y - (b / a) * x, x) } }
    if m == 0 { return None; }
    let (g, x, _) = e((a % m) as i128, m as i128);
    if g != 1 { return None; }
    Some((((x % m as i128) + m as i128) % m as i128) as u128)
}
/// The K shadow projections of a slice: residues of each block on each cylinder (the tomographic
/// projections at K "angles"). Each shadow alone is lossy; enough of them reconstruct the slice.
pub fn project_shadows(slice: &[u8]) -> Vec<Vec<u64>> {
    CYLINDERS.iter().map(|&p| blocks(slice).iter().map(|&b| (b % p as u128) as u64).collect()).collect()
}
/// Reconstruct the exact slice from a SUBSET of shadows (the inverse Radon / CRT). `Held` if the
/// subset's joint modulus doesn't cover a block (under-sampled angles -> ill-posed).
pub fn reconstruct(shadows: &[Vec<u64>], subset: &[usize], orig_len: usize) -> Result<Vec<u8>, Held> {
    if subset.is_empty() { return Err(Held::InsufficientShadows); }
    let range = 1u128 << (8 * BLOCK as u32);
    let nb = shadows[subset[0]].len();
    let mut out = Vec::with_capacity(nb * BLOCK);
    for bi in 0..nb {
        let (mut r, mut m) = (0u128, 1u128);
        for &si in subset {
            let p = CYLINDERS[si] as u128;
            let s = shadows[si][bi] as u128;
            let inv = mod_inv(m % p, p).ok_or(Held::InsufficientShadows)?;
            let diff = (((s as i128 - r as i128) % p as i128) + p as i128) % p as i128;
            r += m * ((diff as u128 * inv) % p);
            m *= p;
            if m >= range { break; }
        }
        if m < range { return Err(Held::InsufficientShadows); }
        for i in (0..BLOCK).rev() { out.push(((r >> (8 * i as u32)) & 0xFF) as u8); }
    }
    out.truncate(orig_len);
    Ok(out)
}

// ---------------------------------------------------------------- metatagged deterministic sim
/// A metatagged particle on the bounded torus grid.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Particle { pub x: u16, pub y: u16, pub vx: i16, pub vy: i16, pub tag: u32 }

/// The simulated-universe STATE = a bounded grid of metatagged particles. Its serialization is the
/// "slice"; the transition rule is constant-velocity torus drift (deterministic + REVERSIBLE, so
/// retrodiction and extrapolation are the same operation at +/-t).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SimState { pub w: u16, pub h: u16, pub particles: Vec<Particle> }

impl SimState {
    /// Canonical serialization (particles sorted by tag) -> the slice bytes.
    pub fn serialize(&self) -> Vec<u8> {
        let mut ps = self.particles.clone();
        ps.sort_by_key(|p| p.tag);
        let mut out = Vec::with_capacity(6 + ps.len() * 12);
        out.extend_from_slice(&self.w.to_be_bytes());
        out.extend_from_slice(&self.h.to_be_bytes());
        out.extend_from_slice(&(ps.len() as u16).to_be_bytes());
        for p in &ps {
            out.extend_from_slice(&p.x.to_be_bytes());
            out.extend_from_slice(&p.y.to_be_bytes());
            out.extend_from_slice(&p.vx.to_be_bytes());
            out.extend_from_slice(&p.vy.to_be_bytes());
            out.extend_from_slice(&p.tag.to_be_bytes());
        }
        out
    }
    pub fn deserialize(bytes: &[u8]) -> Option<SimState> {
        if bytes.len() < 6 { return None; }
        let w = u16::from_be_bytes([bytes[0], bytes[1]]);
        let h = u16::from_be_bytes([bytes[2], bytes[3]]);
        let n = u16::from_be_bytes([bytes[4], bytes[5]]) as usize;
        if w == 0 || h == 0 { return None; }
        let mut particles = Vec::with_capacity(n);
        for i in 0..n {
            let o = 6 + i * 12;
            if o + 12 > bytes.len() { return None; }
            particles.push(Particle {
                x: u16::from_be_bytes([bytes[o], bytes[o + 1]]),
                y: u16::from_be_bytes([bytes[o + 2], bytes[o + 3]]),
                vx: i16::from_be_bytes([bytes[o + 4], bytes[o + 5]]),
                vy: i16::from_be_bytes([bytes[o + 6], bytes[o + 7]]),
                tag: u32::from_be_bytes([bytes[o + 8], bytes[o + 9], bytes[o + 10], bytes[o + 11]]),
            });
        }
        Some(SimState { w, h, particles })
    }
    /// Deterministic, reversible transition by `dt` steps (constant-velocity torus drift).
    pub fn at(&self, dt: i64) -> SimState {
        let particles = self.particles.iter().map(|p| Particle {
            x: (p.x as i64 + dt * p.vx as i64).rem_euclid(self.w as i64) as u16,
            y: (p.y as i64 + dt * p.vy as i64).rem_euclid(self.h as i64) as u16,
            vx: p.vx, vy: p.vy, tag: p.tag,
        }).collect();
        SimState { w: self.w, h: self.h, particles }
    }
    pub fn step_forward(&self) -> SimState { self.at(1) }
    pub fn step_backward(&self) -> SimState { self.at(-1) }
}

// ---------------------------------------------------------------- the slice-time harness
#[derive(Debug, Clone)]
pub struct SliceTimeReport {
    pub reconstructed_ok: bool,
    /// (dt, byte_identical) for each forward step +1..=horizon
    pub forward_exact: Vec<(i64, bool)>,
    /// (dt, byte_identical) for each backward step -1..=-horizon
    pub backward_exact: Vec<(i64, bool)>,
}
impl SliceTimeReport {
    pub fn all_exact(&self) -> bool {
        self.reconstructed_ok
            && self.forward_exact.iter().all(|(_, b)| *b)
            && self.backward_exact.iter().all(|(_, b)| *b)
    }
    /// Hot-path HBP rows (json=0).
    pub fn to_hbp(&self) -> Vec<String> {
        let mut rows = vec![format!(
            "SLICETIME|reconstructed={}|forward_exact={}/{}|backward_exact={}/{}|all_exact={}|fire=0|json=0",
            self.reconstructed_ok,
            self.forward_exact.iter().filter(|(_, b)| *b).count(), self.forward_exact.len(),
            self.backward_exact.iter().filter(|(_, b)| *b).count(), self.backward_exact.len(),
            self.all_exact()
        )];
        rows.push("SLICETIME-BOUNDARY|zero_loss=determined-component-only|new_entropy=detected-not-predicted|schrodinger=free-because-simulated-not-quantum|json=0".to_string());
        rows
    }
}

/// The deterministic regime (steps 1-6): reconstruct t0 from shadows (multi-cylinder / Q-PRISM),
/// extrapolate the KNOWN rule forward and backward, and compare each predicted slice byte-for-byte
/// against generated truth. Exact in the determined regime; `Held` if shadows are insufficient.
pub fn run_deterministic(truth0: &SimState, horizon: i64, subset: &[usize]) -> Result<SliceTimeReport, Held> {
    let slice0 = truth0.serialize();
    let shadows = project_shadows(&slice0);
    let recon = reconstruct(&shadows, subset, slice0.len())?;
    let reconstructed_ok = recon == slice0;
    let state0 = SimState::deserialize(&recon).ok_or(Held::BadState)?;
    let (mut forward_exact, mut backward_exact) = (Vec::new(), Vec::new());
    for k in 1..=horizon {
        forward_exact.push((k, state0.at(k).serialize() == truth0.at(k).serialize()));
        backward_exact.push((-k, state0.at(-k).serialize() == truth0.at(-k).serialize()));
    }
    Ok(SliceTimeReport { reconstructed_ok, forward_exact, backward_exact })
}

/// The HONEST BOUNDARY (step 6): given a truth slice that GAINED entropy the rule cannot fix (an
/// injected particle, an external input), the deterministic prediction from `state0` will NOT match.
/// Returns Err(PredictionDivergedNewEntropy) - the break is DETECTED, never silently "predicted".
pub fn verify_or_detect_entropy(state0: &SimState, observed_future: &SimState, dt: i64) -> Result<(), Held> {
    if state0.at(dt).serialize() == observed_future.serialize() {
        Ok(()) // determined: 0-loss
    } else {
        Err(Held::PredictionDivergedNewEntropy) // new entropy: caught, not manufactured
    }
}

// ================================================================ WATCHER GATE (yin-yang)
// The classical consistency layer. Black-absorb (project shadows) -> white-emit (reconstruct) must
// round-trip to the original; every EXTRA cylinder beyond the sufficient subset is CONSISTENCY-
// CHECKED (over-determination = the hallucination-detection budget); anything that disagrees is
// Held, never emitted. Returns a VERIFIED CLASSICAL clone (a representation copy; no-cloning
// respected - NOT a physical quantum clone) or Held. Round-trip + consistency + Shannon-ledger =
// MEASURED; GNN-forward / reverse-GNN = DESIGN scaffold. Universe-analogue (boundary encodes bulk,
// forward/reverse gate state, extra observation reduces ambiguity, undetermined stays Held) = DESIGN
// - no claim this IS the universe engine.

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Verdict { VerifiedClone(Vec<u8>), Held(HeldReason) }
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HeldReason { InsufficientShadows, WatcherDisagreement }

/// The watcher gate: reconstruct from `sufficient`; then every `cross` cylinder must agree with the
/// reconstruction (catches a corrupted / hallucinated shadow). All agree -> verified classical clone.
pub fn watcher_gate(shadows: &[Vec<u64>], sufficient: &[usize], cross: &[usize], orig_len: usize) -> Verdict {
    let recon = match reconstruct(shadows, sufficient, orig_len) {
        Ok(r) => r,
        Err(_) => return Verdict::Held(HeldReason::InsufficientShadows),
    };
    for &ci in cross {
        let recomputed: Vec<u64> = blocks(&recon).iter().map(|&b| (b % CYLINDERS[ci] as u128) as u64).collect();
        if recomputed != shadows[ci] { return Verdict::Held(HeldReason::WatcherDisagreement); }
    }
    Verdict::VerifiedClone(recon)
}

/// FNV-1a 64-bit fold -> an 8-byte addressing digest (NOT a crypto hash).
fn digest8(slice: &[u8]) -> [u8; 8] {
    let mut h = 0xcbf29ce484222325u64;
    for &b in slice { h = (h ^ b as u64).wrapping_mul(0x100000001b3); }
    h.to_be_bytes()
}
fn hex8(b: &[u8; 8]) -> String { b.iter().map(|x| format!("{:02x}", x)).collect() }

/// The OMNIBIT PIXEL - one pixel carrying the representation STACK as SELECTORS + CHECKS, not raw
/// payload: position/tick/frequency, an 8-byte fold digest, the residual selector bits, the signed
/// over-determination margin (>=0 = redundancy; never literal negative storage), and the watcher
/// bound (# of consistency checks passed). It does NOT reconstruct the slice alone - it selects and
/// checks against the shared atlas. Not a magical all-information bit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OmnibitPixel {
    pub x: u16, pub y: u16, pub tick: u32, pub freq: u16,
    pub digest8: [u8; 8],
    pub residual_bits: u8,
    pub capacity_margin_bits: i32,
    pub watcher_bound: u8,
}
impl OmnibitPixel {
    pub fn of(slice: &[u8], x: u16, y: u16, tick: u32, freq: u16, n_cylinders: usize, watcher_bound: u8) -> Self {
        let roof = (n_cylinders.min(CYLINDERS.len()) as f64) * 25.0; // ~log2(2^25) per cylinder
        let block_bits = (8 * BLOCK) as f64;
        OmnibitPixel {
            x, y, tick, freq, digest8: digest8(slice),
            residual_bits: (block_bits - roof).max(0.0).ceil() as u8,
            capacity_margin_bits: (roof - block_bits) as i32,
            watcher_bound,
        }
    }
    pub fn to_hbp(&self) -> String {
        format!("OMNIBITPIXEL|x={}|y={}|tick={}|freq={}|digest8={}|residual_bits={}|capacity_margin_bits={}|watcher_bound={}|payload=selector-not-raw|json=0",
            self.x, self.y, self.tick, self.freq, hex8(&self.digest8), self.residual_bits, self.capacity_margin_bits, self.watcher_bound)
    }
}

// ================================================================ tests
#[cfg(test)]
mod tests {
    use super::*;

    fn sample_state() -> SimState {
        SimState {
            w: 64, h: 48,
            particles: vec![
                Particle { x: 3, y: 5, vx: 2, vy: -1, tag: 1 },
                Particle { x: 60, y: 1, vx: -3, vy: 4, tag: 2 },
                Particle { x: 31, y: 24, vx: 1, vy: 1, tag: 3 },
                Particle { x: 0, y: 47, vx: 5, vy: -2, tag: 4 },
            ],
        }
    }

    #[test]
    fn serialize_roundtrips() {
        let s = sample_state();
        assert_eq!(SimState::deserialize(&s.serialize()).unwrap(), {
            let mut c = s.clone(); c.particles.sort_by_key(|p| p.tag); c
        });
    }

    #[test]
    fn rule_is_reversible() {
        let s = sample_state();
        for k in 1..=20 {
            assert_eq!(s.at(k).at(-k).serialize(), s.serialize(), "at(+{k}).at(-{k}) must be identity");
        }
    }

    #[test]
    fn deterministic_forward_and_backward_are_byte_identical() {
        let s = sample_state();
        // enough shadows (all 4 cylinders) -> full state reconstructed -> every slice exact both ways
        let report = run_deterministic(&s, 25, &[0, 1, 2, 3]).unwrap();
        assert!(report.reconstructed_ok);
        assert!(report.all_exact(), "determined regime: past AND future slices are 0-loss");
        assert_eq!(report.forward_exact.len(), 25);
        assert_eq!(report.backward_exact.len(), 25);
    }

    #[test]
    fn two_shadows_already_suffice_one_holds() {
        let s = sample_state();
        assert!(run_deterministic(&s, 5, &[0, 1]).unwrap().all_exact()); // 2 cylinders cover a 48-bit block
        assert!(matches!(run_deterministic(&s, 5, &[0]), Err(Held::InsufficientShadows))); // 1 -> Shannon Held
    }

    #[test]
    fn new_entropy_breaks_0loss_and_is_detected_not_predicted() {
        let s = sample_state();
        // the DETERMINED future is 0-loss:
        assert!(verify_or_detect_entropy(&s, &s.at(7), 7).is_ok());
        // but a future that GAINED entropy (an injected particle) is DETECTED, never manufactured:
        let mut injected = s.at(7);
        injected.particles.push(Particle { x: 9, y: 9, vx: 0, vy: 0, tag: 9999 });
        assert_eq!(verify_or_detect_entropy(&s, &injected, 7), Err(Held::PredictionDivergedNewEntropy));
    }

    #[test]
    fn hbp_report_is_hotpath_json0() {
        let s = sample_state();
        let rows = run_deterministic(&s, 10, &[0, 1, 2, 3]).unwrap().to_hbp();
        assert!(rows.iter().all(|r| r.ends_with("json=0")) && !rows.join("").contains('{'));
        assert!(rows[1].contains("new_entropy=detected-not-predicted"));
    }

    #[test]
    fn watcher_gate_verifies_clean_and_catches_tamper() {
        let slice = b"the yin-yang watcher: black-absorb, white-emit, round-trip or Held";
        let shadows = project_shadows(slice);
        // clean: reconstruct from [0,1], cross-check [2,3] -> verified classical clone
        assert_eq!(watcher_gate(&shadows, &[0, 1], &[2, 3], slice.len()), Verdict::VerifiedClone(slice.to_vec()));
        // tamper a sufficient shadow -> the cross cylinders catch the hallucination -> Held
        let mut tampered = shadows.clone();
        tampered[0][0] = tampered[0][0].wrapping_add(1);
        assert_eq!(watcher_gate(&tampered, &[0, 1], &[2, 3], slice.len()), Verdict::Held(HeldReason::WatcherDisagreement));
        // too few -> Held (Shannon)
        assert_eq!(watcher_gate(&shadows, &[0], &[], slice.len()), Verdict::Held(HeldReason::InsufficientShadows));
    }

    #[test]
    fn omnibit_pixel_is_selector_not_payload() {
        let slice = b"omnibit pixel = selectors + checks, not raw payload";
        let px2 = OmnibitPixel::of(slice, 3, 5, 7, 440, 2, 2);
        assert_eq!(px2.residual_bits, 0); // 2 cylinders (~50 bits) cover the 48-bit block
        assert!(px2.capacity_margin_bits >= 0); // over-determination margin, never negative storage
        let px1 = OmnibitPixel::of(slice, 3, 5, 7, 440, 1, 0);
        assert!(px1.residual_bits > 0); // 1 cylinder under-covers -> residual selector remains
        assert!(px2.to_hbp().contains("payload=selector-not-raw") && px2.to_hbp().ends_with("json=0"));
    }
}
