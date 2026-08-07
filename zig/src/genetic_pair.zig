//! Genetic pair geometry — standalone extract for FSOT-Genetics.
//! Authority twin of fsot-neuron-zig genetic.zig pair laws (no network deps).
//! Seeds only. Zero free parameters. Pin lineage D1D38A.

const seeds = @import("seeds.zig");

/// trinaryPairInteraction(τ_i, τ_j) = τ_i·τ_j·e + (1−|τ_i·τ_j|)·π
pub fn trinaryPairInteraction(tau_i: f64, tau_j: f64) f64 {
    var ti = tau_i;
    var tj = tau_j;
    if (ti < -1) ti = -1;
    if (ti > 1) ti = 1;
    if (tj < -1) tj = -1;
    if (tj > 1) tj = 1;
    const prod = ti * tj;
    return prod * seeds.e + (1.0 - @abs(prod)) * seeds.pi;
}

/// geometricScaleDist = φ · dist^(−1/π)   (dist ≥ 1)
pub fn geometricScaleDist(dist: usize) f64 {
    const d: f64 = @floatFromInt(if (dist < 1) 1 else dist);
    return seeds.phi * @exp(@log(d) * (-1.0 / seeds.pi));
}

pub fn electrostaticTerm(q_i: f64, q_j: f64) f64 {
    return -q_i * q_j * seeds.e;
}

/// env(s) = s / (s + π·e)  — same F08 contact envelope as protein F15
pub fn envScale(dist: usize) f64 {
    const d: f64 = @floatFromInt(if (dist < 1) 1 else dist);
    return d / (d + seeds.pi * seeds.e);
}

/// Full pair weight (neuron genetic.zig / protein F07–F08 spine):
/// geom · (base + 0.15·elec) · (0.35 + 0.65·env)
pub fn fsotPairWeight(spin_i: f64, spin_j: f64, charge_i: f64, charge_j: f64, dist: usize) f64 {
    const base = trinaryPairInteraction(spin_i, spin_j);
    const geom = geometricScaleDist(dist);
    const elec = electrostaticTerm(charge_i, charge_j);
    const env = envScale(dist);
    return geom * (base + 0.15 * elec) * (0.35 + 0.65 * env);
}
