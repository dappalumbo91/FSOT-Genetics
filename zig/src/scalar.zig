//! FSOT scalar S = K*(T1+T2+T3) — matches fsot_nuron/scalar.py (f64 path).

const std = @import("std");
const seeds = @import("seeds.zig");

pub fn clamp(x: f64, lo: f64, hi: f64) f64 {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

/// Full fluid scalar. observed=true applies quirk_mod (neuroscience).
pub fn computeScalar(
    N: f64,
    P: f64,
    D_eff: f64,
    recent_hits: f64,
    delta_psi: f64,
    delta_theta: f64,
    rho: f64,
    scale: f64,
    amplitude: f64,
    trend_bias: f64,
    observed: bool,
) f64 {
    const s = seeds;
    const D = if (D_eff < 1e-6) 1e-6 else D_eff;
    const Nn = if (N < 1e-6) 1e-6 else N;

    const growth = @exp(s.alpha * (1.0 - recent_hits / Nn) * s.gamma / s.phi);
    var t1 = (N * P / @sqrt(D)) *
        @cos((s.psi_con + delta_psi) / s.eta_eff) *
        @exp(-s.alpha * recent_hits / Nn + rho + s.b_in * delta_psi) *
        (1.0 + growth * s.c_eff);
    t1 *= (1.0 + s.p_new * @log(D / 25.0));
    if (observed) {
        t1 *= @exp(s.c_factor * s.p_var) * @cos(delta_psi + s.p_var);
    }

    const t2 = scale * amplitude + trend_bias;

    const valve = s.beta *
        @cos(delta_psi) *
        (N * P / @sqrt(D)) *
        (1.0 + s.chaos * (D - 25.0) / 25.0) *
        (1.0 + s.poof * @cos(s.theta_s + s.pi) + s.suction * @sin(s.theta_s));
    const acoustic = 1.0 +
        (s.a_bleed * std.math.pow(f64, @sin(delta_theta), 2) ) / s.phi +
        (s.a_in * std.math.pow(f64, @cos(delta_theta), 2) ) / s.phi;
    const phase = 1.0 + s.b_in * s.p_var;
    const t3 = valve * acoustic * phase;

    return clamp(s.k * (t1 + t2 + t3), -3.0, 3.0);
}

/// Convenience: neuroscience fold defaults, free delta_psi / hits / rho.
pub fn computeNeuro(delta_psi: f64, recent_hits: f64, rho: f64) f64 {
    return computeScalar(
        seeds.neuro_n_channels,
        seeds.neuro_p,
        seeds.neuro_d_eff,
        recent_hits,
        delta_psi,
        1.0,
        rho,
        1.0,
        1.0,
        0.0,
        true,
    );
}
