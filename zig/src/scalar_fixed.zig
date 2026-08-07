//! FSOT scalar S = K*(T1+T2+T3) in pure fixed-point (no IEEE float ops).
//! Twin structure of scalar.zig; constants from seeds_fixed.

const fixed = @import("fixed.zig");
const seeds = @import("seeds_fixed.zig");
const Fixed = fixed.Fixed;

pub fn clampS(x: Fixed) Fixed {
    return fixed.clamp(x, fixed.fromInt(-3), fixed.fromInt(3));
}

pub fn computeScalar(
    N: Fixed,
    P: Fixed,
    D_eff: Fixed,
    recent_hits: Fixed,
    delta_psi: Fixed,
    delta_theta: Fixed,
    rho: Fixed,
    scale: Fixed,
    amplitude: Fixed,
    trend_bias: Fixed,
    observed: bool,
) Fixed {
    const s = seeds;
    var D = D_eff;
    if (fixed.lt(D, fixed.fromDecimalStr("0.000001"))) D = fixed.fromDecimalStr("0.000001");
    var Nn = N;
    if (fixed.lt(Nn, fixed.fromDecimalStr("0.000001"))) Nn = fixed.fromDecimalStr("0.000001");

    // growth = exp(alpha * (1 - hits/N) * gamma / phi)
    const one = fixed.fromInt(1);
    const hits_over = fixed.div(recent_hits, Nn);
    const inner = fixed.mul(
        fixed.mul(s.alpha, fixed.sub(one, hits_over)),
        fixed.div(s.gamma, s.phi),
    );
    const growth = fixed.exp(inner);

    // t1 = (N*P/sqrt(D)) * cos((psi+dpsi)/eta) * exp(-alpha*hits/N + rho + b_in*dpsi) * (1+growth*c_eff)
    const np = fixed.mul(N, P);
    const sq = fixed.sqrt(D);
    const base = fixed.div(np, sq);
    const cos_arg = fixed.div(fixed.add(s.psi_con, delta_psi), s.eta_eff);
    const cosv = fixed.cos(cos_arg);
    const exp_arg = fixed.add(
        fixed.add(fixed.negate(fixed.mul(s.alpha, hits_over)), rho),
        fixed.mul(s.b_in, delta_psi),
    );
    const expv = fixed.exp(exp_arg);
    var t1 = fixed.mul(fixed.mul(base, cosv), expv);
    t1 = fixed.mul(t1, fixed.add(one, fixed.mul(growth, s.c_eff)));
    // * (1 + p_new * log(D/25))
    const logv = fixed.log(fixed.div(D, fixed.fromInt(25)));
    t1 = fixed.mul(t1, fixed.add(one, fixed.mul(s.p_new, logv)));
    if (observed) {
        const q = fixed.mul(fixed.exp(fixed.mul(s.c_factor, s.p_var)), fixed.cos(fixed.add(delta_psi, s.p_var)));
        t1 = fixed.mul(t1, q);
    }

    const t2 = fixed.add(fixed.mul(scale, amplitude), trend_bias);

    // t3 = beta * cos(dpsi) * (N*P/sqrt(D)) * (1+chaos*(D-25)/25) * (1+poof*cos(theta+pi)+suction*sin(theta))
    // beta ~ 2e-17 → 0 at SCALE 1e12; keep structure for parity of form
    const valve_core = fixed.mul(s.beta, fixed.mul(fixed.cos(delta_psi), base));
    const d25 = fixed.div(fixed.sub(D, fixed.fromInt(25)), fixed.fromInt(25));
    const chaos_term = fixed.add(one, fixed.mul(s.chaos, d25));
    const th = fixed.add(s.theta_s, s.pi);
    const phase_mod = fixed.add(
        one,
        fixed.add(fixed.mul(s.poof, fixed.cos(th)), fixed.mul(s.suction, fixed.sin(th))),
    );
    // acoustic
    const st = fixed.sin(delta_theta);
    const ct = fixed.cos(delta_theta);
    const acoustic = fixed.add(
        one,
        fixed.add(
            fixed.div(fixed.mul(s.a_bleed, fixed.mul(st, st)), s.phi),
            fixed.div(fixed.mul(s.a_in, fixed.mul(ct, ct)), s.phi),
        ),
    );
    const phase = fixed.add(one, fixed.mul(s.b_in, s.p_var));
    const t3 = fixed.mul(fixed.mul(fixed.mul(valve_core, chaos_term), phase_mod), fixed.mul(acoustic, phase));

    const raw = fixed.mul(s.k, fixed.add(fixed.add(t1, t2), t3));
    return clampS(raw);
}

pub fn computeNeuro(delta_psi: Fixed, recent_hits: Fixed, rho: Fixed) Fixed {
    return computeScalar(
        seeds.neuro_n_channels,
        seeds.neuro_p,
        seeds.neuro_d_eff,
        recent_hits,
        delta_psi,
        fixed.fromInt(1),
        rho,
        fixed.fromInt(1),
        fixed.fromInt(1),
        0,
        true,
    );
}

pub fn selfTest() bool {
    if (!fixed.selfTest()) return false;
    const s0 = computeNeuro(fixed.fromDecimalStr("0.1"), 0, fixed.fromInt(1));
    // f64 reference ~0.4323; fixed should be in same neighborhood
    const lo = fixed.fromDecimalStr("0.30");
    const hi = fixed.fromDecimalStr("0.60");
    if (fixed.lt(s0, lo) or fixed.gt(s0, hi)) return false;
    return true;
}
