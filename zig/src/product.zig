//! Genetics product residual law — bare-metal authority twin of
//! scripts/residual_physics_refine.py + msa_template_fuse residual weights.
//!
//! Archive:
//!   residual_r = 1 + |S_domain| * P_NEW
//!   force_channel *= residual_r
//!
//! Named pin domains only (D1D38A table):
//!   bond   → Physical_Chemistry  D=8  δψ=0.5
//!   clash  → Chemistry           D=8  δψ=0.6
//!   anchor → Biochemistry        D=13 hits=1 δψ=0.35
//!
//! Zero free parameters. Measured template Cα is host/kernel data — not invented.

const scalar = @import("scalar.zig");
const seeds = @import("seeds.zig");
const codon = @import("codon.zig");
const trit = @import("trit.zig");

pub const ResidualChannels = struct {
    r_bond: f64,
    r_clash: f64,
    r_anchor: f64,
    s_physchem: f64,
    s_chem: f64,
    s_biochem: f64,
};

/// residual = 1 + |S| · P_NEW  (pin P_NEW from seeds)
pub fn residualScale(s: f64) f64 {
    const a = if (s < 0) -s else s;
    return 1.0 + a * seeds.p_new;
}

/// Domain scalar: N=P=1, rho=scale=amp=1 — matches vendor domain_scalar.
fn domainS(d_eff: f64, hits: f64, delta_psi: f64, delta_theta: f64, observed: bool) f64 {
    return scalar.computeScalar(
        1.0,
        1.0,
        d_eff,
        hits,
        delta_psi,
        delta_theta,
        1.0,
        1.0,
        1.0,
        0.0,
        observed,
    );
}

pub fn productResiduals() ResidualChannels {
    // Physical_Chemistry: D=8, hits=0, δψ=0.5, δθ=1, observed
    const s_pc = domainS(8.0, 0.0, 0.5, 1.0, true);
    // Chemistry: D=8, hits=0, δψ=0.6, δθ=1, observed
    const s_ch = domainS(8.0, 0.0, 0.6, 1.0, true);
    // Biochemistry: D=13, hits=1, δψ=0.35, δθ=1, observed
    const s_bc = domainS(13.0, 1.0, 0.35, 1.0, true);
    return .{
        .s_physchem = s_pc,
        .s_chem = s_ch,
        .s_biochem = s_bc,
        .r_bond = residualScale(s_pc),
        .r_clash = residualScale(s_ch),
        .r_anchor = residualScale(s_bc),
    };
}

/// Seed-closed multi-template ensemble knobs (matches Python MULTI_TOP_K / POWER).
pub const multi_top_k: usize = blk: {
    // round(φ³) ≈ 4
    const v = seeds.phi * seeds.phi * seeds.phi;
    break :blk @intFromFloat(@round(v));
};
pub const multi_power: f64 = seeds.phi * seeds.phi * seeds.phi * seeds.phi * seeds.phi * seeds.phi;

/// Mini product cell: DNA ORF → AA string + residual channels finite & in band.
pub fn productCellSelfTest() struct { ok: bool, aa_len: u32 } {
    // Short ORF fragment (insulin A-chain start region coding sample) — codon gate only
    const dna = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCC";
    var aa_count: u32 = 0;
    var i: usize = 0;
    while (i + 2 < dna.len) : (i += 3) {
        const aa = codon.dnaToAa(dna[i], dna[i + 1], dna[i + 2]);
        if (aa == '*') break;
        aa_count += 1;
    }
    const r = productResiduals();
    const finite = r.r_bond == r.r_bond and r.r_clash == r.r_clash and r.r_anchor == r.r_anchor;
    // Residuals must be slightly above 1 (law) and below 2 (sanity on pin S)
    const band = r.r_bond > 1.0 and r.r_bond < 2.0 and r.r_clash > 1.0 and r.r_clash < 2.0 and r.r_anchor > 1.0 and r.r_anchor < 2.0;
    // Cross-check known host residual magnitudes (~1.09–1.13 from Python pin)
    const near_py = r.r_bond > 1.05 and r.r_bond < 1.20 and r.r_clash > 1.05 and r.r_clash < 1.25;
    const atg = codon.primaryTrip('A', 'T', 'G');
    const codon_ok = atg[0] == 1 and atg[1] == -1 and atg[2] == 1;
    const trit_ok = trit.pair(1, -1) == -1 and trit.consensus(1, 1) == 1;
    return .{
        .ok = finite and band and near_py and codon_ok and trit_ok and aa_count >= 10,
        .aa_len = aa_count,
    };
}

/// Host residual self-test (no I/O).
pub fn residualSelfTest() bool {
    const r = productResiduals();
    if (r.r_bond != r.r_bond) return false;
    // |S|·P_NEW small → residual slightly > 1
    if (!(r.r_bond > 1.05 and r.r_bond < 1.20)) return false;
    if (!(r.r_clash > 1.05 and r.r_clash < 1.25)) return false;
    if (!(r.r_anchor > 1.05 and r.r_anchor < 1.20)) return false;
    // multi-template seeds
    if (multi_top_k < 3 or multi_top_k > 6) return false;
    if (!(multi_power > 10.0 and multi_power < 30.0)) return false;
    return true;
}

/// Ideal CA–CA spacing used by structure engine / residual physics (Å).
pub const CA_CA: f64 = 3.8;

/// One residual-weighted bond step on a 3-point chain (parity with Python).
/// X = [(0,0,0), (4.0,0,0), (8.0,0,0)] stretched vs CA_CA; one gradient step.
pub const PhysSample = struct {
    bond_len: f64,
    end_x: f64,
};

pub fn residualPhysicsParitySample() PhysSample {
    const r = productResiduals();
    // Match residual_physics_refine: lr=0.08, W_ANCHOR0=0.05 * r_anchor
    const lr: f64 = 0.08;
    const w_anchor = 0.05 * r.r_anchor;
    const x0 = [_][3]f64{
        .{ 0.0, 0.0, 0.0 },
        .{ 4.0, 0.0, 0.0 },
        .{ 8.0, 0.0, 0.0 },
    };
    var x = x0;
    // one iteration
    var g = [_][3]f64{ .{ 0, 0, 0 }, .{ 0, 0, 0 }, .{ 0, 0, 0 } };
    // anchor
    var i: usize = 0;
    while (i < 3) : (i += 1) {
        g[i][0] = w_anchor * (x[i][0] - x0[i][0]);
        g[i][1] = w_anchor * (x[i][1] - x0[i][1]);
        g[i][2] = w_anchor * (x[i][2] - x0[i][2]);
    }
    // bonds
    i = 0;
    while (i < 2) : (i += 1) {
        const dx = x[i + 1][0] - x[i][0];
        const dy = x[i + 1][1] - x[i][1];
        const dz = x[i + 1][2] - x[i][2];
        const L = @sqrt(dx * dx + dy * dy + dz * dz) + 1e-9;
        const scale = ((L - CA_CA) / L) * r.r_bond;
        const fx = scale * dx;
        const fy = scale * dy;
        const fz = scale * dz;
        g[i][0] -= fx;
        g[i][1] -= fy;
        g[i][2] -= fz;
        g[i + 1][0] += fx;
        g[i + 1][1] += fy;
        g[i + 1][2] += fz;
    }
    i = 0;
    while (i < 3) : (i += 1) {
        x[i][0] -= lr * g[i][0];
        x[i][1] -= lr * g[i][1];
        x[i][2] -= lr * g[i][2];
    }
    // center
    var mx: f64 = 0;
    i = 0;
    while (i < 3) : (i += 1) mx += x[i][0];
    mx /= 3.0;
    i = 0;
    while (i < 3) : (i += 1) x[i][0] -= mx;
    const bdx = x[1][0] - x[0][0];
    const bdy = x[1][1] - x[0][1];
    const bdz = x[1][2] - x[0][2];
    const bl = @sqrt(bdx * bdx + bdy * bdy + bdz * bdz);
    return .{ .bond_len = bl, .end_x = x[2][0] };
}
