//! Codon ORF statistics on fixed lattice (expression / spin — pure fixed math).
//! PRIMARY map itself remains exact integer trits in codon.zig.

const fixed = @import("fixed.zig");
const codon = @import("codon.zig");
const seeds_f = @import("seeds_fixed.zig");
const Fixed = fixed.Fixed;

/// expression = φ^spin · e^{|q|/(π·n)} · (1 + γ·aromatic)
pub fn geneExpression(residues: []const codon.Residue) Fixed {
    if (residues.len == 0) return fixed.fromInt(1);
    var trit_sum: i64 = 0;
    var n_trits: i64 = 0;
    var q: i32 = 0;
    var arom_n: i64 = 0;
    for (residues) |r| {
        trit_sum += r.trip[0] + r.trip[1] + r.trip[2];
        n_trits += 3;
        q += codon.aaCharge(r.aa);
        if (codon.aaIsAromatic(r.aa)) arom_n += 1;
    }
    const spin = fixed.div(fixed.fromInt(trit_sum), fixed.fromInt(n_trits));
    const n = fixed.fromInt(@intCast(residues.len));
    const arom = fixed.div(fixed.fromInt(arom_n), n);
    // phi^spin = exp(spin * log(phi)); spin~0 → 1
    const phi_pow = if (fixed.abs(spin) < fixed.fromDecimalStr("0.000001"))
        fixed.fromInt(1)
    else
        fixed.exp(fixed.mul(spin, fixed.log(seeds_f.phi)));
    const aq = if (q < 0) -q else q;
    const exp_q = if (aq == 0)
        fixed.fromInt(1)
    else
        fixed.exp(fixed.div(fixed.fromInt(aq), fixed.mul(seeds_f.pi, n)));
    const arom_term = fixed.add(fixed.fromInt(1), fixed.mul(seeds_f.gamma, arom));
    var raw = fixed.mul(fixed.mul(phi_pow, exp_q), arom_term);
    raw = fixed.clamp(raw, fixed.fromDecimalStr("0.05"), fixed.fromInt(3));
    return raw;
}

pub fn meanSpin(residues: []const codon.Residue) Fixed {
    if (residues.len == 0) return 0;
    var s: i64 = 0;
    var n: i64 = 0;
    for (residues) |r| {
        s += r.trip[0] + r.trip[1] + r.trip[2];
        n += 3;
    }
    return fixed.div(fixed.fromInt(s), fixed.fromInt(n));
}

pub fn selfTest() bool {
    const scn = "ATGAAATTTCGTTATTGG";
    var res: [codon.MAX_RESIDUES]codon.Residue = undefined;
    const nr = codon.decodeOrf(scn, &res);
    if (nr != 6) return false;
    const spin = meanSpin(res[0..nr]);
    if (fixed.abs(spin) > 2) return false; // ~0
    const expr = geneExpression(res[0..nr]);
    // Python ~1.43285
    const lo = fixed.fromDecimalStr("1.40");
    const hi = fixed.fromDecimalStr("1.47");
    if (fixed.lt(expr, lo) or fixed.gt(expr, hi)) return false;
    return true;
}
