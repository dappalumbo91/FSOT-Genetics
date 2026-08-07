//! FSOT 64-codon trinary foundation — bare-metal authority twin of
//! data/64_codon_trinary_map.txt + chemical_codon.py.
//!
//! Law (PRIMARY): A,G → +1 ; C,T → −1  per base → codon ∈ {-1,+1}³
//! Law (SECONDARY): A=+1, T=−1, G/C=0
//! Then: DNA codon → AA (IUPAC) → process / charge for gene programs.
//!
//! This is not optional decoration. Neuron genetics and W_ij rest on it.

const trit = @import("trit.zig");
const seeds = @import("seeds.zig");

pub const Trit = trit.Trit;
pub const CodonTrip = [3]Trit;

/// Encode codon string "ATG" as packed u15 key: 5 bits per base (A=1,C=2,G=3,T=4).
pub fn baseCode(b: u8) u8 {
    return switch (b) {
        'A', 'a' => 1,
        'C', 'c' => 2,
        'G', 'g' => 3,
        'T', 't', 'U', 'u' => 4,
        else => 0,
    };
}

pub fn codonKey(c0: u8, c1: u8, c2: u8) u16 {
    return (@as(u16, baseCode(c0)) << 8) | (@as(u16, baseCode(c1)) << 4) | @as(u16, baseCode(c2));
}

/// PRIMARY map — exact authority of 64_codon_trinary_map.txt
pub fn primaryTrip(c0: u8, c1: u8, c2: u8) CodonTrip {
    return trit.codonPrimary(c0, c1, c2);
}

/// SECONDARY map: A=+1, T=-1, G/C=0
pub fn secondaryBase(b: u8) Trit {
    return switch (b) {
        'A', 'a' => 1,
        'T', 't', 'U', 'u' => -1,
        'G', 'g', 'C', 'c' => 0,
        else => 0,
    };
}

pub fn secondaryTrip(c0: u8, c1: u8, c2: u8) CodonTrip {
    return .{ secondaryBase(c0), secondaryBase(c1), secondaryBase(c2) };
}

/// IUPAC DNA → AA (standard genetic code). '*' = stop.
pub fn dnaToAa(c0: u8, c1: u8, c2: u8) u8 {
    // Normalize to upper
    const a = upperBase(c0);
    const b = upperBase(c1);
    const c = upperBase(c2);
    // Compact nested switch on first base
    return switch (a) {
        'T' => switch (b) {
            'T' => switch (c) {
                'T', 'C' => 'F',
                'A', 'G' => 'L',
                else => '?',
            },
            'C' => 'S', // TCN all S
            'A' => switch (c) {
                'T', 'C' => 'Y',
                'A', 'G' => '*',
                else => '?',
            },
            'G' => switch (c) {
                'T', 'C' => 'C',
                'A' => '*',
                'G' => 'W',
                else => '?',
            },
            else => '?',
        },
        'C' => switch (b) {
            'T' => 'L',
            'C' => 'P',
            'A' => switch (c) {
                'T', 'C' => 'H',
                'A', 'G' => 'Q',
                else => '?',
            },
            'G' => 'R',
            else => '?',
        },
        'A' => switch (b) {
            'T' => switch (c) {
                'T', 'C', 'A' => 'I',
                'G' => 'M',
                else => '?',
            },
            'C' => 'T',
            'A' => switch (c) {
                'T', 'C' => 'N',
                'A', 'G' => 'K',
                else => '?',
            },
            'G' => switch (c) {
                'T', 'C' => 'S',
                'A', 'G' => 'R',
                else => '?',
            },
            else => '?',
        },
        'G' => switch (b) {
            'T' => 'V',
            'C' => 'A',
            'A' => switch (c) {
                'T', 'C' => 'D',
                'A', 'G' => 'E',
                else => '?',
            },
            'G' => 'G',
            else => '?',
        },
        else => '?',
    };
}

fn upperBase(b: u8) u8 {
    return switch (b) {
        'a' => 'A',
        'c' => 'C',
        'g' => 'G',
        't', 'u' => 'T',
        'U' => 'T',
        else => b,
    };
}

pub fn aaCharge(aa: u8) i8 {
    return switch (aa) {
        'R', 'H', 'K' => 1,
        'D', 'E' => -1,
        else => 0,
    };
}

pub fn aaIsAromatic(aa: u8) bool {
    return aa == 'F' or aa == 'Y' or aa == 'W';
}

pub fn aaIsHydrophobic(aa: u8) bool {
    return aa == 'A' or aa == 'I' or aa == 'L' or aa == 'M' or aa == 'F' or aa == 'V' or aa == 'W';
}

pub const Residue = struct {
    c0: u8,
    c1: u8,
    c2: u8,
    trip: CodonTrip,
    aa: u8,
};

pub const MAX_RESIDUES: usize = 16; // short ORFs fit

/// Decode DNA ORF (length multiple of 3) → residues with PRIMARY trinary.
pub fn decodeOrf(dna: []const u8, out: *[MAX_RESIDUES]Residue) usize {
    var n: usize = 0;
    var i: usize = 0;
    while (i + 2 < dna.len and n < MAX_RESIDUES) : (i += 3) {
        const a = upperBase(dna[i]);
        const b = upperBase(dna[i + 1]);
        const c = upperBase(dna[i + 2]);
        if (baseCode(a) == 0 or baseCode(b) == 0 or baseCode(c) == 0) continue;
        out[n] = .{
            .c0 = a,
            .c1 = b,
            .c2 = c,
            .trip = primaryTrip(a, b, c),
            .aa = dnaToAa(a, b, c),
        };
        n += 1;
    }
    return n;
}

/// expression = φ^spin · e^{|q|/(π·n)} · (1 + γ·aromatic)  — seeds only.
pub fn geneExpression(residues: []const Residue) f64 {
    if (residues.len == 0) return 1.0;
    var trit_sum: f64 = 0;
    var n_trits: f64 = 0;
    var q: i32 = 0;
    var arom_n: f64 = 0;
    for (residues) |r| {
        trit_sum += @as(f64, @floatFromInt(r.trip[0]));
        trit_sum += @as(f64, @floatFromInt(r.trip[1]));
        trit_sum += @as(f64, @floatFromInt(r.trip[2]));
        n_trits += 3;
        q += aaCharge(r.aa);
        if (aaIsAromatic(r.aa)) arom_n += 1;
    }
    const spin = trit_sum / n_trits;
    const n = @as(f64, @floatFromInt(residues.len));
    const arom = arom_n / n;
    const raw = stdPow(seeds.phi, spin) *
        @exp(@as(f64, @floatFromInt(if (q < 0) -q else q)) / (seeds.pi * n)) *
        (1.0 + seeds.gamma * arom);
    if (raw < 0.05) return 0.05;
    if (raw > 3.0) return 3.0;
    return raw;
}

fn stdPow(base: f64, expv: f64) f64 {
    if (base <= 0) return 0;
    return @exp(@log(base) * expv);
}

pub fn meanSpin(residues: []const Residue) f64 {
    if (residues.len == 0) return 0;
    var s: f64 = 0;
    var n: f64 = 0;
    for (residues) |r| {
        s += @as(f64, @floatFromInt(r.trip[0] + r.trip[1] + r.trip[2]));
        n += 3;
    }
    return s / n;
}

pub fn chargeBalance(residues: []const Residue) i32 {
    var q: i32 = 0;
    for (residues) |r| q += aaCharge(r.aa);
    return q;
}

pub fn aromaticFraction(residues: []const Residue) f64 {
    if (residues.len == 0) return 0;
    var a: f64 = 0;
    for (residues) |r| {
        if (aaIsAromatic(r.aa)) a += 1;
    }
    return a / @as(f64, @floatFromInt(residues.len));
}

/// F01-style protein phase (charge, polarity, volume) ∈ {-1,0,+1}³
/// Exact twin of genetic_genotype.aa_trinary_phase.
pub fn aaTrinaryPhase(aa: u8) CodonTrip {
    const ch = aaCharge(aa);
    const c: Trit = if (ch > 0) 1 else if (ch < 0) -1 else 0;
    const polar = aa == 'S' or aa == 'T' or aa == 'N' or aa == 'Q' or aa == 'Y' or aa == 'C';
    const hydro = aaIsHydrophobic(aa);
    const p: Trit = if (polar) 1 else if (hydro) -1 else 0;
    const large = aa == 'I' or aa == 'L' or aa == 'M' or aa == 'F' or aa == 'W' or aa == 'Y' or aa == 'R' or aa == 'K' or aa == 'Q' or aa == 'E';
    const small = aa == 'A' or aa == 'G' or aa == 'S' or aa == 'T' or aa == 'C';
    const v: Trit = if (large) 1 else if (small) -1 else 0;
    return .{ c, p, v };
}

/// Verify all 64 codons: PRIMARY matches AG/CT law; ATG = [1,-1,1]; TTT = [-1,-1,-1].
pub fn selfTest() bool {
    const bases = [_]u8{ 'A', 'C', 'G', 'T' };
    var n: u32 = 0;
    for (bases) |b0| {
        for (bases) |b1| {
            for (bases) |b2| {
                const t = primaryTrip(b0, b1, b2);
                const e0: Trit = if (b0 == 'A' or b0 == 'G') 1 else -1;
                const e1: Trit = if (b1 == 'A' or b1 == 'G') 1 else -1;
                const e2: Trit = if (b2 == 'A' or b2 == 'G') 1 else -1;
                if (t[0] != e0 or t[1] != e1 or t[2] != e2) return false;
                n += 1;
            }
        }
    }
    if (n != 64) return false;

    const atg = primaryTrip('A', 'T', 'G');
    if (atg[0] != 1 or atg[1] != -1 or atg[2] != 1) return false;
    const ttt = primaryTrip('T', 'T', 'T');
    if (ttt[0] != -1 or ttt[1] != -1 or ttt[2] != -1) return false;
    if (dnaToAa('A', 'T', 'G') != 'M') return false;
    if (dnaToAa('T', 'T', 'T') != 'F') return false;
    if (dnaToAa('T', 'A', 'A') != '*') return false;

    // SCN ORF ground path
    const scn = "ATGAAATTTCGTTATTGG";
    var res: [MAX_RESIDUES]Residue = undefined;
    const nr = decodeOrf(scn, &res);
    if (nr != 6) return false;
    if (res[0].aa != 'M') return false;
    const spin = meanSpin(res[0..nr]);
    if (@abs(spin) > 1e-9) return false; // SCN spin = 0
    const expr = geneExpression(res[0..nr]);
    // Python: ~1.43285
    if (expr < 1.42 or expr > 1.45) return false;

    // secondary ATG: A=+1,T=-1,G=0 → [1,-1,0]
    const sec = secondaryTrip('A', 'T', 'G');
    if (sec[0] != 1 or sec[1] != -1 or sec[2] != 0) return false;

    // aa phase ground truth (Python): M (0,-1,1) K (1,0,1) D (-1,0,0)
    const m = aaTrinaryPhase('M');
    if (m[0] != 0 or m[1] != -1 or m[2] != 1) return false;
    const k = aaTrinaryPhase('K');
    if (k[0] != 1 or k[1] != 0 or k[2] != 1) return false;
    const d = aaTrinaryPhase('D');
    if (d[0] != -1 or d[1] != 0 or d[2] != 0) return false;
    return true;
}
