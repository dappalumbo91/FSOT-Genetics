//! Neuron genotype: DNA ORF → 64-codon trinary → channel genes → phenotype.
//! Full port of fsot_nuron/genetic_genotype.py + cell_types ORF path.
//! Authority: codon.zig (PRIMARY map) + seeds — zero free-fit weights.

const codon = @import("codon.zig");
const seeds = @import("seeds.zig");
const bio_probe = @import("bio_probe.zig");
const cell_types = @import("cell_types.zig");

pub const GeneName = enum(u8) { scn = 0, kcn = 1, cacna = 2, leak = 3 };

/// Canonical short ORFs (DNA) — same strings as Python CHANNEL_GENE_ORFS.
pub const ORF_SCN = "ATGAAATTTCGTTATTGG";
pub const ORF_KCN = "ATGCTGGTTTCATCTTAG";
pub const ORF_CACNA = "ATGGATGAGTGTTATTGA";
pub const ORF_LEAK = "ATGGGTGCAAGCTCTTAA";

/// Cell-type ORF overrides (Python CELL_TYPES).
pub const ORF_PYR_SCN = "ATGAAACGGTTCTATTGG";
pub const ORF_PYR_CACNA = "ATGGATGAGTGCTACTGA";
pub const ORF_PYR_KCN = "ATGCTGGTTTCCAGTTAG";
pub const ORF_PYR_LEAK = "ATGGGTGCAAGCTCTTAA";

pub const ORF_PV_SCN = "ATGAAATTTAAGCGTTGG";
pub const ORF_PV_KCN = "ATGAAACGTGTTTCGTAG";
pub const ORF_PV_CACNA = "ATGGCTGCATCCTCTTAG";
pub const ORF_PV_LEAK = "ATGGGTGCTAACTCTTAA";

pub const ORF_SST_SCN = "ATGAAATTCCGCTATTGA";
pub const ORF_SST_KCN = "ATGCTGGTTACATCTTAA";
pub const ORF_SST_CACNA = "ATGGATGACTGCTATTGA";
pub const ORF_SST_LEAK = "ATGGCAGCAAGCTCTTAA";

pub const ORF_VIP_SCN = "ATGAAACAGTTCTATTAA";
pub const ORF_VIP_KCN = "ATGTTGGTTTCTTCTTAA";
pub const ORF_VIP_CACNA = "ATGGAGGAGTGTTCTTGA";
pub const ORF_VIP_LEAK = "ATGGGTTCAAACTCTTAA";

pub const GeneProgram = struct {
    name: GeneName,
    spin: f64 = 0,
    expression: f64 = 1,
    charge_balance: i32 = 0,
    aromatic_fraction: f64 = 0,
    mean_trinary: [3]f64 = .{ 0, 0, 0 },
    n_codons: u8 = 0,
};

pub const Phenotype = struct {
    n_channels: f64 = 4,
    p_props: f64 = 3,
    d_eff: f64 = 13,
    fire_threshold: f64 = 1.05,
    refractory_steps: f64 = 12,
    adapt_step: f64 = 0.7,
    adapt_gain: f64 = 0.02,
    adapt_decay: f64 = 0.988,
    resting_bias: f64 = 0.46,
    fi_stim: f64 = 0.5,
    vrest_mV: f64 = -70,
    composite_spin: f64 = 0,
    composite_charge: f64 = 0,
    scn_expression: f64 = 1,
    kcn_expression: f64 = 1,
    cacna_expression: f64 = 1,
    leak_expression: f64 = 1,
};

pub const NeuronGenotype = struct {
    unit_id: u32 = 0,
    cell_type: cell_types.CellType = .pyr,
    synapse_sign: i8 = 1,
    genes: [4]GeneProgram = undefined,
    phenotype: Phenotype = .{},
    composite_spin: f64 = 0,
    composite_charge: f64 = 0,
};

fn baseOrf(name: GeneName) []const u8 {
    return switch (name) {
        .scn => ORF_SCN,
        .kcn => ORF_KCN,
        .cacna => ORF_CACNA,
        .leak => ORF_LEAK,
    };
}

fn cellOrf(ct: cell_types.CellType, name: GeneName) []const u8 {
    return switch (ct) {
        .pyr => switch (name) {
            .scn => ORF_PYR_SCN,
            .kcn => ORF_PYR_KCN,
            .cacna => ORF_PYR_CACNA,
            .leak => ORF_PYR_LEAK,
        },
        .pv => switch (name) {
            .scn => ORF_PV_SCN,
            .kcn => ORF_PV_KCN,
            .cacna => ORF_PV_CACNA,
            .leak => ORF_PV_LEAK,
        },
        .sst => switch (name) {
            .scn => ORF_SST_SCN,
            .kcn => ORF_SST_KCN,
            .cacna => ORF_SST_CACNA,
            .leak => ORF_SST_LEAK,
        },
        .vip => switch (name) {
            .scn => ORF_VIP_SCN,
            .kcn => ORF_VIP_KCN,
            .cacna => ORF_VIP_CACNA,
            .leak => ORF_VIP_LEAK,
        },
    };
}

fn expressionBias(ct: cell_types.CellType, name: GeneName) f64 {
    // Python CELL_TYPES expression_bias
    return switch (ct) {
        .pyr => switch (name) {
            .scn => 1.15,
            .kcn => 0.95,
            .cacna => 1.10,
            .leak => 1.0,
        },
        .pv => switch (name) {
            .scn => 1.25,
            .kcn => 1.45,
            .cacna => 0.70,
            .leak => 0.95,
        },
        .sst => switch (name) {
            .scn => 0.95,
            .kcn => 1.05,
            .cacna => 1.35,
            .leak => 1.05,
        },
        .vip => switch (name) {
            .scn => 1.05,
            .kcn => 0.90,
            .cacna => 1.15,
            .leak => 1.10,
        },
    };
}

/// Deterministic diversity: flip purine/pyrimidine class at index (same primary trit class).
/// Trinary-preserving base flips (purine↔purine, pyrimidine↔pyrimidine).
/// Number of sites scales with unit_id so cohort diversity is genetic, not free noise.
pub fn mutateOrf(dna_in: []const u8, unit_id: u32, locus: u32, out: *[32]u8) []const u8 {
    const n = @min(dna_in.len, 32);
    @memcpy(out[0..n], dna_in[0..n]);
    if (n == 0) return out[0..0];
    const n_mut: u32 = 1 + (unit_id % 4); // 1..4 sites
    var m: u32 = 0;
    while (m < n_mut) : (m += 1) {
        const idx = (unit_id *% 3 +% locus *% 5 +% m *% 7) % @as(u32, @intCast(n));
        const b = out[idx];
        out[idx] = switch (b) {
            'A', 'a' => 'G',
            'G', 'g' => 'A',
            'C', 'c' => 'T',
            'T', 't', 'U', 'u' => 'C',
            else => b,
        };
    }
    return out[0..n];
}

pub fn buildGeneProgram(name: GeneName, dna: []const u8) GeneProgram {
    var res: [codon.MAX_RESIDUES]codon.Residue = undefined;
    const nr = codon.decodeOrf(dna, &res);
    var g: GeneProgram = .{ .name = name };
    if (nr == 0) return g;
    g.n_codons = @intCast(nr);
    g.spin = codon.meanSpin(res[0..nr]);
    g.expression = codon.geneExpression(res[0..nr]);
    g.charge_balance = codon.chargeBalance(res[0..nr]);
    g.aromatic_fraction = codon.aromaticFraction(res[0..nr]);
    var m0: f64 = 0;
    var m1: f64 = 0;
    var m2: f64 = 0;
    for (res[0..nr]) |r| {
        m0 += @floatFromInt(r.trip[0]);
        m1 += @floatFromInt(r.trip[1]);
        m2 += @floatFromInt(r.trip[2]);
    }
    const nf: f64 = @floatFromInt(nr);
    g.mean_trinary = .{ m0 / nf, m1 / nf, m2 / nf };
    return g;
}

pub fn phenotypeFromGenes(genes: *const [4]GeneProgram) Phenotype {
    const scn = genes[@intFromEnum(GeneName.scn)].expression;
    const kcn = genes[@intFromEnum(GeneName.kcn)].expression;
    const ca = genes[@intFromEnum(GeneName.cacna)].expression;
    const leak = genes[@intFromEnum(GeneName.leak)].expression;
    const mean_expr = (scn + kcn + ca + leak) / 4.0;

    var d_eff = seeds.neuro_d_eff + seeds.phi * (ca - 1.0) * 0.35 + seeds.gamma * (scn - kcn) * 0.2;
    if (d_eff < 8) d_eff = 8;
    if (d_eff > 20) d_eff = 20;

    var fire_threshold = 1.05 - 0.12 * (scn - 1.0) + 0.06 * (kcn - 1.0);
    if (fire_threshold < 0.85) fire_threshold = 0.85;
    if (fire_threshold > 1.25) fire_threshold = 1.25;

    var refractory_steps = 12.0 * (0.85 + 0.30 * kcn);
    if (refractory_steps < 4) refractory_steps = 4;
    if (refractory_steps > 40) refractory_steps = 40;

    var adapt_step = 0.7 * (0.6 + 0.8 * ca);
    if (adapt_step < 0) adapt_step = 0;
    if (adapt_step > 8) adapt_step = 8;

    const adapt_gain = 0.02 * (0.7 + 0.6 * ca);
    var adapt_decay = 0.988 - 0.004 * (ca - 1.0);
    if (adapt_decay < 0.96) adapt_decay = 0.96;
    if (adapt_decay > 0.995) adapt_decay = 0.995;

    var resting_bias = seeds.resting_s + 0.02 * (scn - leak) * seeds.eta_eff;
    if (resting_bias < 0.30) resting_bias = 0.30;
    if (resting_bias > 0.65) resting_bias = 0.65;

    var fi_stim = 0.50 * (0.85 + 0.25 * scn - 0.10 * kcn);
    if (fi_stim < 0.25) fi_stim = 0.25;
    if (fi_stim > 0.95) fi_stim = 0.95;

    var vrest = -70.0 - 3.0 * (leak - 1.0) + 1.5 * (scn - 1.0);
    if (vrest < -80) vrest = -80;
    if (vrest > -55) vrest = -55;

    var sum_spin_w: f64 = 0;
    var sum_expr: f64 = 0;
    var sum_q: f64 = 0;
    var i: usize = 0;
    while (i < 4) : (i += 1) {
        sum_spin_w += genes[i].spin * genes[i].expression;
        sum_expr += genes[i].expression;
        sum_q += @as(f64, @floatFromInt(genes[i].charge_balance)) * genes[i].expression;
    }
    const composite_spin = if (sum_expr > 1e-12) sum_spin_w / sum_expr else 0;
    const composite_charge = sum_q / 4.0;

    return .{
        .n_channels = seeds.neuro_n_channels * (0.75 + 0.25 * mean_expr),
        .p_props = seeds.neuro_p,
        .d_eff = d_eff,
        .fire_threshold = fire_threshold,
        .refractory_steps = refractory_steps,
        .adapt_step = adapt_step,
        .adapt_gain = adapt_gain,
        .adapt_decay = adapt_decay,
        .resting_bias = resting_bias,
        .fi_stim = fi_stim,
        .vrest_mV = vrest,
        .composite_spin = composite_spin,
        .composite_charge = composite_charge,
        .scn_expression = scn,
        .kcn_expression = kcn,
        .cacna_expression = ca,
        .leak_expression = leak,
    };
}

fn applyClassNudge(ct: cell_types.CellType, ph: *Phenotype) void {
    switch (ct) {
        .pv => {
            ph.refractory_steps = @max(3.0, ph.refractory_steps * 0.45);
            ph.adapt_step *= 0.35;
            ph.fire_threshold = @max(0.85, ph.fire_threshold - 0.04);
            ph.fi_stim = @min(0.95, ph.fi_stim * 1.15);
        },
        .sst => {
            ph.adapt_step = @min(10.0, ph.adapt_step * 1.4);
            ph.refractory_steps *= 1.05;
        },
        .vip => {
            ph.fi_stim *= 0.9;
            ph.d_eff = @min(20.0, ph.d_eff + 0.3 * seeds.phi);
        },
        .pyr => {
            ph.adapt_step *= 1.05;
        },
    }
}

/// Build full genotype for one unit of a cell type (codon path, optional diversity).
pub fn buildCellTypeGenotype(unit_id: u32, ct: cell_types.CellType, diversity: bool) NeuronGenotype {
    var gt: NeuronGenotype = .{
        .unit_id = unit_id,
        .cell_type = ct,
        .synapse_sign = cell_types.specOf(ct).sign,
    };
    const names = [_]GeneName{ .scn, .kcn, .cacna, .leak };
    var mutbuf: [32]u8 = undefined;
    var i: usize = 0;
    while (i < 4) : (i += 1) {
        const name = names[i];
        const base = cellOrf(ct, name);
        const dna = if (diversity)
            mutateOrf(base, unit_id, @intCast(i), &mutbuf)
        else
            base;
        var g = buildGeneProgram(name, dna);
        // class expression bias
        const bias = expressionBias(ct, name);
        g.expression = @max(0.05, @min(3.5, g.expression * bias));
        gt.genes[i] = g;
    }
    gt.phenotype = phenotypeFromGenes(&gt.genes);
    applyClassNudge(ct, &gt.phenotype);
    gt.composite_spin = gt.phenotype.composite_spin;
    gt.composite_charge = gt.phenotype.composite_charge;
    return gt;
}

pub fn buildNeuronGenotype(unit_id: u32, diversity: bool) NeuronGenotype {
    // default Pyr template without cell-type ORF overrides (generic CHANNEL_GENE_ORFS)
    var gt: NeuronGenotype = .{
        .unit_id = unit_id,
        .cell_type = .pyr,
        .synapse_sign = 1,
    };
    const names = [_]GeneName{ .scn, .kcn, .cacna, .leak };
    var mutbuf: [32]u8 = undefined;
    var i: usize = 0;
    while (i < 4) : (i += 1) {
        const name = names[i];
        const base = baseOrf(name);
        const dna = if (diversity)
            mutateOrf(base, unit_id, @intCast(i), &mutbuf)
        else
            base;
        gt.genes[i] = buildGeneProgram(name, dna);
    }
    gt.phenotype = phenotypeFromGenes(&gt.genes);
    gt.composite_spin = gt.phenotype.composite_spin;
    gt.composite_charge = gt.phenotype.composite_charge;
    return gt;
}

pub fn phenotypeToUnitParams(ph: Phenotype) bio_probe.UnitParams {
    return .{
        .d_eff = ph.d_eff,
        .fire_thr = ph.fire_threshold,
        .ref_steps = @intFromFloat(@round(ph.refractory_steps)),
        .adapt_gain = ph.adapt_gain,
        .adapt_decay = ph.adapt_decay,
        .adapt_step = ph.adapt_step,
        .fi_stim = ph.fi_stim,
    };
}

pub fn selfTest() bool {
    if (!codon.selfTest()) return false;

    // Base gene programs (no diversity) — match Python CHANNEL_GENE_ORFS
    const scn = buildGeneProgram(.scn, ORF_SCN);
    if (@abs(scn.spin) > 1e-9) return false;
    if (scn.expression < 1.42 or scn.expression > 1.45) return false;
    if (scn.charge_balance != 2) return false;

    const kcn = buildGeneProgram(.kcn, ORF_KCN);
    if (kcn.spin > -0.20 or kcn.spin < -0.24) return false;

    const gt = buildCellTypeGenotype(0, .pyr, false);
    // Python Pyr spin ~0.1045, charge ~0.123, ref ~13.44, fi ~0.597, d_eff ~13.375
    if (@abs(gt.composite_spin - 0.1044915460095421) > 0.02) return false;
    if (@abs(gt.phenotype.refractory_steps - 13.4419) > 0.5) return false;
    if (@abs(gt.phenotype.fi_stim - 0.59726) > 0.05) return false;
    if (@abs(gt.phenotype.d_eff - 13.3754) > 0.15) return false;

    const pv = buildCellTypeGenotype(0, .pv, false);
    if (pv.phenotype.refractory_steps >= gt.phenotype.refractory_steps) return false; // PV faster
    if (pv.synapse_sign >= 0) return false;
    return true;
}
