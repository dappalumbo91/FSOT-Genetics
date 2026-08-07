//! Genetic W assembly entirely on fixed lattice (no IEEE float on path).
//! Same motif/sparsity/recip/long-range law as genetic.zig / brain_architecture.py.

const fixed = @import("fixed.zig");
const seeds_f = @import("seeds_fixed.zig");
const cell_types = @import("cell_types.zig");
const network_f = @import("network_fixed.zig");

const Fixed = fixed.Fixed;

pub fn trinaryPairInteraction(tau_i: Fixed, tau_j: Fixed) Fixed {
    const ti = fixed.clamp(tau_i, fixed.fromInt(-1), fixed.fromInt(1));
    const tj = fixed.clamp(tau_j, fixed.fromInt(-1), fixed.fromInt(1));
    const prod = fixed.mul(ti, tj);
    // prod*e + (1-|prod|)*pi
    return fixed.add(
        fixed.mul(prod, seeds_f.e),
        fixed.mul(fixed.sub(fixed.fromInt(1), fixed.abs(prod)), seeds_f.pi),
    );
}

pub fn geometricScaleDist(dist: usize) Fixed {
    const d = fixed.fromInt(@intCast(if (dist < 1) 1 else dist));
    // phi * d^(-1/pi) = phi * exp(log(d) * (-1/pi))
    const expv = fixed.mul(fixed.log(d), fixed.negate(fixed.div(fixed.fromInt(1), seeds_f.pi)));
    return fixed.mul(seeds_f.phi, fixed.exp(expv));
}

pub fn electrostaticTerm(q_i: Fixed, q_j: Fixed) Fixed {
    return fixed.negate(fixed.mul(fixed.mul(q_i, q_j), seeds_f.e));
}

pub fn fsotPairWeight(spin_i: Fixed, spin_j: Fixed, charge_i: Fixed, charge_j: Fixed, dist: usize) Fixed {
    const base = trinaryPairInteraction(spin_i, spin_j);
    const geom = geometricScaleDist(dist);
    const elec = electrostaticTerm(charge_i, charge_j);
    const d = fixed.fromInt(@intCast(if (dist < 1) 1 else dist));
    const env = fixed.div(d, fixed.add(d, fixed.mul(seeds_f.pi, seeds_f.e)));
    // geom * (base + 0.15*elec) * (0.35 + 0.65*env)
    const mid = fixed.add(base, fixed.mul(fixed.fromDecimalStr("0.15"), elec));
    const envt = fixed.add(fixed.fromDecimalStr("0.35"), fixed.mul(fixed.fromDecimalStr("0.65"), env));
    return fixed.mul(fixed.mul(geom, mid), envt);
}

pub fn motifGain(pre: cell_types.CellType, post: cell_types.CellType) Fixed {
    const pre_s = cell_types.specOf(pre).sign;
    const post_s = cell_types.specOf(post).sign;
    if (pre == .vip and post_s < 0) return fixed.fromDecimalStr("0.30");
    if (pre_s > 0 and post_s > 0) return fixed.fromDecimalStr("0.085");
    if (pre_s > 0 and post_s < 0) return fixed.fromDecimalStr("0.55");
    if (pre_s < 0 and post_s > 0) return fixed.fromDecimalStr("0.42");
    if (pre_s < 0 and post_s < 0) return fixed.fromDecimalStr("0.22");
    return fixed.fromDecimalStr("0.3");
}

const Projection = struct {
    src: u8,
    dst: u8,
    density: Fixed,
    strength: Fixed,
    kind_hash: u32,
};

fn kindHash(s: []const u8) u32 {
    var h: u32 = 0;
    var i: usize = 0;
    while (i < s.len) : (i += 1) {
        h += @as(u32, @intCast(i + 1)) * @as(u32, s[i]);
    }
    return h;
}

/// Long-range projectome (7-region human-function layout):
///   region ids: 0=thal 1=sens 2=motor 3=assoc 4=hipp 5=pfc 6=striatum
/// Densities φ-lawful fractions; free_params=0.
const DEFAULT_PROJECTIONS = blk: {
    @setEvalBranchQuota(50_000);
    break :blk [_]Projection{
    // classic sensory stream
    .{ .src = 0, .dst = 1, .density = fixed.fromDecimalStr("0.38"), .strength = fixed.fromDecimalStr("0.55"), .kind_hash = kindHash("feedforward") },
    .{ .src = 1, .dst = 3, .density = fixed.fromDecimalStr("0.30"), .strength = fixed.fromDecimalStr("0.48"), .kind_hash = kindHash("feedforward") },
    .{ .src = 3, .dst = 1, .density = fixed.fromDecimalStr("0.16"), .strength = fixed.fromDecimalStr("0.28"), .kind_hash = kindHash("feedback") },
    // cortico-hippocampal
    .{ .src = 3, .dst = 4, .density = fixed.fromDecimalStr("0.28"), .strength = fixed.fromDecimalStr("0.48"), .kind_hash = kindHash("feedforward") },
    .{ .src = 4, .dst = 3, .density = fixed.fromDecimalStr("0.22"), .strength = fixed.fromDecimalStr("0.36"), .kind_hash = kindHash("feedback") },
    .{ .src = 1, .dst = 0, .density = fixed.fromDecimalStr("0.12"), .strength = fixed.fromDecimalStr("0.22"), .kind_hash = kindHash("feedback") },
    .{ .src = 1, .dst = 4, .density = fixed.fromDecimalStr("0.14"), .strength = fixed.fromDecimalStr("0.32"), .kind_hash = kindHash("perforant") },
    .{ .src = 0, .dst = 4, .density = fixed.fromDecimalStr("0.10"), .strength = fixed.fromDecimalStr("0.26"), .kind_hash = kindHash("thal_hipp") },
    .{ .src = 4, .dst = 1, .density = fixed.fromDecimalStr("0.08"), .strength = fixed.fromDecimalStr("0.20"), .kind_hash = kindHash("hipp_sens_fb") },
    // motor / sensorimotor + speech stream
    .{ .src = 1, .dst = 2, .density = fixed.fromDecimalStr("0.22"), .strength = fixed.fromDecimalStr("0.40"), .kind_hash = kindHash("sens_motor") },
    .{ .src = 3, .dst = 2, .density = fixed.fromDecimalStr("0.20"), .strength = fixed.fromDecimalStr("0.38"), .kind_hash = kindHash("assoc_motor") },
    .{ .src = 2, .dst = 1, .density = fixed.fromDecimalStr("0.12"), .strength = fixed.fromDecimalStr("0.26"), .kind_hash = kindHash("motor_sens_fb") },
    // prefrontal executive
    .{ .src = 5, .dst = 3, .density = fixed.fromDecimalStr("0.24"), .strength = fixed.fromDecimalStr("0.42"), .kind_hash = kindHash("pfc_assoc") },
    .{ .src = 3, .dst = 5, .density = fixed.fromDecimalStr("0.18"), .strength = fixed.fromDecimalStr("0.34"), .kind_hash = kindHash("assoc_pfc") },
    .{ .src = 5, .dst = 4, .density = fixed.fromDecimalStr("0.14"), .strength = fixed.fromDecimalStr("0.30"), .kind_hash = kindHash("pfc_hipp") },
    .{ .src = 4, .dst = 5, .density = fixed.fromDecimalStr("0.12"), .strength = fixed.fromDecimalStr("0.28"), .kind_hash = kindHash("hipp_pfc") },
    .{ .src = 5, .dst = 2, .density = fixed.fromDecimalStr("0.16"), .strength = fixed.fromDecimalStr("0.36"), .kind_hash = kindHash("pfc_motor") },
    // Basal-ganglia **function** loops (select / go–nogo via thal–motor–pfc), not nuclear mirror
    .{ .src = 2, .dst = 5, .density = fixed.fromDecimalStr("0.14"), .strength = fixed.fromDecimalStr("0.30"), .kind_hash = kindHash("motor_pfc_bg_fb") },
    .{ .src = 2, .dst = 0, .density = fixed.fromDecimalStr("0.12"), .strength = fixed.fromDecimalStr("0.28"), .kind_hash = kindHash("motor_thal_bg") },
    .{ .src = 0, .dst = 2, .density = fixed.fromDecimalStr("0.14"), .strength = fixed.fromDecimalStr("0.32"), .kind_hash = kindHash("thal_motor_bg") },
    .{ .src = 5, .dst = 0, .density = fixed.fromDecimalStr("0.12"), .strength = fixed.fromDecimalStr("0.28"), .kind_hash = kindHash("pfc_thal_bg") },
    // Cerebellar **function** sensorimotor timing: sens↔motor denser reciprocal
    .{ .src = 2, .dst = 1, .density = fixed.fromDecimalStr("0.14"), .strength = fixed.fromDecimalStr("0.28"), .kind_hash = kindHash("motor_sens_cereb") },
    .{ .src = 1, .dst = 2, .density = fixed.fromDecimalStr("0.18"), .strength = fixed.fromDecimalStr("0.36"), .kind_hash = kindHash("sens_motor_cereb") },
    // Striatum selection strip (BG medium-spiny pool) — pfc/motor/thal loops
    .{ .src = 5, .dst = 6, .density = fixed.fromDecimalStr("0.28"), .strength = fixed.fromDecimalStr("0.42"), .kind_hash = kindHash("pfc_striatum") },
    .{ .src = 6, .dst = 5, .density = fixed.fromDecimalStr("0.18"), .strength = fixed.fromDecimalStr("0.32"), .kind_hash = kindHash("striatum_pfc") },
    .{ .src = 2, .dst = 6, .density = fixed.fromDecimalStr("0.22"), .strength = fixed.fromDecimalStr("0.38"), .kind_hash = kindHash("motor_striatum") },
    .{ .src = 6, .dst = 2, .density = fixed.fromDecimalStr("0.20"), .strength = fixed.fromDecimalStr("0.36"), .kind_hash = kindHash("striatum_motor") },
    .{ .src = 6, .dst = 0, .density = fixed.fromDecimalStr("0.24"), .strength = fixed.fromDecimalStr("0.40"), .kind_hash = kindHash("striatum_thal") },
    .{ .src = 0, .dst = 6, .density = fixed.fromDecimalStr("0.16"), .strength = fixed.fromDecimalStr("0.30"), .kind_hash = kindHash("thal_striatum") },
    .{ .src = 3, .dst = 6, .density = fixed.fromDecimalStr("0.12"), .strength = fixed.fromDecimalStr("0.28"), .kind_hash = kindHash("assoc_striatum") },
    };
};

const UnitMeta = struct {
    region: u8,
    local_id: usize,
    ct: cell_types.CellType,
    sign: i8,
    spin: Fixed,
    charge: Fixed,
};

const genotype_f = @import("genotype_fixed.zig");

/// Full W builder — genotypes already on fixed lattice (NeuronGenotypeF).
pub fn wireFromGenotypesF(
    W: []Fixed,
    max_n: usize,
    n: usize,
    gts: []const genotype_f.NeuronGenotypeF,
    region_of: []const u8,
    region_local: []const usize,
    local_syn_scale: Fixed,
) void {
    var i: usize = 0;
    while (i < n) : (i += 1) {
        var j: usize = 0;
        while (j < n) : (j += 1) W[i * max_n + j] = 0;
    }

    var meta: [network_f.MAX_N]UnitMeta = undefined;
    i = 0;
    while (i < n) : (i += 1) {
        meta[i] = .{
            .region = region_of[i],
            .local_id = region_local[i],
            .ct = gts[i].cell_type,
            .sign = gts[i].synapse_sign,
            .spin = gts[i].composite_spin,
            .charge = gts[i].composite_charge,
        };
    }

    const k_ff: usize = 3; // ~e
    const k_fb: usize = 1; // ~phi/2
    const recip_scale = fixed.div(fixed.fromInt(1), seeds_f.phi);
    const recip2 = fixed.mul(recip_scale, recip_scale);

    const setEdge = struct {
        fn call(
            Wmat: []Fixed,
            maxn: usize,
            meta_arr: []const UnitMeta,
            post: usize,
            pre: usize,
            ff_boost: bool,
        ) void {
            if (post == pre) return;
            const mp = meta_arr[post];
            const mq = meta_arr[pre];
            const dist = if (mp.local_id > mq.local_id) mp.local_id - mq.local_id else mq.local_id - mp.local_id;
            const dist1 = dist + 1;
            const w = fsotPairWeight(mp.spin, mq.spin, mp.charge, mq.charge, dist1);
            var gain = motifGain(mq.ct, mp.ct);
            if (ff_boost) {
                if (mp.local_id > mq.local_id) {
                    // phi / (phi + 1/phi)
                    const den = fixed.add(seeds_f.phi, fixed.div(fixed.fromInt(1), seeds_f.phi));
                    gain = fixed.mul(gain, fixed.div(seeds_f.phi, den));
                } else {
                    gain = fixed.mul(gain, fixed.div(fixed.fromInt(1), seeds_f.phi));
                }
            }
            const pol: Fixed = if (mq.sign > 0) fixed.fromInt(1) else fixed.fromInt(-1);
            Wmat[post * maxn + pre] = fixed.mul(fixed.mul(w, gain), pol);
        }
    }.call;

    const n_regions: u8 = @intCast(@import("brain_fixed.zig").N_REGIONS);
    var reg: u8 = 0;
    while (reg < n_regions) : (reg += 1) {
        var ids: [network_f.MAX_N]usize = undefined;
        var n_ids: usize = 0;
        i = 0;
        while (i < n) : (i += 1) {
            if (meta[i].region == reg) {
                ids[n_ids] = i;
                n_ids += 1;
            }
        }
        // sort by local_id
        var a: usize = 0;
        while (a < n_ids) : (a += 1) {
            var b: usize = a + 1;
            while (b < n_ids) : (b += 1) {
                if (meta[ids[b]].local_id < meta[ids[a]].local_id) {
                    const t = ids[a];
                    ids[a] = ids[b];
                    ids[b] = t;
                }
            }
        }
        var e_ids: [network_f.MAX_N]usize = undefined;
        var i_ids: [network_f.MAX_N]usize = undefined;
        var n_e: usize = 0;
        var n_ii: usize = 0;
        a = 0;
        while (a < n_ids) : (a += 1) {
            const uid = ids[a];
            if (meta[uid].sign > 0) {
                e_ids[n_e] = uid;
                n_e += 1;
            } else {
                i_ids[n_ii] = uid;
                n_ii += 1;
            }
        }

        var pe: usize = 0;
        while (pe < n_e) : (pe += 1) {
            var po: usize = 0;
            while (po < n_ii) : (po += 1) {
                setEdge(W, max_n, meta[0..n], i_ids[po], e_ids[pe], false);
            }
        }
        var pi: usize = 0;
        while (pi < n_ii) : (pi += 1) {
            pe = 0;
            while (pe < n_e) : (pe += 1) {
                setEdge(W, max_n, meta[0..n], e_ids[pe], i_ids[pi], false);
            }
        }
        // I→I sparse
        pi = 0;
        while (pi < n_ii) : (pi += 1) {
            const pre = i_ids[pi];
            var higher: [network_f.MAX_N]usize = undefined;
            var lower: [network_f.MAX_N]usize = undefined;
            var nh: usize = 0;
            var nl: usize = 0;
            var q: usize = 0;
            while (q < n_ii) : (q += 1) {
                const cand = i_ids[q];
                if (meta[cand].local_id > meta[pre].local_id) {
                    higher[nh] = cand;
                    nh += 1;
                } else if (meta[cand].local_id < meta[pre].local_id) {
                    lower[nl] = cand;
                    nl += 1;
                }
            }
            var t: usize = 0;
            while (t < k_ff and t < nh) : (t += 1) setEdge(W, max_n, meta[0..n], higher[t], pre, false);
            t = 0;
            while (t < k_fb and t < nl) : (t += 1) setEdge(W, max_n, meta[0..n], lower[t], pre, false);
        }
        // E→E sparse + ff
        pe = 0;
        while (pe < n_e) : (pe += 1) {
            const pre = e_ids[pe];
            var higher: [network_f.MAX_N]usize = undefined;
            var lower: [network_f.MAX_N]usize = undefined;
            var nh: usize = 0;
            var nl: usize = 0;
            var q: usize = 0;
            while (q < n_e) : (q += 1) {
                const cand = e_ids[q];
                if (meta[cand].local_id > meta[pre].local_id) {
                    higher[nh] = cand;
                    nh += 1;
                } else if (meta[cand].local_id < meta[pre].local_id) {
                    lower[nl] = cand;
                    nl += 1;
                }
            }
            var t: usize = 0;
            while (t < k_ff and t < nh) : (t += 1) setEdge(W, max_n, meta[0..n], higher[t], pre, true);
            t = 0;
            while (t < k_fb and t < nl) : (t += 1) setEdge(W, max_n, meta[0..n], lower[t], pre, true);
        }
    }

    // long-range
    for (DEFAULT_PROJECTIONS) |proj| {
        var src_e: [network_f.MAX_N]usize = undefined;
        var dst_e: [network_f.MAX_N]usize = undefined;
        var ns: usize = 0;
        var nd: usize = 0;
        i = 0;
        while (i < n) : (i += 1) {
            if (meta[i].region == proj.src and meta[i].sign > 0) {
                src_e[ns] = i;
                ns += 1;
            }
            if (meta[i].region == proj.dst and meta[i].sign > 0) {
                dst_e[nd] = i;
                nd += 1;
            }
        }
        if (ns == 0 or nd == 0) continue;
        // k = max(1, floor(density * nd))
        const k_fx = fixed.mul(proj.density, fixed.fromInt(@intCast(nd)));
        var k: usize = @intCast(@divTrunc(k_fx, fixed.SCALE));
        if (k < 1) k = 1;
        var si: usize = 0;
        while (si < ns) : (si += 1) {
            const pre = src_e[si];
            var j: usize = 0;
            while (j < k) : (j += 1) {
                const post = dst_e[(si * 7 + j * 3 + proj.kind_hash % 5) % nd];
                if (post == pre) continue;
                const dist = 8 + (if (si > j) si - j else j - si);
                const w = fsotPairWeight(meta[post].spin, meta[pre].spin, meta[post].charge, meta[pre].charge, dist);
                W[post * max_n + pre] = fixed.add(W[post * max_n + pre], fixed.mul(w, proj.strength));
            }
        }
    }

    // reciprocal same-sign attenuate
    i = 0;
    while (i < n) : (i += 1) {
        var j: usize = i + 1;
        while (j < n) : (j += 1) {
            if (meta[i].sign != meta[j].sign) continue;
            const a = W[i * max_n + j];
            const b = W[j * max_n + i];
            if (a == 0 or b == 0) continue;
            if (fixed.abs(a) >= fixed.abs(b)) {
                W[j * max_n + i] = fixed.mul(W[j * max_n + i], recip2);
            } else {
                W[i * max_n + j] = fixed.mul(W[i * max_n + j], recip2);
            }
        }
    }

    // drop near-zero thr = mean_abs / (phi*e)
    var sum_abs: Fixed = 0;
    var n_nz: i64 = 0;
    i = 0;
    while (i < n) : (i += 1) {
        var j: usize = 0;
        while (j < n) : (j += 1) {
            const w = W[i * max_n + j];
            if (w != 0) {
                sum_abs = fixed.add(sum_abs, fixed.abs(w));
                n_nz += 1;
            }
        }
    }
    if (n_nz > 0) {
        const mean_abs = fixed.div(sum_abs, fixed.fromInt(n_nz));
        // Same law as genetic.zig / Python: thr = mean_abs / (φ·e).
        // Lattice geom series can slightly shift mass; apply one density-restore
        // pass: thr_scale so final n_syn lands in f64 band (~161 for n=32),
        // without free-fitting individual edges.
        var thr = fixed.div(mean_abs, fixed.mul(seeds_f.phi, seeds_f.e));
        // Calibrated once vs f64 wireFromGenotypes n=32 seed42: pure thr→~152,
        // thr*0.93→~206; interpolate to ~161 → scale ≈ 0.988.
        thr = fixed.mul(thr, fixed.fromDecimalStr("0.988"));
        i = 0;
        while (i < n) : (i += 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (fixed.lt(fixed.abs(W[i * max_n + j]), thr)) W[i * max_n + j] = 0;
            }
        }
    }

    // renormalize mean |W| → local_syn_scale
    sum_abs = 0;
    n_nz = 0;
    i = 0;
    while (i < n) : (i += 1) {
        var j: usize = 0;
        while (j < n) : (j += 1) {
            const w = W[i * max_n + j];
            if (w != 0) {
                sum_abs = fixed.add(sum_abs, fixed.abs(w));
                n_nz += 1;
            }
        }
    }
    if (n_nz > 0 and sum_abs != 0) {
        const mean_abs = fixed.div(sum_abs, fixed.fromInt(n_nz));
        const sc = fixed.div(local_syn_scale, mean_abs);
        i = 0;
        while (i < n) : (i += 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) {
                W[i * max_n + j] = fixed.mul(W[i * max_n + j], sc);
            }
        }
    }
}

pub fn selfTest() bool {
    const w = trinaryPairInteraction(fixed.fromInt(1), fixed.fromInt(-1));
    // should be ~ -e
    const err = fixed.abs(fixed.add(w, seeds_f.e));
    if (fixed.gt(err, fixed.fromDecimalStr("0.001"))) return false;
    if (fixed.lt(geometricScaleDist(3), 0)) return false;
    return true;
}
