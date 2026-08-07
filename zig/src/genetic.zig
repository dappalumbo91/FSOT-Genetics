//! Genetic synaptic weights — full twin of brain_architecture._build_weight_matrix.
//! Spins/charges from codon genotypes only. Seeds only. No free-fit shortcuts.

const seeds = @import("seeds.zig");
const cell_types = @import("cell_types.zig");
const genotype = @import("genotype.zig");
const network = @import("network.zig");

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

pub fn geometricScaleDist(dist: usize) f64 {
    // φ · dist^(-1/π) with dist >= 1
    const d: f64 = @floatFromInt(if (dist < 1) 1 else dist);
    return seeds.phi * @exp(@log(d) * (-1.0 / seeds.pi));
}

pub fn electrostaticTerm(q_i: f64, q_j: f64) f64 {
    return -q_i * q_j * seeds.e;
}

/// Python _fsot_pair_weight(gi, gj, dist)
pub fn fsotPairWeight(spin_i: f64, spin_j: f64, charge_i: f64, charge_j: f64, dist: usize) f64 {
    const base = trinaryPairInteraction(spin_i, spin_j);
    const geom = geometricScaleDist(dist);
    const elec = electrostaticTerm(charge_i, charge_j);
    const d: f64 = @floatFromInt(if (dist < 1) 1 else dist);
    const env = d / (d + seeds.pi * seeds.e);
    return geom * (base + 0.15 * elec) * (0.35 + 0.65 * env);
}

pub fn motifGain(pre: cell_types.CellType, post: cell_types.CellType) f64 {
    const pre_s = cell_types.specOf(pre).sign;
    const post_s = cell_types.specOf(post).sign;
    if (pre == .vip and post_s < 0) return 0.30; // gain_vip_i
    if (pre_s > 0 and post_s > 0) return 0.085; // gain_ee
    if (pre_s > 0 and post_s < 0) return 0.55; // gain_ei
    if (pre_s < 0 and post_s > 0) return 0.42; // gain_ie
    if (pre_s < 0 and post_s < 0) return 0.22; // gain_ii
    return 0.3;
}

pub const Projection = struct {
    src: u8, // RegionId as u8
    dst: u8,
    density: f64,
    strength: f64,
    /// hash input like Python proj.kind string
    kind_hash: u32,
};

/// Default projections (DEFAULT_PROJECTIONS)
pub const DEFAULT_PROJECTIONS = [_]Projection{
    .{ .src = 0, .dst = 1, .density = 0.35, .strength = 0.55, .kind_hash = kindHash("feedforward") },
    .{ .src = 1, .dst = 2, .density = 0.25, .strength = 0.45, .kind_hash = kindHash("feedforward") },
    .{ .src = 2, .dst = 1, .density = 0.12, .strength = 0.25, .kind_hash = kindHash("feedback") },
    .{ .src = 2, .dst = 3, .density = 0.20, .strength = 0.40, .kind_hash = kindHash("feedforward") },
    .{ .src = 3, .dst = 2, .density = 0.15, .strength = 0.30, .kind_hash = kindHash("feedback") },
    .{ .src = 1, .dst = 0, .density = 0.10, .strength = 0.20, .kind_hash = kindHash("feedback") },
};

fn kindHash(s: []const u8) u32 {
    // Python: sum((i+1)*ord(c) for i,c in enumerate(proj.kind))
    var h: u32 = 0;
    var i: usize = 0;
    while (i < s.len) : (i += 1) {
        h += @as(u32, @intCast(i + 1)) * @as(u32, s[i]);
    }
    return h;
}

const UnitMeta = struct {
    region: u8,
    local_id: usize,
    ct: cell_types.CellType,
    sign: i8,
    spin: f64,
    charge: f64,
};

/// Full W builder — Python parity path.
pub fn wireFromGenotypes(
    W: []f64,
    max_n: usize,
    n: usize,
    gts: []const genotype.NeuronGenotype,
    region_of: []const u8,
    region_local: []const usize,
    local_syn_scale: f64,
) void {
    // clear
    var i: usize = 0;
    while (i < n) : (i += 1) {
        var j: usize = 0;
        while (j < n) : (j += 1) W[i * max_n + j] = 0;
    }

    var meta: [network.MAX_N]UnitMeta = undefined;
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

    // k_ff ~ e ≈ 3, k_fb ~ phi/2 ≈ 1
    const k_ff: usize = @max(@as(usize, 2), @as(usize, @intFromFloat(@round(seeds.e))));
    const k_fb: usize = @max(@as(usize, 1), @as(usize, @intFromFloat(@round(seeds.phi / 2.0))));
    const recip_scale = 1.0 / seeds.phi;

    const setEdge = struct {
        fn call(
            Wmat: []f64,
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
                    gain *= seeds.phi / (seeds.phi + 1.0 / seeds.phi);
                } else {
                    gain *= 1.0 / seeds.phi;
                }
            }
            const polarity: f64 = @floatFromInt(mq.sign);
            Wmat[post * maxn + pre] = w * gain * polarity;
        }
    }.call;

    // --- local per region ---
    var reg: u8 = 0;
    while (reg < 4) : (reg += 1) {
        // collect ids in region ordered by local_id
        var ids: [network.MAX_N]usize = undefined;
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
        var e_ids: [network.MAX_N]usize = undefined;
        var i_ids: [network.MAX_N]usize = undefined;
        var n_e: usize = 0;
        var n_i: usize = 0;
        a = 0;
        while (a < n_ids) : (a += 1) {
            const uid = ids[a];
            if (meta[uid].sign > 0) {
                e_ids[n_e] = uid;
                n_e += 1;
            } else {
                i_ids[n_i] = uid;
                n_i += 1;
            }
        }

        // E→I full bipartite
        var pe: usize = 0;
        while (pe < n_e) : (pe += 1) {
            var po: usize = 0;
            while (po < n_i) : (po += 1) {
                setEdge(W, max_n, meta[0..n], i_ids[po], e_ids[pe], false);
            }
        }
        // I→E full bipartite
        var pi: usize = 0;
        while (pi < n_i) : (pi += 1) {
            pe = 0;
            while (pe < n_e) : (pe += 1) {
                setEdge(W, max_n, meta[0..n], e_ids[pe], i_ids[pi], false);
            }
        }
        // I→I sparse directed
        pi = 0;
        while (pi < n_i) : (pi += 1) {
            const pre = i_ids[pi];
            var higher: [network.MAX_N]usize = undefined;
            var lower: [network.MAX_N]usize = undefined;
            var nh: usize = 0;
            var nl: usize = 0;
            var q: usize = 0;
            while (q < n_i) : (q += 1) {
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
            while (t < k_ff and t < nh) : (t += 1) {
                setEdge(W, max_n, meta[0..n], higher[t], pre, false);
            }
            const kfb_i = @max(@as(usize, 1), k_fb);
            t = 0;
            while (t < kfb_i and t < nl) : (t += 1) {
                setEdge(W, max_n, meta[0..n], lower[t], pre, false);
            }
        }
        // E→E sparse directed + ff_boost
        pe = 0;
        while (pe < n_e) : (pe += 1) {
            const pre = e_ids[pe];
            var higher: [network.MAX_N]usize = undefined;
            var lower: [network.MAX_N]usize = undefined;
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
            while (t < k_ff and t < nh) : (t += 1) {
                setEdge(W, max_n, meta[0..n], higher[t], pre, true);
            }
            t = 0;
            while (t < k_fb and t < nl) : (t += 1) {
                setEdge(W, max_n, meta[0..n], lower[t], pre, true);
            }
        }
    }

    // --- long-range projections ---
    for (DEFAULT_PROJECTIONS) |proj| {
        var src_e: [network.MAX_N]usize = undefined;
        var dst_e: [network.MAX_N]usize = undefined;
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
        const k = @max(@as(usize, 1), @as(usize, @intFromFloat(@floor(proj.density * @as(f64, @floatFromInt(nd))))));
        var si: usize = 0;
        while (si < ns) : (si += 1) {
            const pre = src_e[si];
            var j: usize = 0;
            while (j < k) : (j += 1) {
                const post = dst_e[(si * 7 + j * 3 + proj.kind_hash % 5) % nd];
                if (post == pre) continue;
                const dist = 8 + (if (si > j) si - j else j - si);
                const w = fsotPairWeight(meta[post].spin, meta[pre].spin, meta[post].charge, meta[pre].charge, dist);
                W[post * max_n + pre] += w * proj.strength;
            }
        }
    }

    // Reciprocal attenuation same-sign pairs (1/φ²)
    const scale = recip_scale * recip_scale;
    i = 0;
    while (i < n) : (i += 1) {
        var j: usize = i + 1;
        while (j < n) : (j += 1) {
            if (meta[i].sign != meta[j].sign) continue;
            const a = W[i * max_n + j];
            const b = W[j * max_n + i];
            if (a == 0 or b == 0) continue;
            const aa = if (a < 0) -a else a;
            const bb = if (b < 0) -b else b;
            if (aa >= bb) {
                W[j * max_n + i] *= scale;
            } else {
                W[i * max_n + j] *= scale;
            }
        }
    }

    // Drop near-zero: thr = mean_abs / (φ·e)
    var sum_abs: f64 = 0;
    var n_nz: f64 = 0;
    i = 0;
    while (i < n) : (i += 1) {
        var j: usize = 0;
        while (j < n) : (j += 1) {
            const w = W[i * max_n + j];
            if (w != 0) {
                sum_abs += if (w < 0) -w else w;
                n_nz += 1;
            }
        }
    }
    if (n_nz > 0) {
        const mean_abs = sum_abs / n_nz;
        const thr = mean_abs / (seeds.phi * seeds.e);
        i = 0;
        while (i < n) : (i += 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) {
                const w = W[i * max_n + j];
                const aw = if (w < 0) -w else w;
                if (aw < thr) W[i * max_n + j] = 0;
            }
        }
    }

    // Renormalize mean |W| → local_syn_scale
    sum_abs = 0;
    n_nz = 0;
    i = 0;
    while (i < n) : (i += 1) {
        var j: usize = 0;
        while (j < n) : (j += 1) {
            const w = W[i * max_n + j];
            if (w != 0) {
                sum_abs += if (w < 0) -w else w;
                n_nz += 1;
            }
        }
    }
    if (n_nz > 0 and sum_abs > 1e-12) {
        const mean_abs = sum_abs / n_nz;
        const sc = local_syn_scale / mean_abs;
        i = 0;
        while (i < n) : (i += 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) W[i * max_n + j] *= sc;
        }
    }
}

pub fn selfTest() bool {
    const w = trinaryPairInteraction(1.0, -1.0);
    if (@abs(w + seeds.e) > 1e-9) return false;
    if (geometricScaleDist(3) <= 0) return false;
    if (motifGain(.pyr, .pv) < 0.5) return false;
    const g0 = genotype.buildCellTypeGenotype(0, .pyr, false);
    const g1 = genotype.buildCellTypeGenotype(1, .pv, false);
    const pw = fsotPairWeight(g0.composite_spin, g1.composite_spin, g0.composite_charge, g1.composite_charge, 2);
    if (pw != pw) return false;
    return true;
}
