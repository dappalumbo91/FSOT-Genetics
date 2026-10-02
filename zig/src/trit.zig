//! FSOT trinary substrate — Zig bare-metal oracle (matches Python trinary_substrate).
//! Ontology: T = {-1, 0, +1}. Integer packing is transport only.

pub const Trit = i8; // -1 | 0 | +1

pub fn asTrit(x: i32) Trit {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

pub fn neg(t: Trit) Trit {
    return asTrit(-@as(i32, t));
}

pub fn pair(a: Trit, b: Trit) Trit {
    return asTrit(@as(i32, a) * @as(i32, b));
}

pub fn sumSat(a: Trit, b: Trit) Trit {
    return asTrit(@as(i32, a) + @as(i32, b));
}

pub fn consensus(a: Trit, b: Trit) Trit {
    return if (a == b) a else 0;
}

pub fn fromS(s: f32, lo: f32, hi: f32) Trit {
    if (s < lo) return -1;
    if (s > hi) return 1;
    return 0;
}

/// A,G → +1; C,T,U → −1 (U is a pyrimidine and maps exactly like T;
/// same as codon_core `nt_primary` in FSOT-Genetics Rust, TRIT_SPEC §3).
pub fn basePrimary(b: u8) Trit {
    return switch (b) {
        'A', 'a', 'G', 'g' => 1,
        'C', 'c', 'T', 't', 'U', 'u' => -1,
        else => 0,
    };
}

pub fn codonPrimary(c0: u8, c1: u8, c2: u8) [3]Trit {
    return .{ basePrimary(c0), basePrimary(c1), basePrimary(c2) };
}

// --- T1 packing: 2 bits/trit  00=0, 01=+1, 11=-1 ---
// LAYOUT LABEL: T1 is a legacy sign/magnitude IN-MEMORY compute layout
// (bit0 = non-zero, bit1 = sign; 10 invalid). It is NOT the interchange format.
// The canonical storage/wire layout (TRIT_SPEC §2a, used by FSOT-GPU /
// FSOT-Quantum fsot_lib/trinary.py `pack_u64`, codon_core and the CUDA kernels)
// is code = t + 1:  -1 -> 00, 0 -> 01, +1 -> 10, 11 invalid.
// The same bit pattern means different trits in the two layouts (01 = +1 in T1,
// 0 in canonical), so convert with toCanonical/fromCanonical (single trit) or
// t1WordToCanonical/canonicalWordToT1 (32-trit u64) at every file/wire boundary.
// T1 is kept unchanged so existing T1 data stays readable.

const t1_pos: u8 = 0b01;
const t1_neg: u8 = 0b11;
const t1_zero: u8 = 0b00;

pub fn packT1(t: Trit) u8 {
    if (t > 0) return t1_pos;
    if (t < 0) return t1_neg;
    return t1_zero;
}

pub fn unpackT1(bits: u8) ?Trit {
    return switch (bits & 0b11) {
        t1_pos => 1,
        t1_neg => -1,
        t1_zero => 0,
        else => null,
    };
}

/// Canonical v1 blob tag (TRIT_SPEC §2a): prefix a stored canonical trit blob with this byte.
pub const canonical_v1_tag: u8 = 0xF1;

/// Trit -> canonical 2-bit code (code = t + 1).
pub fn toCanonical(t: Trit) u8 {
    if (t > 0) return 2;
    if (t < 0) return 0;
    return 1;
}

/// Canonical 2-bit code -> trit; code 3 (bits 11) is invalid and returns null.
pub fn fromCanonical(code: u8) ?Trit {
    return switch (code & 0b11) {
        0 => -1,
        1 => 0,
        2 => 1,
        else => null,
    };
}

/// Re-encode n (<= 32) T1 lanes of a u64 as canonical lanes. Null on an invalid T1 lane (10).
pub fn t1WordToCanonical(t1: u64, n: u8) ?u64 {
    var out: u64 = 0;
    var i: u8 = 0;
    while (i < n and i < 32) : (i += 1) {
        const bits: u8 = @truncate(t1 >> @intCast(2 * @as(u32, i)));
        const t = unpackT1(bits) orelse return null;
        out |= @as(u64, toCanonical(t)) << @intCast(2 * @as(u32, i));
    }
    return out;
}

/// Re-encode n (<= 32) canonical lanes of a u64 as T1 lanes. Null on an invalid canonical lane (11).
pub fn canonicalWordToT1(canon: u64, n: u8) ?u64 {
    var out: u64 = 0;
    var i: u8 = 0;
    while (i < n and i < 32) : (i += 1) {
        const code: u8 = @truncate(canon >> @intCast(2 * @as(u32, i)));
        const t = fromCanonical(code) orelse return null;
        out |= @as(u64, packT1(t)) << @intCast(2 * @as(u32, i));
    }
    return out;
}

/// Parallel trit word: up to 32 trits in a u64 carrier (2 bits each).
pub const TritWord = struct {
    n: u8,
    pack: u64,

    pub fn fromTrits(trits: []const Trit) TritWord {
        var p: u64 = 0;
        const n: u8 = @intCast(@min(trits.len, 32));
        var i: u8 = 0;
        while (i < n) : (i += 1) {
            p |= @as(u64, packT1(trits[i])) << @intCast(2 * i);
        }
        return .{ .n = n, .pack = p };
    }

    pub fn get(self: TritWord, i: u8) ?Trit {
        if (i >= self.n) return null;
        const bits: u8 = @truncate(self.pack >> @intCast(2 * i));
        return unpackT1(bits);
    }
};

/// Parallel pairwise multiply across two words (min length), saturating field.
pub fn pairWords(a: TritWord, b: TritWord) TritWord {
    const n: u8 = @min(a.n, b.n);
    var out: [32]Trit = undefined;
    var i: u8 = 0;
    while (i < n) : (i += 1) {
        const ta = a.get(i) orelse 0;
        const tb = b.get(i) orelse 0;
        out[i] = pair(ta, tb);
    }
    return TritWord.fromTrits(out[0..n]);
}

pub const SelfTestResult = struct {
    ok: bool,
    fails: u32,
};

/// Host- and freestanding-safe self test (no allocator).
pub fn selfTest() SelfTestResult {
    var fails: u32 = 0;

    if (pair(1, -1) != -1) fails += 1;
    if (pair(1, 1) != 1) fails += 1;
    if (sumSat(1, 1) != 1) fails += 1;
    if (sumSat(1, -1) != 0) fails += 1;
    if (consensus(1, 1) != 1) fails += 1;
    if (consensus(1, -1) != 0) fails += 1;
    if (neg(1) != -1) fails += 1;

    if (fromS(-0.5, -0.4, 0.4) != -1) fails += 1;
    if (fromS(0.0, -0.4, 0.4) != 0) fails += 1;
    if (fromS(0.9, -0.4, 0.4) != 1) fails += 1;

    const atg = codonPrimary('A', 'T', 'G');
    // A=+1, T=-1, G=+1
    if (atg[0] != 1 or atg[1] != -1 or atg[2] != 1) fails += 1;

    // RNA: U maps like T (cross-language vector shared with codon_core): AUG -> (+1,-1,+1)
    const aug = codonPrimary('A', 'U', 'G');
    if (aug[0] != 1 or aug[1] != -1 or aug[2] != 1) fails += 1;
    const aug_l = codonPrimary('a', 'u', 'g');
    if (aug_l[0] != 1 or aug_l[1] != -1 or aug_l[2] != 1) fails += 1;

    // canonical <-> trit, and T1 word <-> canonical word (TRIT_SPEC 2a / 2c)
    if (toCanonical(-1) != 0 or toCanonical(0) != 1 or toCanonical(1) != 2) fails += 1;
    if (fromCanonical(0) != -1 or fromCanonical(1) != 0 or fromCanonical(2) != 1) fails += 1;
    if (fromCanonical(3) != null) fails += 1;
    {
        const tv = [_]Trit{ -1, 0, 1, 1, -1 };
        const tw = TritWord.fromTrits(&tv);
        // canonical lanes: 00, 01, 10, 10, 00 (lane 0 = LSB)
        const expect_canon: u64 = 0b00_10_10_01_00;
        const c = t1WordToCanonical(tw.pack, tw.n) orelse 0xFFFF;
        if (c != expect_canon) fails += 1;
        const back = canonicalWordToT1(expect_canon, tw.n) orelse 0xFFFF;
        if (back != tw.pack) fails += 1;
        if (canonicalWordToT1(0b11, 1) != null) fails += 1;
        if (t1WordToCanonical(0b10, 1) != null) fails += 1;
    }

    // packing round-trip
    const ts = [_]Trit{ 1, -1, 0, 1 };
    const w = TritWord.fromTrits(&ts);
    var i: u8 = 0;
    while (i < 4) : (i += 1) {
        const g = w.get(i) orelse {
            fails += 1;
            continue;
        };
        if (g != ts[i]) fails += 1;
    }

    // parallel pair of words
    const wa = TritWord.fromTrits(&[_]Trit{ 1, 1, -1 });
    const wb = TritWord.fromTrits(&[_]Trit{ 1, -1, -1 });
    const wp = pairWords(wa, wb);
    if (wp.get(0) != 1) fails += 1;
    if (wp.get(1) != -1) fails += 1;
    if (wp.get(2) != 1) fails += 1;

    return .{ .ok = fails == 0, .fails = fails };
}

test "trit selfTest (incl. RNA U and canonical layout)" {
    const r = selfTest();
    try @import("std").testing.expect(r.ok);
}

test "AUG primary trits match codon_core" {
    const aug = codonPrimary('A', 'U', 'G');
    try @import("std").testing.expectEqual(@as(Trit, 1), aug[0]);
    try @import("std").testing.expectEqual(@as(Trit, -1), aug[1]);
    try @import("std").testing.expectEqual(@as(Trit, 1), aug[2]);
}
