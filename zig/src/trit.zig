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

/// A,G → +1; C,T → −1
pub fn basePrimary(b: u8) Trit {
    return switch (b) {
        'A', 'a', 'G', 'g' => 1,
        'C', 'c', 'T', 't' => -1,
        else => 0,
    };
}

pub fn codonPrimary(c0: u8, c1: u8, c2: u8) [3]Trit {
    return .{ basePrimary(c0), basePrimary(c1), basePrimary(c2) };
}

// --- T1 packing: 2 bits/trit  00=0, 01=+1, 11=-1 ---

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
