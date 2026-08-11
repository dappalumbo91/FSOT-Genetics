//! Host product self-test — no QEMU. Run: zig build host
const std = @import("std");
const trit = @import("trit.zig");
const codon = @import("codon.zig");
const scalar = @import("scalar.zig");
const seeds = @import("seeds.zig");
const product = @import("product.zig");

pub fn main() void {
    std.debug.print("FSOT_GENETICS_HOST pin=D1D38A\n", .{});

    // Trit ops
    const t_ok = trit.pair(1, -1) == -1 and trit.sumSat(1, 1) == 1 and trit.consensus(0, 1) == 0;
    std.debug.print("FSOT_TRIT {s}\n", .{if (t_ok) "PASS" else "FAIL"});

    // Codon ATG
    const atg = codon.primaryTrip('A', 'T', 'G');
    const c_ok = atg[0] == 1 and atg[1] == -1 and atg[2] == 1 and codon.dnaToAa('A', 'T', 'G') == 'M';
    std.debug.print("FSOT_CODON {s} ATG=[+1,-1,+1] AA=M\n", .{if (c_ok) "PASS" else "FAIL"});

    // Scalar finite
    const s = scalar.computeNeuro(0.7, 1.0, 1.0);
    const s_ok = s == s and s > -3.1 and s < 3.1;
    std.debug.print("FSOT_SCALAR {s} S={d:.6}\n", .{ if (s_ok) "PASS" else "FAIL", s });

    // Residual product channels
    const r_ok = product.residualSelfTest();
    const r = product.productResiduals();
    std.debug.print("FSOT_RESIDUAL {s}\n", .{if (r_ok) "PASS" else "FAIL"});
    std.debug.print("  r_bond={d:.6} r_clash={d:.6} r_anchor={d:.6}\n", .{ r.r_bond, r.r_clash, r.r_anchor });
    std.debug.print("  multi_top_k={d} multi_power={d:.4} P_NEW={d:.6}\n", .{ product.multi_top_k, product.multi_power, seeds.p_new });

    const cell = product.productCellSelfTest();
    std.debug.print("FSOT_PRODUCT_CELL {s} aa_len={d}\n", .{ if (cell.ok) "PASS" else "FAIL", cell.aa_len });

    const all = t_ok and c_ok and s_ok and r_ok and cell.ok;
    std.debug.print("FSOT_STAGE_GENETICS_{s}\n", .{if (all) "OK" else "FAIL"});
    if (!all) std.process.exit(1);
}
