//! Host product self-test — no QEMU. Run: zig build host
//! Machine-readable PARITY lines for scripts/parity_zig_python.py
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
    std.debug.print("PARITY trit_pair={d}\n", .{@as(i32, trit.pair(1, -1))});
    std.debug.print("PARITY trit_consensus={d}\n", .{@as(i32, trit.consensus(1, 1))});

    // Codon ATG
    const atg = codon.primaryTrip('A', 'T', 'G');
    const c_ok = atg[0] == 1 and atg[1] == -1 and atg[2] == 1 and codon.dnaToAa('A', 'T', 'G') == 'M';
    std.debug.print("FSOT_CODON {s} ATG=[+1,-1,+1] AA=M\n", .{if (c_ok) "PASS" else "FAIL"});
    std.debug.print("PARITY codon_atg0={d} codon_atg1={d} codon_atg2={d} aa_atg={c}\n", .{
        @as(i32, atg[0]),
        @as(i32, atg[1]),
        @as(i32, atg[2]),
        codon.dnaToAa('A', 'T', 'G'),
    });

    // DNA fragment → AA (same string as product.zig productCellSelfTest)
    const dna = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCC";
    var aa_buf: [64]u8 = undefined;
    var aa_n: usize = 0;
    var di: usize = 0;
    while (di + 2 < dna.len and aa_n < aa_buf.len) : (di += 3) {
        const aa = codon.dnaToAa(dna[di], dna[di + 1], dna[di + 2]);
        if (aa == '*') break;
        aa_buf[aa_n] = aa;
        aa_n += 1;
    }
    std.debug.print("PARITY aa_translate={s}\n", .{aa_buf[0..aa_n]});

    // Scalar finite
    const s = scalar.computeNeuro(0.7, 1.0, 1.0);
    const s_ok = s == s and s > -3.1 and s < 3.1;
    std.debug.print("FSOT_SCALAR {s} S={d:.12}\n", .{ if (s_ok) "PASS" else "FAIL", s });
    std.debug.print("PARITY neuro_S={d:.12}\n", .{s});

    // Residual product channels
    const r_ok = product.residualSelfTest();
    const r = product.productResiduals();
    std.debug.print("FSOT_RESIDUAL {s}\n", .{if (r_ok) "PASS" else "FAIL"});
    std.debug.print("  r_bond={d:.12} r_clash={d:.12} r_anchor={d:.12}\n", .{ r.r_bond, r.r_clash, r.r_anchor });
    std.debug.print("  multi_top_k={d} multi_power={d:.12} P_NEW={d:.12}\n", .{ product.multi_top_k, product.multi_power, seeds.p_new });
    std.debug.print("PARITY P_NEW={d:.12}\n", .{seeds.p_new});
    std.debug.print("PARITY K={d:.12}\n", .{seeds.k});
    std.debug.print("PARITY PHI={d:.12}\n", .{seeds.phi});
    std.debug.print("PARITY S_Physical_Chemistry={d:.12}\n", .{r.s_physchem});
    std.debug.print("PARITY S_Chemistry={d:.12}\n", .{r.s_chem});
    std.debug.print("PARITY S_Biochemistry={d:.12}\n", .{r.s_biochem});
    std.debug.print("PARITY r_bond={d:.12}\n", .{r.r_bond});
    std.debug.print("PARITY r_clash={d:.12}\n", .{r.r_clash});
    std.debug.print("PARITY r_anchor={d:.12}\n", .{r.r_anchor});
    std.debug.print("PARITY multi_top_k={d}\n", .{product.multi_top_k});
    std.debug.print("PARITY multi_power={d:.12}\n", .{product.multi_power});

    // Residual physics one-step on a 3-CA stick (bond stretch) — geometry twin check
    const phys = product.residualPhysicsParitySample();
    std.debug.print("PARITY phys_bond_len={d:.12}\n", .{phys.bond_len});
    std.debug.print("PARITY phys_end_x={d:.12}\n", .{phys.end_x});

    const cell = product.productCellSelfTest();
    std.debug.print("FSOT_PRODUCT_CELL {s} aa_len={d}\n", .{ if (cell.ok) "PASS" else "FAIL", cell.aa_len });
    std.debug.print("PARITY aa_len={d}\n", .{cell.aa_len});

    const all = t_ok and c_ok and s_ok and r_ok and cell.ok;
    std.debug.print("FSOT_STAGE_GENETICS_{s}\n", .{if (all) "OK" else "FAIL"});
    if (!all) std.process.exit(1);
}
