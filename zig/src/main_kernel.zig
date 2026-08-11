//! FSOT Genetics freestanding kernel — Multiboot1 + COM1 serial.
//! Pattern: fsot-neuron-zig main_kernel (QEMU -kernel).
//! Product cell: trit + codon + scalar residual law (no Python).

const trit = @import("trit.zig");
const codon = @import("codon.zig");
const scalar = @import("scalar.zig");
const product = @import("product.zig");
const serial = @import("serial.zig");

const MULTIBOOT_MAGIC: u32 = 0x1BADB002;
const MULTIBOOT_FLAGS: u32 = 0x00000003;
const MULTIBOOT_CHECKSUM: u32 = 0 -% (MULTIBOOT_MAGIC +% MULTIBOOT_FLAGS);

export const multiboot_header align(4) linksection(".multiboot") = [_]u32{
    MULTIBOOT_MAGIC,
    MULTIBOOT_FLAGS,
    MULTIBOOT_CHECKSUM,
};

var stack_bytes: [256 * 1024]u8 align(16) = undefined;

export fn _start() callconv(.c) noreturn {
    const stack_top = @intFromPtr(&stack_bytes) + stack_bytes.len;
    asm volatile (
        \\mov %[sp], %%esp
        \\mov %[sp], %%ebp
        :
        : [sp] "r" (stack_top),
        : .{ .memory = true }
    );
    kmain();
}

fn enableFpu() void {
    asm volatile (
        \\mov %%cr0, %%eax
        \\and $0xFFFFFFF3, %%eax
        \\or  $0x2, %%eax
        \\mov %%eax, %%cr0
        \\mov %%cr4, %%eax
        \\or  $0x600, %%eax
        \\mov %%eax, %%cr4
        \\fninit
        ::: .{ .eax = true, .memory = true }
    );
}

fn kmain() noreturn {
    serial.init();
    serial.write("FSOT_GENETICS_KERNEL pin=D1D38A\n");
    serial.write("product residual cell (bare metal)\n");
    enableFpu();
    serial.write("FPU enabled\n");

    // Trit
    const t_ok = trit.pair(1, -1) == -1 and trit.consensus(1, 1) == 1;
    if (t_ok) serial.write("FSOT_TRIT PASS\n") else serial.write("FSOT_TRIT FAIL\n");

    // Codon ATG
    const atg = codon.primaryTrip('A', 'T', 'G');
    const c_ok = (atg[0] == 1 and atg[1] == -1 and atg[2] == 1) and (codon.dnaToAa('A', 'T', 'G') == 'M');
    if (c_ok) serial.write("FSOT_CODON PASS ATG=[+1,-1,+1] AA=M\n") else serial.write("FSOT_CODON FAIL\n");

    // Scalar
    const s = scalar.computeNeuro(0.7, 1.0, 1.0);
    const s_ok = s == s and s > -3.1 and s < 3.1;
    if (s_ok) {
        serial.write("FSOT_SCALAR PASS S=");
        serial.writeF64_3(s);
        serial.write("\n");
    } else {
        serial.write("FSOT_SCALAR FAIL\n");
    }

    // Residual product channels
    const r_ok = product.residualSelfTest();
    const r = product.productResiduals();
    if (r_ok) {
        serial.write("FSOT_RESIDUAL PASS\n");
        serial.write("  r_bond=");
        serial.writeF64_3(r.r_bond);
        serial.write(" r_clash=");
        serial.writeF64_3(r.r_clash);
        serial.write(" r_anchor=");
        serial.writeF64_3(r.r_anchor);
        serial.write("\n");
    } else {
        serial.write("FSOT_RESIDUAL FAIL\n");
    }

    const cell = product.productCellSelfTest();
    if (cell.ok) {
        serial.write("FSOT_PRODUCT_CELL PASS aa_len=");
        serial.writeU32(cell.aa_len);
        serial.write("\n");
    } else {
        serial.write("FSOT_PRODUCT_CELL FAIL\n");
    }

    const all = t_ok and c_ok and s_ok and r_ok and cell.ok;
    if (all) {
        serial.write("FSOT_STAGE_GENETICS_OK\n");
    } else {
        serial.write("FSOT_STAGE_GENETICS_FAIL\n");
    }

    // Halt
    while (true) {
        asm volatile ("hlt");
    }
}
