//! COM1 UART (0x3F8) — freestanding console for QEMU -serial file/stdio.

const COM1: u16 = 0x3F8;

fn outb(port: u16, value: u8) void {
    asm volatile ("outb %[val], %[port]"
        :
        : [val] "{al}" (value),
          [port] "N{dx}" (port),
    );
}

fn inb(port: u16) u8 {
    return asm volatile ("inb %[port], %[ret]"
        : [ret] "={al}" (-> u8),
        : [port] "N{dx}" (port),
    );
}

pub fn init() void {
    outb(COM1 + 1, 0x00);
    outb(COM1 + 3, 0x80);
    outb(COM1 + 0, 0x03);
    outb(COM1 + 1, 0x00);
    outb(COM1 + 3, 0x03);
    outb(COM1 + 2, 0xC7);
    outb(COM1 + 4, 0x0B);
}

fn isTransmitEmpty() bool {
    return (inb(COM1 + 5) & 0x20) != 0;
}

pub fn putc(c: u8) void {
    while (!isTransmitEmpty()) {}
    outb(COM1, c);
}

pub fn write(s: []const u8) void {
    for (s) |c| {
        if (c == '\n') putc('\r');
        putc(c);
    }
}

pub fn writeU32(n: u32) void {
    var buf: [10]u8 = undefined;
    var x = n;
    var i: usize = 0;
    if (x == 0) {
        putc('0');
        return;
    }
    while (x > 0) : (i += 1) {
        buf[i] = @intCast('0' + (x % 10));
        x /= 10;
    }
    while (i > 0) {
        i -= 1;
        putc(buf[i]);
    }
}

/// Non-negative or signed f64 with 3 decimals.
pub fn writeF64_3(x: f64) void {
    if (x != x) {
        write("nan");
        return;
    }
    var v = x;
    if (v < 0) {
        putc('-');
        v = -v;
    }
    const whole: u32 = @intFromFloat(v);
    writeU32(whole);
    putc('.');
    var frac = v - @as(f64, @floatFromInt(whole));
    var d: u32 = 0;
    while (d < 3) : (d += 1) {
        frac *= 10.0;
        const digit: u32 = @intFromFloat(frac);
        putc(@intCast('0' + digit));
        frac -= @as(f64, @floatFromInt(digit));
    }
}
