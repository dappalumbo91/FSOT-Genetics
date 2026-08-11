const std = @import("std");

pub fn build(b: *std.Build) void {
    const target_host = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // --- Host product self-test (fast) ---
    const host = b.addExecutable(.{
        .name = "fsot_genetics_host",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main_host.zig"),
            .target = target_host,
            .optimize = optimize,
        }),
    });
    const run_host = b.addRunArtifact(host);
    const host_step = b.step("host", "Build+run genetics product residual self-test on host");
    host_step.dependOn(&b.addInstallArtifact(host, .{}).step);
    host_step.dependOn(&run_host.step);

    // --- Freestanding Multiboot kernel (QEMU -kernel, i386) ---
    // Same path as fsot-neuron-zig: 32-bit freestanding Multiboot1.
    const kernel_target = b.resolveTargetQuery(.{
        .cpu_arch = .x86,
        .os_tag = .freestanding,
        .abi = .none,
    });

    const kernel = b.addExecutable(.{
        .name = "fsot_genetics_kernel",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main_kernel.zig"),
            .target = kernel_target,
            .optimize = .ReleaseSafe,
            .code_model = .kernel,
            .red_zone = false,
        }),
    });
    kernel.entry = .{ .symbol_name = "_start" };
    kernel.setLinkerScript(b.path("linker.ld"));
    kernel.pie = false;
    kernel.link_eh_frame_hdr = false;
    b.installArtifact(kernel);

    const kernel_step = b.step("kernel", "Build freestanding QEMU genetics kernel");
    kernel_step.dependOn(&b.addInstallArtifact(kernel, .{}).step);
}
