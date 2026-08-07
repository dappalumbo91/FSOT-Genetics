//! FSOT seed constants as fixed-point (same decimal authority as seeds.zig).
//! No free parameters; SCALE lattice only.

const fixed = @import("fixed.zig");
const Fixed = fixed.Fixed;

pub const pi: Fixed = fixed.fromDecimalStr("3.141592653589");
pub const e: Fixed = fixed.fromDecimalStr("2.718281828459");
pub const phi: Fixed = fixed.fromDecimalStr("1.618033988749");
pub const gamma: Fixed = fixed.fromDecimalStr("0.577215664901");
pub const g_catalan: Fixed = fixed.fromDecimalStr("0.915965594177");

pub const alpha: Fixed = fixed.fromDecimalStr("0.000808293741");
pub const psi_con: Fixed = fixed.fromDecimalStr("0.632120558828");
pub const eta_eff: Fixed = fixed.fromDecimalStr("0.466942206924");
pub const beta: Fixed = fixed.fromDecimalStr("0.00000000000000002620866911333223"); // ~2.62e-17 ≈ 0 at 1e12
pub const chaos: Fixed = fixed.fromDecimalStr("-0.331024182610");
pub const theta_s: Fixed = fixed.fromDecimalStr("0.290896540545");
pub const poof: Fixed = fixed.fromDecimalStr("0.153482214894");

pub const c_eff: Fixed = fixed.fromDecimalStr("0.957702202620");
pub const p_var: Fixed = fixed.fromDecimalStr("0.957987122672");
pub const b_in: Fixed = fixed.fromDecimalStr("0.787940792276");
pub const a_in: Fixed = fixed.fromDecimalStr("1.666853845004");
pub const a_bleed: Fixed = fixed.fromDecimalStr("1.046973630587");
pub const suction: Fixed = fixed.fromDecimalStr("0.147033985428");
pub const p_new: Fixed = fixed.fromDecimalStr("0.300302276670");
pub const c_factor: Fixed = fixed.fromDecimalStr("0.287600151819");
pub const k: Fixed = fixed.fromDecimalStr("0.420221664160");

pub const neuro_d_eff: Fixed = fixed.fromInt(13);
pub const neuro_n_channels: Fixed = fixed.fromInt(4);
pub const neuro_p: Fixed = fixed.fromInt(3);
pub const resting_s: Fixed = fixed.fromDecimalStr("0.46");
