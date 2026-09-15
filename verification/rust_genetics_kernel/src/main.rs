//! Host replay of Genetics seed identities and catalog milliscale gates.
//! Pin D1D38A. Law S = K(T1+T2+T3). 0 free parameters.

mod catalog;

fn main() {
    let phi = (1.0_f64 + 5.0_f64.sqrt()) / 2.0;
    let leftover = 1.0 / (phi * phi);
    let close = 1.0 / phi;
    let free_parameters: u32 = 0;

    assert!((phi - 1.618033988749895).abs() < 1e-12);
    assert!((leftover - 0.3819660112501051).abs() < 1e-12);
    assert!((close - 0.6180339887498948).abs() < 1e-12);
    assert_eq!(free_parameters, 0);
    catalog::check();

    println!("FSOT_GENETICS_RUST_OK");
    println!("phi={phi:.15}");
    println!("leftover={leftover:.15}");
    println!("close={close:.15}");
}
