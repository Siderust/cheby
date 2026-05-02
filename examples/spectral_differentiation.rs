fn main() {
    let d = cheby::spectral::chebyshev_differentiation_matrix(8);
    println!("{}x{}", d.rows(), d.cols());
}
