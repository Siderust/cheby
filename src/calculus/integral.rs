//! Integral helpers.

use crate::core::{ChebyScalar, ChebySeriesOn, ChebyTime};

pub use crate::core::IntegrateWith;

impl<T, X, const N: usize> ChebySeriesOn<T, X, N>
where
    T: ChebyScalar + IntegrateWith<X>,
    X: ChebyTime,
{
    /// Evaluate the definite integral from the domain start to `x`.
    pub fn evaluate_integral_from_start(
        &self,
        x: X,
    ) -> Result<<T as IntegrateWith<X>>::Integral, crate::core::ChebyError> {
        let antiderivative = self.series().integral(T::zero());
        let tau = self.domain().normalize(x);
        let normalized = antiderivative.evaluate(tau) - antiderivative.evaluate(-1.0);
        Ok(normalized.scale_integral(self.domain().half_width()))
    }
}
