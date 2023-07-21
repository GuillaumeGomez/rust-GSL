//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::{Value, VectorF64};
use ffi::FFI;

/// This function constructs a Gaussian kernel parameterized by `alpha` and stores the output in
/// `kernel`. The parameter `order` specifies the derivative order, with `0` corresponding to a
/// Gaussian, `1` corresponding to a first derivative Gaussian, and so on. If `normalize` is set to
/// `true`, then the kernel will be normalized to sum to one on output. If `normalize` is set to
/// `false`, no normalization is performed.
#[doc(alias = "gsl_filter_gaussian_kernel")]
pub fn gaussian_kernel(
    alpha: f64,
    order: usize,
    normalize: bool,
    kernel: &mut VectorF64,
) -> Result<(), Value> {
    let ret = unsafe {
        sys::gsl_filter_gaussian_kernel(alpha, order, normalize as _, kernel.unwrap_unique())
    };
    result_handler!(ret, ())
}
