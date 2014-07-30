/*
 * A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
 */

pub mod Cblas {
    pub mod Level1 {
        pub fn sdsdot(N: i32, alpha: f32, x: f32, incx: i32, y: f32, incy: i32) -> f32 {
            unsafe { ::ffi::cblas_sdsdot(N, alpha, &x, incx, &y, incy) }
        }
    }
}