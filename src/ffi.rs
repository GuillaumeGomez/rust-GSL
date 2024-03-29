//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#[allow(clippy::upper_case_acronyms)]
pub trait FFI<T> {
    fn wrap(r: *mut T) -> Self;
    fn soft_wrap(r: *mut T) -> Self;
    fn unwrap_shared(&self) -> *const T;
    fn unwrap_unique(&mut self) -> *mut T;
}
