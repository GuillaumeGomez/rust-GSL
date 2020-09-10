//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub trait FFI<T> {
    fn wrap(r: *mut T) -> Self;
    fn soft_wrap(r: *mut T) -> Self;
    fn unwrap_shared(&Self) -> *const T;
    fn unwrap_unique(&mut Self) -> *mut T;
}
