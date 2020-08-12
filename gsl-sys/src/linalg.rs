//! Matrix operations and linear algebra.

use libc::{c_double, c_float, c_int, size_t};

use super::{gsl_complex, gsl_complex_float, gsl_permutation};

extern "C" {
    // Vector functions
    pub fn gsl_vector_alloc(size: size_t) -> *mut gsl_vector;
    pub fn gsl_vector_calloc(size: size_t) -> *mut gsl_vector;
    pub fn gsl_vector_free(vector: *mut gsl_vector);
    pub fn gsl_vector_get(vector: *const gsl_vector, i: size_t) -> c_double;
    pub fn gsl_vector_set(vector: *mut gsl_vector, i: size_t, x: c_double);
    pub fn gsl_vector_set_all(vector: *mut gsl_vector, x: c_double);
    pub fn gsl_vector_set_zero(vector: *mut gsl_vector);
    pub fn gsl_vector_set_basis(vector: *mut gsl_vector, i: size_t);
    pub fn gsl_vector_memcpy(dest: *mut gsl_vector, src: *const gsl_vector) -> c_int;
    pub fn gsl_vector_swap(v: *mut gsl_vector, w: *mut gsl_vector) -> c_int;
    pub fn gsl_vector_swap_elements(vector: *mut gsl_vector, i: size_t, j: size_t) -> c_int;
    pub fn gsl_vector_reverse(vector: *mut gsl_vector) -> c_int;
    pub fn gsl_vector_add(dest: *mut gsl_vector, src: *const gsl_vector) -> c_int;
    pub fn gsl_vector_sub(dest: *mut gsl_vector, src: *const gsl_vector) -> c_int;
    pub fn gsl_vector_mul(dest: *mut gsl_vector, src: *const gsl_vector) -> c_int;
    pub fn gsl_vector_div(dest: *mut gsl_vector, src: *const gsl_vector) -> c_int;
    pub fn gsl_vector_scale(dest: *mut gsl_vector, x: c_double) -> c_int;
    pub fn gsl_vector_add_constant(dest: *mut gsl_vector, x: c_double) -> c_int;
    pub fn gsl_vector_max(vector: *const gsl_vector) -> c_double;
    pub fn gsl_vector_min(vector: *const gsl_vector) -> c_double;
    pub fn gsl_vector_minmax(
        vector: *const gsl_vector,
        min_out: *mut c_double,
        max_out: *mut c_double,
    );
    pub fn gsl_vector_max_index(vector: *const gsl_vector) -> size_t;
    pub fn gsl_vector_min_index(vector: *const gsl_vector) -> size_t;
    pub fn gsl_vector_minmax_index(vector: *const gsl_vector, imin: *mut size_t, imax: *mut size_t);
    pub fn gsl_vector_isnull(vector: *const gsl_vector) -> c_int;
    pub fn gsl_vector_ispos(vector: *const gsl_vector) -> c_int;
    pub fn gsl_vector_isneg(vector: *const gsl_vector) -> c_int;
    pub fn gsl_vector_isnonneg(vector: *const gsl_vector) -> c_int;
    pub fn gsl_vector_equal(u: *const gsl_vector, v: *const gsl_vector) -> c_int;

    // Vector views
    pub fn gsl_vector_subvector(v: *mut gsl_vector, offset: size_t, n: size_t) -> gsl_vector_view;
    pub fn gsl_vector_subvector_with_stride(
        v: *mut gsl_vector,
        offset: size_t,
        stride: size_t,
        n: size_t,
    ) -> gsl_vector_view;
    //pub fn gsl_vector_complex_real(v: *mut gsl_vector_complex) -> gsl_vector_view;
    //pub fn gsl_vector_complex_imag(v: *mut gsl_vector_complex) -> gsl_vector_view;
    pub fn gsl_vector_view_array(base: *mut c_double, n: size_t) -> gsl_vector_view;
    pub fn gsl_vector_view_array_with_stride(
        base: *mut c_double,
        stride: size_t,
        n: size_t,
    ) -> gsl_vector_view;

    // VectorComplex functions
    pub fn gsl_vector_complex_alloc(size: size_t) -> *mut gsl_vector_complex;
    pub fn gsl_vector_complex_calloc(size: size_t) -> *mut gsl_vector_complex;
    pub fn gsl_vector_complex_free(vector: *mut gsl_vector_complex);
    pub fn gsl_vector_complex_get(vector: *const gsl_vector_complex, i: size_t) -> gsl_complex;
    pub fn gsl_vector_complex_set(vector: *mut gsl_vector_complex, i: size_t, x: gsl_complex);
    pub fn gsl_vector_complex_set_all(vector: *mut gsl_vector_complex, x: gsl_complex);
    pub fn gsl_vector_complex_set_zero(vector: *mut gsl_vector_complex);
    pub fn gsl_vector_complex_set_basis(vector: *mut gsl_vector_complex, i: size_t);
    pub fn gsl_vector_complex_memcpy(
        dest: *mut gsl_vector_complex,
        src: *const gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_vector_complex_swap(v: *mut gsl_vector_complex, w: *mut gsl_vector_complex)
        -> c_int;
    pub fn gsl_vector_complex_swap_elements(
        vector: *mut gsl_vector_complex,
        i: size_t,
        j: size_t,
    ) -> c_int;
    pub fn gsl_vector_complex_reverse(vector: *mut gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_add(
        dest: *mut gsl_vector_complex,
        src: *const gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_vector_complex_sub(
        dest: *mut gsl_vector_complex,
        src: *const gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_vector_complex_mul(
        dest: *mut gsl_vector_complex,
        src: *const gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_vector_complex_div(
        dest: *mut gsl_vector_complex,
        src: *const gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_vector_complex_scale(dest: *mut gsl_vector_complex, x: gsl_complex) -> c_int;
    pub fn gsl_vector_complex_add_constant(dest: *mut gsl_vector_complex, x: gsl_complex) -> c_int;
    pub fn gsl_vector_complex_isnull(vector: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_ispos(vector: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_isneg(vector: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_isnonneg(vector: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_equal(
        u: *const gsl_vector_complex,
        v: *const gsl_vector_complex,
    ) -> c_int;

    // VectorFloat functions
    pub fn gsl_vector_float_alloc(size: size_t) -> *mut gsl_vector_float;
    pub fn gsl_vector_float_calloc(size: size_t) -> *mut gsl_vector_float;
    pub fn gsl_vector_float_free(vector: *mut gsl_vector_float);
    pub fn gsl_vector_float_get(vector: *const gsl_vector_float, i: size_t) -> c_float;
    pub fn gsl_vector_float_set(vector: *mut gsl_vector_float, i: size_t, x: c_float);
    pub fn gsl_vector_float_set_all(vector: *mut gsl_vector_float, x: c_float);
    pub fn gsl_vector_float_set_zero(vector: *mut gsl_vector_float);
    pub fn gsl_vector_float_set_basis(vector: *mut gsl_vector_float, i: size_t);
    pub fn gsl_vector_float_memcpy(
        dest: *mut gsl_vector_float,
        src: *const gsl_vector_float,
    ) -> c_int;
    pub fn gsl_vector_float_swap(v: *mut gsl_vector_float, w: *mut gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_swap_elements(
        vector: *mut gsl_vector_float,
        i: size_t,
        j: size_t,
    ) -> c_int;
    pub fn gsl_vector_float_reverse(vector: *mut gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_add(dest: *mut gsl_vector_float, src: *const gsl_vector_float)
        -> c_int;
    pub fn gsl_vector_float_sub(dest: *mut gsl_vector_float, src: *const gsl_vector_float)
        -> c_int;
    pub fn gsl_vector_float_mul(dest: *mut gsl_vector_float, src: *const gsl_vector_float)
        -> c_int;
    pub fn gsl_vector_float_div(dest: *mut gsl_vector_float, src: *const gsl_vector_float)
        -> c_int;
    pub fn gsl_vector_float_scale(dest: *mut gsl_vector_float, x: c_float) -> c_int;
    pub fn gsl_vector_float_add_constant(dest: *mut gsl_vector_float, x: c_float) -> c_int;
    pub fn gsl_vector_float_max(vector: *const gsl_vector_float) -> c_float;
    pub fn gsl_vector_float_min(vector: *const gsl_vector_float) -> c_float;
    pub fn gsl_vector_float_minmax(
        vector: *const gsl_vector_float,
        min_out: *mut c_float,
        max_out: *mut c_float,
    );
    pub fn gsl_vector_float_max_index(vector: *const gsl_vector_float) -> size_t;
    pub fn gsl_vector_float_min_index(vector: *const gsl_vector_float) -> size_t;
    pub fn gsl_vector_float_minmax_index(
        vector: *const gsl_vector_float,
        imin: *mut size_t,
        imax: *mut size_t,
    );
    pub fn gsl_vector_float_isnull(vector: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_ispos(vector: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_isneg(vector: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_isnonneg(vector: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_equal(u: *const gsl_vector_float, v: *const gsl_vector_float) -> c_int;

    // VectorComplexFloat functions
    pub fn gsl_vector_complex_float_alloc(size: size_t) -> *mut gsl_vector_complex_float;
    pub fn gsl_vector_complex_float_calloc(size: size_t) -> *mut gsl_vector_complex_float;
    pub fn gsl_vector_complex_float_free(vector: *mut gsl_vector_complex_float);
    pub fn gsl_vector_complex_float_get(
        vector: *const gsl_vector_complex_float,
        i: size_t,
    ) -> gsl_complex_float;
    pub fn gsl_vector_complex_float_set(
        vector: *mut gsl_vector_complex_float,
        i: size_t,
        x: gsl_complex_float,
    );
    pub fn gsl_vector_complex_float_set_all(
        vector: *mut gsl_vector_complex_float,
        x: gsl_complex_float,
    );
    pub fn gsl_vector_complex_float_set_zero(vector: *mut gsl_vector_complex_float);
    pub fn gsl_vector_complex_float_set_basis(vector: *mut gsl_vector_complex_float, i: size_t);
    pub fn gsl_vector_complex_float_memcpy(
        dest: *mut gsl_vector_complex_float,
        src: *const gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_vector_complex_float_swap(
        v: *mut gsl_vector_complex_float,
        w: *mut gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_vector_complex_float_swap_elements(
        vector: *mut gsl_vector_complex_float,
        i: size_t,
        j: size_t,
    ) -> c_int;
    pub fn gsl_vector_complex_float_reverse(vector: *mut gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_add(
        dest: *mut gsl_vector_complex_float,
        src: *const gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_vector_complex_float_sub(
        dest: *mut gsl_vector_complex_float,
        src: *const gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_vector_complex_float_mul(
        dest: *mut gsl_vector_complex_float,
        src: *const gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_vector_complex_float_div(
        dest: *mut gsl_vector_complex_float,
        src: *const gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_vector_complex_float_scale(
        dest: *mut gsl_vector_complex_float,
        x: gsl_complex_float,
    ) -> c_int;
    pub fn gsl_vector_complex_float_add_constant(
        dest: *mut gsl_vector_complex_float,
        x: gsl_complex_float,
    ) -> c_int;
    pub fn gsl_vector_complex_float_isnull(vector: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_ispos(vector: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_isneg(vector: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_isnonneg(vector: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_equal(
        u: *const gsl_vector_complex_float,
        v: *const gsl_vector_complex_float,
    ) -> c_int;

    // Matrix functions
    pub fn gsl_matrix_alloc(size1: size_t, size2: size_t) -> *mut gsl_matrix;
    pub fn gsl_matrix_calloc(size1: size_t, size2: size_t) -> *mut gsl_matrix;
    pub fn gsl_matrix_free(m: *mut gsl_matrix);
    pub fn gsl_matrix_get(m: *const gsl_matrix, i: size_t, j: size_t) -> c_double;
    pub fn gsl_matrix_set(m: *mut gsl_matrix, i: size_t, j: size_t, x: c_double);
    pub fn gsl_matrix_set_all(m: *mut gsl_matrix, x: c_double);
    pub fn gsl_matrix_set_zero(m: *mut gsl_matrix);
    pub fn gsl_matrix_set_identity(m: *mut gsl_matrix);
    pub fn gsl_matrix_memcpy(dest: *mut gsl_matrix, src: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_swap(m: *mut gsl_matrix, w: *mut gsl_matrix) -> c_int;
    pub fn gsl_matrix_get_row(vector: *mut gsl_vector, m: *const gsl_matrix, i: size_t) -> c_int;
    pub fn gsl_matrix_get_col(vector: *mut gsl_vector, m: *const gsl_matrix, j: size_t) -> c_int;
    pub fn gsl_matrix_set_row(m: *mut gsl_matrix, i: size_t, v: *const gsl_vector) -> c_int;
    pub fn gsl_matrix_set_col(m: *mut gsl_matrix, j: size_t, v: *const gsl_vector) -> c_int;
    pub fn gsl_matrix_swap_rows(m: *mut gsl_matrix, i: size_t, j: size_t) -> c_int;
    pub fn gsl_matrix_swap_columns(m: *mut gsl_matrix, i: size_t, j: size_t) -> c_int;
    pub fn gsl_matrix_swap_rowcol(m: *mut gsl_matrix, i: size_t, j: size_t) -> c_int;
    pub fn gsl_matrix_CBLAS_TRANSPOSE_t_memcpy(
        dest: *mut gsl_matrix,
        src: *const gsl_matrix,
    ) -> c_int;
    pub fn gsl_matrix_CBLAS_TRANSPOSE_t(m: *mut gsl_matrix) -> c_int;
    pub fn gsl_matrix_add(dest: *mut gsl_matrix, src: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_sub(dest: *mut gsl_matrix, src: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_mul_elements(dest: *mut gsl_matrix, src: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_div_elements(dest: *mut gsl_matrix, src: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_scale(dest: *mut gsl_matrix, x: c_double) -> c_int;
    pub fn gsl_matrix_add_constant(dest: *mut gsl_matrix, x: c_double) -> c_int;
    pub fn gsl_matrix_max(m: *const gsl_matrix) -> c_double;
    pub fn gsl_matrix_min(m: *const gsl_matrix) -> c_double;
    pub fn gsl_matrix_minmax(m: *const gsl_matrix, min_out: *mut c_double, max_out: *mut c_double);
    pub fn gsl_matrix_max_index(m: *const gsl_matrix, imax: *mut size_t, jmax: *mut size_t);
    pub fn gsl_matrix_min_index(m: *const gsl_matrix, imin: *mut size_t, jmin: *mut size_t);
    pub fn gsl_matrix_minmax_index(
        m: *const gsl_matrix,
        imin: *mut size_t,
        jmin: *mut size_t,
        imax: *mut size_t,
        jmax: *mut size_t,
    );
    pub fn gsl_matrix_isnull(m: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_ispos(m: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_isneg(m: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_isnonneg(m: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_equal(u: *const gsl_matrix, v: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_transpose_memcpy(dest: *mut gsl_matrix, src: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_transpose(m: *mut gsl_matrix) -> c_int;

    // Matrix views
    pub fn gsl_matrix_submatrix(
        m: *mut gsl_matrix,
        k1: size_t,
        k2: size_t,
        n1: size_t,
        n2: size_t,
    ) -> gsl_matrix_view;
    pub fn gsl_matrix_view_array(base: *mut c_double, n1: size_t, n2: size_t) -> gsl_matrix_view;
    pub fn gsl_matrix_view_array_with_tda(
        base: *mut c_double,
        n1: size_t,
        n2: size_t,
        tda: size_t,
    ) -> gsl_matrix_view;
    pub fn gsl_matrix_view_vector(v: *mut gsl_vector, n1: size_t, n2: size_t) -> gsl_matrix_view;
    pub fn gsl_matrix_view_vector_with_tda(
        v: *mut gsl_vector,
        n1: size_t,
        n2: size_t,
        tda: size_t,
    ) -> gsl_matrix_view;

    // MatrixFloat functions
    pub fn gsl_matrix_float_alloc(size1: size_t, size2: size_t) -> *mut gsl_matrix_float;
    pub fn gsl_matrix_float_calloc(size1: size_t, size2: size_t) -> *mut gsl_matrix_float;
    pub fn gsl_matrix_float_free(m: *mut gsl_matrix_float);
    pub fn gsl_matrix_float_get(m: *const gsl_matrix_float, i: size_t, j: size_t) -> c_float;
    pub fn gsl_matrix_float_set(m: *mut gsl_matrix_float, i: size_t, j: size_t, x: c_float);
    pub fn gsl_matrix_float_set_all(m: *mut gsl_matrix_float, x: c_float);
    pub fn gsl_matrix_float_set_zero(m: *mut gsl_matrix_float);
    pub fn gsl_matrix_float_set_identity(m: *mut gsl_matrix_float);
    pub fn gsl_matrix_float_memcpy(
        dest: *mut gsl_matrix_float,
        src: *const gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_matrix_float_swap(m: *mut gsl_matrix_float, w: *mut gsl_matrix_float) -> c_int;
    pub fn gsl_matrix_float_get_row(
        vector: *mut gsl_vector_float,
        m: *const gsl_matrix_float,
        i: size_t,
    ) -> c_int;
    pub fn gsl_matrix_float_get_col(
        vector: *mut gsl_vector_float,
        m: *const gsl_matrix_float,
        j: size_t,
    ) -> c_int;
    pub fn gsl_matrix_float_set_row(
        m: *mut gsl_matrix_float,
        i: size_t,
        v: *const gsl_vector_float,
    ) -> c_int;
    pub fn gsl_matrix_float_set_col(
        m: *mut gsl_matrix_float,
        j: size_t,
        v: *const gsl_vector_float,
    ) -> c_int;
    pub fn gsl_matrix_float_swap_rows(m: *mut gsl_matrix_float, i: size_t, j: size_t) -> c_int;
    pub fn gsl_matrix_float_swap_columns(m: *mut gsl_matrix_float, i: size_t, j: size_t) -> c_int;
    pub fn gsl_matrix_float_swap_rowcol(m: *mut gsl_matrix_float, i: size_t, j: size_t) -> c_int;
    pub fn gsl_matrix_float_CBLAS_TRANSPOSE_t_memcpy(
        dest: *mut gsl_matrix_float,
        src: *const gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_matrix_float_CBLAS_TRANSPOSE_t(m: *mut gsl_matrix_float) -> c_int;
    pub fn gsl_matrix_float_add(dest: *mut gsl_matrix_float, src: *const gsl_matrix_float)
        -> c_int;
    pub fn gsl_matrix_float_sub(dest: *mut gsl_matrix_float, src: *const gsl_matrix_float)
        -> c_int;
    pub fn gsl_matrix_float_mul_elements(
        dest: *mut gsl_matrix_float,
        src: *const gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_matrix_float_div_elements(
        dest: *mut gsl_matrix_float,
        src: *const gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_matrix_float_scale(dest: *mut gsl_matrix_float, x: c_float) -> c_int;
    pub fn gsl_matrix_float_add_constant(dest: *mut gsl_matrix_float, x: c_float) -> c_int;
    pub fn gsl_matrix_float_max(m: *const gsl_matrix_float) -> c_float;
    pub fn gsl_matrix_float_min(m: *const gsl_matrix_float) -> c_float;
    pub fn gsl_matrix_float_minmax(
        m: *const gsl_matrix_float,
        min_out: *mut c_float,
        max_out: *mut c_float,
    );
    pub fn gsl_matrix_float_max_index(
        m: *const gsl_matrix_float,
        imax: *mut size_t,
        jmax: *mut size_t,
    );
    pub fn gsl_matrix_float_min_index(
        m: *const gsl_matrix_float,
        imin: *mut size_t,
        jmin: *mut size_t,
    );
    pub fn gsl_matrix_float_minmax_index(
        m: *const gsl_matrix_float,
        imin: *mut size_t,
        jmin: *mut size_t,
        imax: *mut size_t,
        jmax: *mut size_t,
    );
    pub fn gsl_matrix_float_isnull(m: *const gsl_matrix_float) -> c_int;
    pub fn gsl_matrix_float_ispos(m: *const gsl_matrix_float) -> c_int;
    pub fn gsl_matrix_float_isneg(m: *const gsl_matrix_float) -> c_int;
    pub fn gsl_matrix_float_isnonneg(m: *const gsl_matrix_float) -> c_int;
    pub fn gsl_matrix_float_equal(u: *const gsl_matrix_float, v: *const gsl_matrix_float) -> c_int;
    pub fn gsl_matrix_float_transpose_memcpy(
        dest: *mut gsl_matrix_float,
        src: *const gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_matrix_float_transpose(m: *mut gsl_matrix_float) -> c_int;

    // MatrixComplex functions
    pub fn gsl_matrix_complex_alloc(size1: size_t, size2: size_t) -> *mut gsl_matrix_complex;
    pub fn gsl_matrix_complex_calloc(size1: size_t, size2: size_t) -> *mut gsl_matrix_complex;
    pub fn gsl_matrix_complex_free(m: *mut gsl_matrix_complex);
    pub fn gsl_matrix_complex_get(
        m: *const gsl_matrix_complex,
        i: size_t,
        j: size_t,
    ) -> gsl_complex;
    pub fn gsl_matrix_complex_set(m: *mut gsl_matrix_complex, i: size_t, j: size_t, x: gsl_complex);
    pub fn gsl_matrix_complex_set_all(m: *mut gsl_matrix_complex, x: gsl_complex);
    pub fn gsl_matrix_complex_set_zero(m: *mut gsl_matrix_complex);
    pub fn gsl_matrix_complex_set_identity(m: *mut gsl_matrix_complex);
    pub fn gsl_matrix_complex_memcpy(
        dest: *mut gsl_matrix_complex,
        src: *const gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_matrix_complex_swap(m: *mut gsl_matrix_complex, w: *mut gsl_matrix_complex)
        -> c_int;
    pub fn gsl_matrix_complex_get_row(
        vector: *mut gsl_vector_complex,
        m: *const gsl_matrix_complex,
        i: size_t,
    ) -> c_int;
    pub fn gsl_matrix_complex_get_col(
        vector: *mut gsl_vector_complex,
        m: *const gsl_matrix_complex,
        j: size_t,
    ) -> c_int;
    pub fn gsl_matrix_complex_set_row(
        m: *mut gsl_matrix_complex,
        i: size_t,
        v: *const gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_matrix_complex_set_col(
        m: *mut gsl_matrix_complex,
        j: size_t,
        v: *const gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_matrix_complex_swap_rows(m: *mut gsl_matrix_complex, i: size_t, j: size_t) -> c_int;
    pub fn gsl_matrix_complex_swap_columns(
        m: *mut gsl_matrix_complex,
        i: size_t,
        j: size_t,
    ) -> c_int;
    pub fn gsl_matrix_complex_swap_rowcol(
        m: *mut gsl_matrix_complex,
        i: size_t,
        j: size_t,
    ) -> c_int;
    pub fn gsl_matrix_complex_CBLAS_TRANSPOSE_t_memcpy(
        dest: *mut gsl_matrix_complex,
        src: *const gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_matrix_complex_CBLAS_TRANSPOSE_t(m: *mut gsl_matrix_complex) -> c_int;
    pub fn gsl_matrix_complex_add(
        dest: *mut gsl_matrix_complex,
        src: *const gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_matrix_complex_sub(
        dest: *mut gsl_matrix_complex,
        src: *const gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_matrix_complex_mul_elements(
        dest: *mut gsl_matrix_complex,
        src: *const gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_matrix_complex_div_elements(
        dest: *mut gsl_matrix_complex,
        src: *const gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_matrix_complex_scale(dest: *mut gsl_matrix_complex, x: gsl_complex) -> c_int;
    pub fn gsl_matrix_complex_add_constant(dest: *mut gsl_matrix_complex, x: gsl_complex) -> c_int;
    pub fn gsl_matrix_complex_isnull(m: *const gsl_matrix_complex) -> c_int;
    pub fn gsl_matrix_complex_ispos(m: *const gsl_matrix_complex) -> c_int;
    pub fn gsl_matrix_complex_isneg(m: *const gsl_matrix_complex) -> c_int;
    pub fn gsl_matrix_complex_isnonneg(m: *const gsl_matrix_complex) -> c_int;
    pub fn gsl_matrix_complex_equal(
        u: *const gsl_matrix_complex,
        v: *const gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_matrix_complex_transpose_memcpy(
        dest: *mut gsl_matrix_complex,
        src: *const gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_matrix_complex_transpose(m: *mut gsl_matrix_complex) -> c_int;

    // MatrixComplexFloat functions
    pub fn gsl_matrix_complex_float_alloc(
        size1: size_t,
        size2: size_t,
    ) -> *mut gsl_matrix_complex_float;
    pub fn gsl_matrix_complex_float_calloc(
        size1: size_t,
        size2: size_t,
    ) -> *mut gsl_matrix_complex_float;
    pub fn gsl_matrix_complex_float_free(m: *mut gsl_matrix_complex_float);
    pub fn gsl_matrix_complex_float_get(
        m: *const gsl_matrix_complex_float,
        i: size_t,
        j: size_t,
    ) -> gsl_complex_float;
    pub fn gsl_matrix_complex_float_set(
        m: *mut gsl_matrix_complex_float,
        i: size_t,
        j: size_t,
        x: gsl_complex_float,
    );
    pub fn gsl_matrix_complex_float_set_all(m: *mut gsl_matrix_complex_float, x: gsl_complex_float);
    pub fn gsl_matrix_complex_float_set_zero(m: *mut gsl_matrix_complex_float);
    pub fn gsl_matrix_complex_float_set_identity(m: *mut gsl_matrix_complex_float);
    pub fn gsl_matrix_complex_float_memcpy(
        dest: *mut gsl_matrix_complex_float,
        src: *const gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_swap(
        m: *mut gsl_matrix_complex_float,
        w: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_get_row(
        vector: *mut gsl_vector_complex_float,
        m: *const gsl_matrix_complex_float,
        i: size_t,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_get_col(
        vector: *mut gsl_vector_complex_float,
        m: *const gsl_matrix_complex_float,
        j: size_t,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_set_row(
        m: *mut gsl_matrix_complex_float,
        i: size_t,
        v: *const gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_set_col(
        m: *mut gsl_matrix_complex_float,
        j: size_t,
        v: *const gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_swap_rows(
        m: *mut gsl_matrix_complex_float,
        i: size_t,
        j: size_t,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_swap_columns(
        m: *mut gsl_matrix_complex_float,
        i: size_t,
        j: size_t,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_swap_rowcol(
        m: *mut gsl_matrix_complex_float,
        i: size_t,
        j: size_t,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_CBLAS_TRANSPOSE_t_memcpy(
        dest: *mut gsl_matrix_complex_float,
        src: *const gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_CBLAS_TRANSPOSE_t(m: *mut gsl_matrix_complex_float) -> c_int;
    pub fn gsl_matrix_complex_float_add(
        dest: *mut gsl_matrix_complex_float,
        src: *const gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_sub(
        dest: *mut gsl_matrix_complex_float,
        src: *const gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_mul_elements(
        dest: *mut gsl_matrix_complex_float,
        src: *const gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_div_elements(
        dest: *mut gsl_matrix_complex_float,
        src: *const gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_scale(
        dest: *mut gsl_matrix_complex_float,
        x: gsl_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_add_constant(
        dest: *mut gsl_matrix_complex_float,
        x: gsl_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_isnull(m: *const gsl_matrix_complex_float) -> c_int;
    pub fn gsl_matrix_complex_float_ispos(m: *const gsl_matrix_complex_float) -> c_int;
    pub fn gsl_matrix_complex_float_isneg(m: *const gsl_matrix_complex_float) -> c_int;
    pub fn gsl_matrix_complex_float_isnonneg(m: *const gsl_matrix_complex_float) -> c_int;
    pub fn gsl_matrix_complex_float_equal(
        u: *const gsl_matrix_complex_float,
        v: *const gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_transpose_memcpy(
        dest: *mut gsl_matrix_complex_float,
        src: *const gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_matrix_complex_float_transpose(m: *mut gsl_matrix_complex_float) -> c_int;

    // Real Symmetric Matrices
    pub fn gsl_eigen_symm_alloc(n: size_t) -> *mut gsl_eigen_symm_workspace;
    pub fn gsl_eigen_symm_free(w: *mut gsl_eigen_symm_workspace);
    pub fn gsl_eigen_symm(
        A: *mut gsl_matrix,
        eval: *mut gsl_vector,
        w: *mut gsl_eigen_symm_workspace,
    ) -> c_int;
    pub fn gsl_eigen_symmv_alloc(n: size_t) -> *mut gsl_eigen_symmv_workspace;
    pub fn gsl_eigen_symmv_free(w: *mut gsl_eigen_symmv_workspace);
    pub fn gsl_eigen_symmv(
        A: *mut gsl_matrix,
        eval: *mut gsl_vector,
        evec: *mut gsl_matrix,
        w: *mut gsl_eigen_symmv_workspace,
    ) -> c_int;
    // Complex Hermitian Matrices
    pub fn gsl_eigen_herm_alloc(n: size_t) -> *mut gsl_eigen_herm_workspace;
    pub fn gsl_eigen_herm_free(w: *mut gsl_eigen_herm_workspace);
    pub fn gsl_eigen_herm(
        A: *mut gsl_matrix_complex,
        eval: *mut gsl_vector,
        w: *mut gsl_eigen_herm_workspace,
    ) -> c_int;
    pub fn gsl_eigen_hermv_alloc(n: size_t) -> *mut gsl_eigen_hermv_workspace;
    pub fn gsl_eigen_hermv_free(w: *mut gsl_eigen_hermv_workspace);
    pub fn gsl_eigen_hermv(
        A: *mut gsl_matrix_complex,
        eval: *mut gsl_vector,
        evec: *mut gsl_matrix_complex,
        w: *mut gsl_eigen_hermv_workspace,
    ) -> c_int;
    // Real Nonsymmetric Matrices
    pub fn gsl_eigen_nonsymm_alloc(n: size_t) -> *mut gsl_eigen_nonsymm_workspace;
    pub fn gsl_eigen_nonsymm_free(w: *mut gsl_eigen_nonsymm_workspace);
    pub fn gsl_eigen_nonsymm_params(
        compute_t: c_int,
        balance: c_int,
        w: *mut gsl_eigen_nonsymm_workspace,
    );
    pub fn gsl_eigen_nonsymm(
        A: *mut gsl_matrix,
        eval: *mut gsl_vector_complex,
        w: *mut gsl_eigen_nonsymm_workspace,
    ) -> c_int;
    pub fn gsl_eigen_nonsymm_Z(
        A: *mut gsl_matrix,
        eval: *mut gsl_vector_complex,
        z: *mut gsl_matrix,
        w: *mut gsl_eigen_nonsymm_workspace,
    ) -> c_int;
    pub fn gsl_eigen_nonsymmv_alloc(n: size_t) -> *mut gsl_eigen_nonsymmv_workspace;
    pub fn gsl_eigen_nonsymmv_free(w: *mut gsl_eigen_nonsymmv_workspace);
    pub fn gsl_eigen_nonsymmv_params(balance: c_int, w: *mut gsl_eigen_nonsymmv_workspace);
    pub fn gsl_eigen_nonsymmv(
        A: *mut gsl_matrix,
        eval: *mut gsl_vector_complex,
        evec: *mut gsl_matrix_complex,
        w: *mut gsl_eigen_nonsymmv_workspace,
    ) -> c_int;
    pub fn gsl_eigen_nonsymmv_Z(
        A: *mut gsl_matrix,
        eval: *mut gsl_vector_complex,
        evec: *mut gsl_matrix_complex,
        z: *mut gsl_matrix,
        w: *mut gsl_eigen_nonsymmv_workspace,
    ) -> c_int;
    // Real Generalized Symmetric-Definite Eigensystems
    pub fn gsl_eigen_gensymm_alloc(n: size_t) -> *mut gsl_eigen_gensymm_workspace;
    pub fn gsl_eigen_gensymm_free(w: *mut gsl_eigen_gensymm_workspace);
    pub fn gsl_eigen_gensymm(
        A: *mut gsl_matrix,
        B: *mut gsl_matrix,
        eval: *mut gsl_vector,
        w: *mut gsl_eigen_gensymm_workspace,
    ) -> c_int;
    pub fn gsl_eigen_gensymmv_alloc(n: size_t) -> *mut gsl_eigen_gensymmv_workspace;
    pub fn gsl_eigen_gensymmv_free(w: *mut gsl_eigen_gensymmv_workspace);
    pub fn gsl_eigen_gensymmv(
        A: *mut gsl_matrix,
        B: *mut gsl_matrix,
        eval: *mut gsl_vector,
        evec: *mut gsl_matrix,
        w: *mut gsl_eigen_gensymmv_workspace,
    ) -> c_int;
    // Complex Generalized Hermitian-Definite Eigensystems
    pub fn gsl_eigen_genherm_alloc(n: size_t) -> *mut gsl_eigen_genherm_workspace;
    pub fn gsl_eigen_genherm_free(w: *mut gsl_eigen_genherm_workspace);
    pub fn gsl_eigen_genherm(
        A: *mut gsl_matrix_complex,
        B: *mut gsl_matrix_complex,
        eval: *mut gsl_vector,
        w: *mut gsl_eigen_genherm_workspace,
    ) -> c_int;
    pub fn gsl_eigen_genhermv_alloc(n: size_t) -> *mut gsl_eigen_genhermv_workspace;
    pub fn gsl_eigen_genhermv_free(w: *mut gsl_eigen_genhermv_workspace);
    pub fn gsl_eigen_genhermv(
        A: *mut gsl_matrix_complex,
        B: *mut gsl_matrix_complex,
        eval: *mut gsl_vector,
        evec: *mut gsl_matrix_complex,
        w: *mut gsl_eigen_genhermv_workspace,
    ) -> c_int;
    // Real Generalized Nonsymmetric Eigensystems
    pub fn gsl_eigen_gen_alloc(n: size_t) -> *mut gsl_eigen_gen_workspace;
    pub fn gsl_eigen_gen_free(w: *mut gsl_eigen_gen_workspace);
    pub fn gsl_eigen_gen_params(
        compute_s: c_int,
        compute_t: c_int,
        balance: c_int,
        w: *mut gsl_eigen_gen_workspace,
    );
    pub fn gsl_eigen_gen(
        A: *mut gsl_matrix,
        B: *mut gsl_matrix,
        alpha: *mut gsl_vector_complex,
        beta: *mut gsl_vector,
        w: *mut gsl_eigen_gen_workspace,
    ) -> c_int;
    pub fn gsl_eigen_gen_QZ(
        A: *mut gsl_matrix,
        B: *mut gsl_matrix,
        alpha: *mut gsl_vector_complex,
        beta: *mut gsl_vector,
        Q: *mut gsl_matrix,
        Z: *mut gsl_matrix,
        w: *mut gsl_eigen_gen_workspace,
    ) -> c_int;
    pub fn gsl_eigen_genv_alloc(n: size_t) -> *mut gsl_eigen_genv_workspace;
    pub fn gsl_eigen_genv_free(w: *mut gsl_eigen_genv_workspace);
    pub fn gsl_eigen_genv(
        A: *mut gsl_matrix,
        B: *mut gsl_matrix,
        alpha: *mut gsl_vector_complex,
        beta: *mut gsl_vector,
        evec: *mut gsl_matrix_complex,
        w: *mut gsl_eigen_genv_workspace,
    ) -> c_int;
    pub fn gsl_eigen_genv_QZ(
        A: *mut gsl_matrix,
        B: *mut gsl_matrix,
        alpha: *mut gsl_vector_complex,
        beta: *mut gsl_vector,
        evec: *mut gsl_matrix_complex,
        Q: *mut gsl_matrix,
        Z: *mut gsl_matrix,
        w: *mut gsl_eigen_genv_workspace,
    ) -> c_int;
    // Sorting Eigenvalues and Eigenvectors
    pub fn gsl_eigen_symmv_sort(
        eval: *mut gsl_vector,
        evec: *mut gsl_matrix,
        sort_type: c_int,
    ) -> c_int;
    pub fn gsl_eigen_hermv_sort(
        eval: *mut gsl_vector,
        evec: *mut gsl_matrix_complex,
        sort_type: c_int,
    ) -> c_int;
    pub fn gsl_eigen_nonsymmv_sort(
        eval: *mut gsl_vector_complex,
        evec: *mut gsl_matrix_complex,
        sort_type: c_int,
    ) -> c_int;
    pub fn gsl_eigen_gensymmv_sort(
        eval: *mut gsl_vector,
        evec: *mut gsl_matrix,
        sort_type: c_int,
    ) -> c_int;
    pub fn gsl_eigen_genhermv_sort(
        eval: *mut gsl_vector,
        evec: *mut gsl_matrix_complex,
        sort_type: c_int,
    ) -> c_int;
    pub fn gsl_eigen_genv_sort(
        alpha: *mut gsl_vector_complex,
        beta: *mut gsl_vector,
        evec: *mut gsl_matrix_complex,
        sort_type: c_int,
    ) -> c_int;

    // LU Decomposition
    pub fn gsl_linalg_LU_decomp(
        a: *mut gsl_matrix,
        p: *mut gsl_permutation,
        signum: *mut c_int,
    ) -> c_int;
    pub fn gsl_linalg_complex_LU_decomp(
        a: *mut gsl_matrix_complex,
        p: *mut gsl_permutation,
        signum: *mut c_int,
    ) -> c_int;
    pub fn gsl_linalg_LU_solve(
        lu: *const gsl_matrix,
        p: *const gsl_permutation,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_complex_LU_solve(
        lu: *const gsl_matrix_complex,
        p: *const gsl_permutation,
        b: *const gsl_vector_complex,
        x: *mut gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_linalg_LU_svx(
        lu: *const gsl_matrix,
        p: *const gsl_permutation,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_complex_LU_svx(
        lu: *const gsl_matrix_complex,
        p: *const gsl_permutation,
        x: *mut gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_linalg_LU_refine(
        a: *const gsl_matrix,
        lu: *const gsl_matrix,
        p: *const gsl_permutation,
        b: *const gsl_vector,
        x: *mut gsl_vector,
        residual: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_complex_LU_refine(
        a: *const gsl_matrix_complex,
        lu: *const gsl_matrix_complex,
        p: *const gsl_permutation,
        b: *const gsl_vector_complex,
        x: *mut gsl_vector_complex,
        residual: *mut gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_linalg_LU_invert(
        lu: *const gsl_matrix,
        p: *const gsl_permutation,
        inverse: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_linalg_complex_LU_invert(
        lu: *const gsl_matrix_complex,
        p: *const gsl_permutation,
        inverse: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_linalg_LU_det(lu: *mut gsl_matrix, signum: c_int) -> c_double;
    pub fn gsl_linalg_complex_LU_det(lu: *mut gsl_matrix_complex, signum: c_int) -> gsl_complex;
    pub fn gsl_linalg_LU_lndet(lu: *mut gsl_matrix) -> c_double;
    pub fn gsl_linalg_complex_LU_lndet(lu: *mut gsl_matrix_complex) -> c_double;
    pub fn gsl_linalg_LU_sgndet(lu: *mut gsl_matrix, signum: c_int) -> c_double;
    pub fn gsl_linalg_complex_LU_sgndet(lu: *mut gsl_matrix_complex, signum: c_int) -> gsl_complex;
    // QR Decomposition
    pub fn gsl_linalg_QR_decomp(a: *mut gsl_matrix, tau: *mut gsl_vector) -> c_int;
    pub fn gsl_linalg_QR_solve(
        qr: *const gsl_matrix,
        tau: *const gsl_vector,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QR_svx(
        qr: *const gsl_matrix,
        tau: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QR_lssolve(
        qr: *const gsl_matrix,
        tau: *const gsl_vector,
        b: *const gsl_vector,
        x: *mut gsl_vector,
        residual: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QR_QTvec(
        qr: *const gsl_matrix,
        tau: *const gsl_vector,
        v: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QR_Qvec(
        qr: *const gsl_matrix,
        tau: *const gsl_vector,
        v: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QR_QTmat(
        qr: *const gsl_matrix,
        tau: *const gsl_vector,
        v: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_linalg_QR_Rsolve(
        qr: *const gsl_matrix,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QR_Rsvx(qr: *const gsl_matrix, x: *mut gsl_vector) -> c_int;
    pub fn gsl_linalg_QR_unpack(
        qr: *const gsl_matrix,
        tau: *const gsl_vector,
        q: *mut gsl_matrix,
        r: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_linalg_QR_QRsolve(
        q: *mut gsl_matrix,
        r: *mut gsl_matrix,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QR_update(
        q: *mut gsl_matrix,
        r: *mut gsl_matrix,
        w: *mut gsl_vector,
        v: *const gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_R_solve(
        r: *const gsl_matrix,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_R_svx(r: *const gsl_matrix, x: *mut gsl_vector) -> c_int;
    // QR Decomposition with Column Pivoting
    pub fn gsl_linalg_QRPT_decomp(
        a: *mut gsl_matrix,
        tau: *mut gsl_vector,
        p: *mut gsl_permutation,
        signum: *mut c_int,
        norm: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QRPT_decomp2(
        a: *const gsl_matrix,
        q: *mut gsl_matrix,
        r: *mut gsl_matrix,
        tau: *mut gsl_vector,
        p: *mut gsl_permutation,
        signum: *mut c_int,
        norm: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QRPT_solve(
        qr: *const gsl_matrix,
        tau: *const gsl_vector,
        p: *const gsl_permutation,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QRPT_svx(
        qr: *const gsl_matrix,
        tau: *const gsl_vector,
        p: *const gsl_permutation,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QRPT_QRsolve(
        q: *const gsl_matrix,
        r: *const gsl_matrix,
        p: *const gsl_permutation,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QRPT_update(
        q: *const gsl_matrix,
        r: *const gsl_matrix,
        p: *const gsl_permutation,
        w: *mut gsl_vector,
        v: *const gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QRPT_Rsolve(
        qr: *const gsl_matrix,
        p: *const gsl_permutation,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_QRPT_Rsvx(
        qr: *const gsl_matrix,
        p: *const gsl_permutation,
        x: *mut gsl_vector,
    ) -> c_int;
    // Singular Value Decomposition
    pub fn gsl_linalg_SV_decomp(
        a: *mut gsl_matrix,
        v: *mut gsl_matrix,
        s: *mut gsl_vector,
        work: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_SV_decomp_mod(
        a: *mut gsl_matrix,
        x: *mut gsl_matrix,
        v: *mut gsl_matrix,
        s: *mut gsl_vector,
        work: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_SV_decomp_jacobi(
        a: *mut gsl_matrix,
        v: *mut gsl_matrix,
        s: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_SV_solve(
        u: *const gsl_matrix,
        v: *const gsl_matrix,
        s: *const gsl_vector,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_SV_leverage(u: *const gsl_matrix, h: *mut gsl_vector) -> c_int;
    // Cholesky Decomposition
    pub fn gsl_linalg_cholesky_decomp(a: *mut gsl_matrix) -> c_int;
    pub fn gsl_linalg_complex_cholesky_decomp(a: *mut gsl_matrix_complex) -> c_int;
    pub fn gsl_linalg_cholesky_solve(
        cholesky: *const gsl_matrix,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_complex_cholesky_solve(
        cholesky: *const gsl_matrix_complex,
        b: *const gsl_vector_complex,
        x: *mut gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_linalg_cholesky_svx(cholesky: *const gsl_matrix, x: *mut gsl_vector) -> c_int;
    pub fn gsl_linalg_complex_cholesky_svx(
        cholesky: *const gsl_matrix_complex,
        x: *mut gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_linalg_cholesky_invert(cholesky: *mut gsl_matrix) -> c_int;
    pub fn gsl_linalg_complex_cholesky_invert(cholesky: *mut gsl_matrix_complex) -> c_int;
    // Tridiagonal Decomposition of Real Symmetric Matrices
    pub fn gsl_linalg_symmtd_decomp(a: *mut gsl_matrix, tau: *mut gsl_vector) -> c_int;
    pub fn gsl_linalg_symmtd_unpack(
        a: *const gsl_matrix,
        tau: *const gsl_vector,
        q: *mut gsl_matrix,
        diag: *mut gsl_vector,
        subdiag: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_symmtd_unpack_T(
        a: *const gsl_matrix,
        diag: *mut gsl_vector,
        subdiag: *mut gsl_vector,
    ) -> c_int;
    // Tridiagonal Decomposition of Hermitian Matrices
    pub fn gsl_linalg_hermtd_decomp(
        a: *mut gsl_matrix_complex,
        tau: *mut gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_linalg_hermtd_unpack(
        a: *const gsl_matrix_complex,
        tau: *const gsl_vector_complex,
        u: *mut gsl_matrix_complex,
        diag: *mut gsl_vector,
        subdiag: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_hermtd_unpack_T(
        a: *const gsl_matrix_complex,
        diag: *mut gsl_vector,
        subdiag: *mut gsl_vector,
    ) -> c_int;
    // Hessenberg Decomposition of Real Matrices
    pub fn gsl_linalg_hessenberg_decomp(a: *mut gsl_matrix, tau: *mut gsl_vector) -> c_int;
    pub fn gsl_linalg_hessenberg_unpack(
        h: *mut gsl_matrix,
        tau: *mut gsl_vector,
        u: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_linalg_hessenberg_unpack_accum(
        h: *mut gsl_matrix,
        tau: *mut gsl_vector,
        v: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_linalg_hessenberg_set_zero(a: *mut gsl_matrix) -> c_int;
    // Hessenberg-Triangular Decomposition of Real Matrices
    pub fn gsl_linalg_hesstri_decomp(
        a: *mut gsl_matrix,
        b: *mut gsl_matrix,
        u: *mut gsl_matrix,
        v: *mut gsl_matrix,
        work: *mut gsl_vector,
    ) -> c_int;
    // Bidiagonalization
    pub fn gsl_linalg_bidiag_decomp(
        a: *mut gsl_matrix,
        tau_u: *mut gsl_vector,
        tau_v: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_bidiag_unpack(
        a: *const gsl_matrix,
        tau_u: *const gsl_vector,
        u: *mut gsl_matrix,
        tau_v: *const gsl_vector,
        v: *mut gsl_matrix,
        diag: *mut gsl_vector,
        superdiag: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_bidiag_unpack2(
        a: *mut gsl_matrix,
        tau_u: *mut gsl_vector,
        tau_v: *mut gsl_vector,
        v: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_linalg_bidiag_unpack_B(
        a: *const gsl_matrix,
        diag: *mut gsl_vector,
        superdiag: *mut gsl_vector,
    ) -> c_int;
    // Givens Rotations
    //pub fn gsl_linalg_givens();
    //pub fn gsl_linalg_givens_gv();
    // Householder Transformations
    pub fn gsl_linalg_householder_transform(v: *mut gsl_vector) -> c_double;
    pub fn gsl_linalg_complex_householder_transform(v: *mut gsl_vector_complex) -> gsl_complex;
    pub fn gsl_linalg_householder_hm(
        tau: c_double,
        v: *const gsl_vector,
        a: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_linalg_complex_householder_hm(
        tau: gsl_complex,
        v: *const gsl_vector_complex,
        a: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_linalg_householder_mh(
        tau: c_double,
        v: *const gsl_vector,
        a: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_linalg_complex_householder_mh(
        tau: gsl_complex,
        v: *const gsl_vector_complex,
        a: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_linalg_householder_hv(
        tau: c_double,
        v: *const gsl_vector,
        w: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_linalg_complex_householder_hv(
        tau: gsl_complex,
        v: *const gsl_vector_complex,
        w: *mut gsl_matrix_complex,
    ) -> c_int;
    // Householder solver for linear systems
    pub fn gsl_linalg_HH_solve(
        a: *mut gsl_matrix,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_HH_svx(a: *mut gsl_matrix, x: *mut gsl_vector) -> c_int;
    // Tridiagonal Systems
    pub fn gsl_linalg_solve_tridiag(
        diag: *const gsl_vector,
        e: *const gsl_vector,
        f: *const gsl_vector,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_solve_symm_tridiag(
        diag: *const gsl_vector,
        e: *const gsl_vector,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_solve_cyc_tridiag(
        diag: *const gsl_vector,
        e: *const gsl_vector,
        f: *const gsl_vector,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_linalg_solve_symm_cyc_tridiag(
        diag: *const gsl_vector,
        e: *const gsl_vector,
        b: *const gsl_vector,
        x: *mut gsl_vector,
    ) -> c_int;
    // Balancing
    pub fn gsl_linalg_balance_matrix(a: *mut gsl_matrix, d: *mut gsl_vector) -> c_int;
}

#[repr(C)]
pub struct gsl_vector_float {
    pub size: size_t,
    pub stride: size_t,
    pub data: *mut c_float,
    pub block: *mut gsl_block_float,
    pub owner: c_int,
}

#[repr(C)]
pub struct gsl_block_float {
    pub size: size_t,
    pub data: *mut c_float,
}

#[repr(C)]
pub struct gsl_vector {
    pub size: size_t,
    pub stride: size_t,
    pub data: *mut c_double,
    pub block: *mut gsl_block,
    pub owner: c_int,
}

#[repr(C)]
pub struct gsl_vector_view {
    pub vector: gsl_vector,
}

#[repr(C)]
pub struct gsl_block {
    pub size: size_t,
    pub data: *mut c_double,
}

#[repr(C)]
pub struct gsl_vector_complex_float {
    pub size: size_t,
    pub stride: size_t,
    pub data: *mut c_float,
    pub block: *mut gsl_block_complex_float,
    pub owner: c_int,
}

#[repr(C)]
pub struct gsl_block_complex_float {
    pub size: size_t,
    pub data: *mut c_float,
}

#[repr(C)]
pub struct gsl_vector_complex {
    pub size: size_t,
    pub stride: size_t,
    pub data: *mut c_double,
    pub block: *mut gsl_block_complex,
    pub owner: c_int,
}

#[repr(C)]
pub struct gsl_block_complex {
    pub size: size_t,
    pub data: *mut c_double,
}

#[repr(C)]
pub struct gsl_matrix {
    pub size1: size_t,
    pub size2: size_t,
    pub tda: size_t,
    pub data: *mut c_double,
    pub block: *mut gsl_block,
    pub owner: c_int,
}

#[repr(C)]
pub struct gsl_matrix_view {
    pub mat: gsl_matrix,
}

#[repr(C)]
pub struct gsl_matrix_float {
    pub size1: size_t,
    pub size2: size_t,
    pub tda: size_t,
    pub data: *mut c_float,
    pub block: *mut gsl_block_float,
    pub owner: c_int,
}

#[repr(C)]
pub struct gsl_matrix_complex {
    pub size1: size_t,
    pub size2: size_t,
    pub tda: size_t,
    pub data: *mut c_double,
    pub block: *mut gsl_block_complex,
    pub owner: c_int,
}

#[repr(C)]
pub struct gsl_matrix_complex_float {
    pub size1: size_t,
    pub size2: size_t,
    pub tda: size_t,
    pub data: *mut c_float,
    pub block: *mut gsl_block_complex_float,
    pub owner: c_int,
}

#[repr(C)]
pub struct gsl_eigen_symm_workspace {
    pub size: size_t,
    pub d: *mut c_double,
    pub sd: *mut c_double,
}

#[repr(C)]
pub struct gsl_eigen_symmv_workspace {
    pub size: size_t,
    pub d: *mut c_double,
    pub sd: *mut c_double,
    pub gc: *mut c_double,
    pub gs: *mut c_double,
}

#[repr(C)]
pub struct gsl_eigen_herm_workspace {
    pub size: size_t,
    pub d: *mut c_double,
    pub sd: *mut c_double,
    pub tau: *mut c_double,
}

#[repr(C)]
pub struct gsl_eigen_hermv_workspace {
    pub size: size_t,
    pub d: *mut c_double,
    pub sd: *mut c_double,
    pub gc: *mut c_double,
    pub gs: *mut c_double,
}

#[repr(C)]
pub struct gsl_eigen_francis_workspace {
    pub size: size_t,           // matrix size
    pub max_iterations: size_t, // max iterations since last eigenvalue found
    pub n_iter: size_t,         // number of iterations since last eigenvalue found
    pub n_evals: size_t,        // number of eigenvalues found so far
    pub compute_t: c_int,       // compute Schur form T = Z^t A Z
    pub H: *mut gsl_matrix,     // pointer to Hessenberg matrix
    pub Z: *mut gsl_matrix,     // pointer to Schur vector matrix
}

#[repr(C)]
pub struct gsl_eigen_nonsymm_workspace {
    pub size: size_t,          // size of matrices
    pub diag: *mut gsl_vector, // diagonal matrix elements from balancing
    pub tau: *mut gsl_vector,  // Householder coefficients
    pub Z: *mut gsl_matrix,    // pointer to Z matrix
    pub do_balance: c_int,     // perform balancing transformation?
    pub n_evals: size_t,       // number of eigenvalues found
    pub francis_workspace_p: *mut gsl_eigen_francis_workspace,
}

#[repr(C)]
pub struct gsl_eigen_nonsymmv_workspace {
    pub size: size_t,           // size of matrices
    pub work: *mut gsl_vector,  // scratch workspace
    pub work2: *mut gsl_vector, // scratch workspace
    pub work3: *mut gsl_vector, // scratch workspace
    pub Z: *mut gsl_matrix,     // pointer to Schur vectors
    pub nonsymm_workspace_p: *mut gsl_eigen_nonsymm_workspace,
}

#[repr(C)]
pub struct gsl_eigen_gensymm_workspace {
    pub size: size_t,
    pub symm_workspace_p: gsl_eigen_symm_workspace,
}

#[repr(C)]
pub struct gsl_eigen_gensymmv_workspace {
    pub size: size_t,
    pub symmv_workspace_p: gsl_eigen_symmv_workspace,
}

#[repr(C)]
pub struct gsl_eigen_genherm_workspace {
    pub size: size_t,
    pub herm_workspace_p: *mut gsl_eigen_herm_workspace,
}

#[repr(C)]
pub struct gsl_eigen_genhermv_workspace {
    pub size: size_t,
    pub hermv_workspace_p: *mut gsl_eigen_hermv_workspace,
}

#[repr(C)]
pub struct gsl_eigen_gen_workspace {
    pub size: size_t,           // size of matrices
    pub work: *mut gsl_vector,  // scratch workspace
    pub n_evals: size_t,        // number of eigenvalues found
    pub max_iterations: size_t, // maximum QZ iterations allowed
    pub n_iter: size_t,         // number of iterations since last eigenvalue found
    pub eshift: c_double,       // exceptional shift counter
    pub need_top: c_int,        // need to compute top index?
    pub atol: c_double,         // tolerance for splitting A matrix
    pub btol: c_double,         // tolerance for splitting B matrix
    pub ascale: c_double,       // scaling factor for shifts
    pub bscale: c_double,       // scaling factor for shifts
    pub H: *mut gsl_matrix,     // pointer to hessenberg matrix
    pub R: *mut gsl_matrix,     // pointer to upper triangular matrix
    pub compute_s: c_int,       // compute generalized Schur form S
    pub compute_t: c_int,       // compute generalized Schur form T
    pub Q: *mut gsl_matrix,     // pointer to left Schur vectors
    pub Z: *mut gsl_matrix,     // pointer to right Schur vectors
}

#[repr(C)]
pub struct gsl_eigen_genv_workspace {
    pub size: size_t,           // size of matrices
    pub work1: *mut gsl_vector, // 1-norm of columns of A
    pub work2: *mut gsl_vector, // 1-norm of columns of B
    pub work3: *mut gsl_vector, // real part of eigenvector
    pub work4: *mut gsl_vector, // imag part of eigenvector
    pub work5: *mut gsl_vector, // real part of back-transformed eigenvector
    pub work6: *mut gsl_vector, // imag part of back-transformed eigenvector
    pub Q: *mut gsl_matrix,     // pointer to left Schur vectors
    pub Z: *mut gsl_matrix,     // pointer to right Schur vectors
    pub gen_workspace_p: *mut gsl_eigen_gen_workspace,
}
