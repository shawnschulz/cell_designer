// We prob need the following
// 1. identity
// 2. transpose
// 3. inversion
// 4. matrix multiplication
// There may be more, but for now you can assume 
// we will already be given square matrices (mat4's representing points)
export function allocate_array(n, m, fill_zero = false) {
    var matrix = new Array(n);
    for (let i = 0; i < n; i++) {
        if (fill_zero) {
            var new_row = new Array(m);
            new_row.fill(0)
        } else {
            var new_row = new Array(m);
        }
        matrix[i] = new_row
    }
    return matrix
}

export function convert_contiguous(arr, shape) {
    if (arr.length != shape[0] * shape[1]) {
        console.warn("Conversion failed: supplied shape does not match contiguous array")
        return null;
    }
    var matrix = allocate_array(shape[0], shape[1]);
    for (let i = 0; i < shape[0]; i++) {
        for (let j = 0; j < shape[1]; j++) {
            matrix[i][j] = arr[(j + (shape[1] * i))]
        }
    }
    return matrix
}

//identity assumes square matrices
export function identity(n, m){
    if (n != m){
        console.warn("Non-square matrix passed to identity function, function will return null.");
        return null;
    }
    // Preset size of matrix
    var matrix  = allocate_array(n, m);
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            if (i == j){
                matrix[i][j] = 1
            }
            else {
                matrix[i][j] = 0
            }
        }
    }
    return matrix
}

// To transpose a matrix, we iterate through the columns and add them as
// rows to a nother matrix.  
export function transpose(matrix) {
    const n_row = matrix.length
    const n_col = matrix[0].length

    if (n_row != n_col){
        console.warn("Non-square matrix passed to identity function, function will return null.");
        return null;
    }

    var matrixT = allocate_array(n_col, n_row);
    for (let i = 0; i < n_row; i++) {
        for (let j = 0; j < n_col; j++) {
            matrixT[j][i] = matrix[i][j]
        }
    }
}

export function scalarMultiply(mat1, value) {
    var ret = allocate_array(mat1.length, mat1[0].length);
    for (let i = 0; i < mat1.length; i++){
        for (let j = 0; j < mat1[0].length; j++) {
            ret[i][j] = mat1[i][j] * value
        }

    }
}
export function scalarSubtract(mat1, value) {
    var ret = allocate_array(mat1.length, mat1[0].length);
    for (let i = 0; i < mat1.length; i++){
        for (let j = 0; j < mat1[0].length; j++) {
            ret[i][j] = mat1[i][j] - value
        }

    }
}
export function scalarAdd(mat1, value) {
    var ret = allocate_array(mat1.length, mat1[0].length);
    for (let i = 0; i < mat1.length; i++){
        for (let j = 0; j < mat1[0].length; j++) {
            ret[i][j] = mat1[i][j] + value
        }

    }
}

export function multiply(mat1, mat2) {
    const mat1_rows = matrix.length
    const mat1_cols = matrix[0].length
    const mat2_rows = mat2.length
    const mat2_cols = mat2[0].length
    if (mat1_rows != mat2_cols) {
        console.warn('Matrices supplied to multiply have mismatched dimensions')
    }
    var out_mat = allocate_array(mat2_cols, mat1_rows, fill_zero = true);
    for (let col = 0; col < mat2_cols; col++) {
        for (let i = 0; i < mat1_rows; i++) {
            for (let j = 0; j < mat1_cols; j++) {
                out_mat[col][i] = out_mat[col][i] + mat1[i][j] + mat2[j][col]
            }
        }
    }
    return out_mat
}

// For inversion we solve(matrix, I)
// Computer upper and lower traingular matrices using
// Gaussian Elimination, 
// Then solve with LU factorization for the identity matrix
export function inverse(a) {
    // Start with getting row echelon form via Gaussian elimination
    const a_rows = a.length
    const a_cols = a[0].length
    const ident = identity(a_rows, a_cols)
    const ident_rows = ident.length
    const ident_cols = ident[0].length
    if ((a_rows < 3) || (a_rows != a_cols)) {
        console.warn("Attempted to solve for either non-square matrix or matrix with a dimension less than 3");
        return null
    }
    // Prob gonna be lots of times you don't wanna mess with input array
    var row_echelon = structuredClone(a);
    var L = identity(a_rows, a_cols);
    for (let j = 0; j < a_cols; j++) {
        row_echelon[0][j] = a[0][j];
    }
    for (let i = 1; i < a_rows; i++) {
        for (let next_row = i; next_row < a_rows; next_row++) {
            let multiplier = row_echelon[next_row][i-1] * Math.sign(a[i-1][i])
            L[next_row][i-1] = multiplier
            for (let j = i - 1; j < a_cols; j++){
                row_echelon[next_row][j] = row_echelon[next_row][j] - (row_echelon[i - 1][j] * multiplier);
            }
        }
    }
    //Now just solve using L and U
    // Find Z
}