setClassUnion(name = 'SMatrix', members = c("dgCMatrix", "dgTMatrix", "matrix"))
setClassUnion("MatrixOrNull", members = c("dgTMatrix", "dgCMatrix", "matrix", "missing"))
