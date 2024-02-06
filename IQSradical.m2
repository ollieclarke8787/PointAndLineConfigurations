needsPackage "MinimalPrimes"
needsPackage "Matrices"

R = QQ[x_1 .. x_6, y_1 .. y_6, z_1 .. z_6]
X = (transpose genericMatrix(R, 6, 3)) | (id_(R^3))
g1 = det submatrix(X,{0,4,6}) * det submatrix(X,{1,5,6}) * det submatrix(X,{2,3,6}) - det submatrix(X,{0,5,6}) * det submatrix(X,{1,3,6}) * det submatrix(X,{2,4,6});
g2 = det submatrix(X,{0,4,6}) * det submatrix(X,{1,5,6}) * det submatrix(X,{2,3,7}) - det submatrix(X,{0,5,6}) * det submatrix(X,{1,3,6}) * det submatrix(X,{2,4,7}); 
g3 = det submatrix(X,{0,4,6}) * det submatrix(X,{1,5,6}) * det submatrix(X,{2,3,8}) - det submatrix(X,{0,5,6}) * det submatrix(X,{1,3,6}) * det submatrix(X,{2,4,8});
g4 = det submatrix(X,{0,4,6}) * det submatrix(X,{1,5,7}) * det submatrix(X,{2,3,7}) - det submatrix(X,{0,5,6}) * det submatrix(X,{1,3,7}) * det submatrix(X,{2,4,7}); 
g5 = det submatrix(X,{0,4,6}) * det submatrix(X,{1,5,7}) * det submatrix(X,{2,3,8}) - det submatrix(X,{0,5,6}) * det submatrix(X,{1,3,7}) * det submatrix(X,{2,4,8}); 
g6 = det submatrix(X,{0,4,6}) * det submatrix(X,{1,5,8}) * det submatrix(X,{2,3,8}) - det submatrix(X,{0,5,6}) * det submatrix(X,{1,3,8}) * det submatrix(X,{2,4,8});
g7 = det submatrix(X,{0,4,7}) * det submatrix(X,{1,5,7}) * det submatrix(X,{2,3,7}) - det submatrix(X,{0,5,7}) * det submatrix(X,{1,3,7}) * det submatrix(X,{2,4,7});
g8 = det submatrix(X,{0,4,7}) * det submatrix(X,{1,5,7}) * det submatrix(X,{2,3,8}) - det submatrix(X,{0,5,7}) * det submatrix(X,{1,3,7}) * det submatrix(X,{2,4,8});
g9 = det submatrix(X,{0,4,7}) * det submatrix(X,{1,5,8}) * det submatrix(X,{2,3,8}) - det submatrix(X,{0,5,7}) * det submatrix(X,{1,3,8}) * det submatrix(X,{2,4,8});
g10 = det submatrix(X,{0,4,8}) * det submatrix(X,{1,5,8}) * det submatrix(X,{2,3,8}) - det submatrix(X,{0,5,8}) * det submatrix(X,{1,3,8}) * det submatrix(X,{2,4,8});
 
G = ideal(det submatrix(X,{0,1,2}), det submatrix(X,{0,4,5}), det submatrix(X,{1,3,5}), det submatrix(X,{2,3,4}),g1,g2,g3,g4,g5,g6,g7,g8,g9,g10);
-- lead term is square free so the ideal G is radical
leadTerm G
