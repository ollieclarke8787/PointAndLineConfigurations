
-- 12 points in a grid point-line configuration

load "gridGens.m2"
-- Defines the ring R with variables: x_1 .. x_12, y_1 .. y_12, z_1 .. z_12
-- Defines the ideal P2 of R - which is a prime ideal, whose generators are
-- taken from the paper:
-- [On the primary decomposition of some determinantal hyperedge ideal - Pfister and Steenpass 2021]


-- matrix of variables with basis
X = transpose genericMatrix(R, 12, 3) | id_(R^3)
dets = (gens P2)_{0 .. 15} -- 3x3 determinantal generators of P2

-- basis are columns 12, 13, 14 of X

-- grid shape:
-- 0  3  6  9
-- 1  4  7  10
-- 2  5  8  11

QSgen = method()
QSgen(List) := L -> (
    -- L should be a list of elements from {1,2,3}
    K := apply(L, x -> x+11);
    sum{
	 (det X_{0, 3, K_0})*(det X_{1, 7,  K_1})*(det X_{2, 11, K_2})*(det X_{4, 5, K_3})*(det X_{6, 8, K_4})*(det X_{9, 10,  K_5}),
	 (det X_{0, 6, K_0})*(det X_{1, 10, K_1})*(det X_{2, 5,  K_2})*(det X_{3, 4, K_3})*(det X_{7, 8, K_4})*(det X_{9, 11,  K_5}),
	 (det X_{0, 9, K_0})*(det X_{1, 4,  K_1})*(det X_{2, 8,  K_2})*(det X_{3, 5, K_3})*(det X_{6, 7, K_4})*(det X_{10, 11, K_5}),
	-(det X_{0, 3, K_0})*(det X_{1, 10, K_1})*(det X_{2, 8,  K_2})*(det X_{4, 5, K_3})*(det X_{6, 7, K_4})*(det X_{9, 11,  K_5}),
	-(det X_{0, 6, K_0})*(det X_{1, 4,  K_1})*(det X_{2, 11, K_2})*(det X_{3, 5, K_3})*(det X_{7, 8, K_4})*(det X_{9, 10,  K_5}),
	-(det X_{0, 9, K_0})*(det X_{1, 7,  K_1})*(det X_{2, 5,  K_2})*(det X_{3, 4, K_3})*(det X_{6, 8, K_4})*(det X_{10, 11, K_5}) 
	}
    )

-- for a general QS-polynomial we swap the role of points 1,2,3 with the points in another column
QSgen(ZZ, List) := (i, L) -> (
    -- i is a column index: 0,1,2,3
    -- L should be a list of elements from {1,2,3}
    f := QSgen L;
    sub(f, flatten for j from 1 to 12 list (
	    if 1 <= j and j <= 3 then {x_j => x_(j+3*i), y_j => y_(j+3*i), z_j => z_(j+3*i)}
	    else if 3*i+1 <= j and j <= 3*i+3 then {x_j => x_(j-3*i), y_j => y_(j-3*i), z_j => z_(j-3*i)}
	    else {x_j => x_j, y_j => y_j, z_j => z_j}
	    )
	)
    )

QSInds = reverse for s in subsets(8, 2) list (
    c1 := s_0;
    c2 := s_1 - s_0 - 1;
    c3 := 7 - s_1;
    toList((c1 : 1) | (c2 : 2) | (c3 : 3))
    )

QSGenMatrix = matrix {for L in QSInds list QSgen L};

P2' = ideal(dets | QSGenMatrix);


print "-- inclusion of P2' inside P2"
for L in QSInds do (
    print(L, zero(QSgen L % P2))
    )

print "-- inclusion of P2 inside P2'"
for f in first entries gens P2 do (
    print zero (f % P2')
    )

end --
-- All generators of each ideal reduce to zero so the ideals are equal

-- check the degree limit that the GB computations needed
print "-- P2 GB computation:"
print gb(P2, DegreeLimit => 1)

print ""
print "-- P2' GB computation:"
print gb(P2', DegreeLimit => 1)

-- GB computation is needed up to degree 12

-- Check that the general QS polynomial belong to P2 (= P2'):
-- the following will give produce an error is one of the QS polynomials doesn't belong
S = (S := set {1,2,3}; for i from 1 to 5 do S = (S ** set {1,2,3}) / splice; S);

for i from 1 to 3 do (
    for s in keys S do (
	sList := toList s;
	ans := zero (QSgen(i, sList) % P2);
	print ("G_3,4(c_"| toString (i+1) | "," | toString s | ") in P2 : " | toString ans);
	if not ans then error("QS polynomial not in ideal!");
	)
    )


end --

i0 = n -> (
    R := QQ[x_1 .. x_n, y_1 .. y_n, z_1 .. z_n];
    X := transpose genericMatrix(R, n, 3);
    ideal for s in subsets(n, 3) list det(X_s)
    )

for i from 3 to 20 do (
    print(i, dim i0 i);
    )
dim (i0 3)
dim (ring i0 3)

R = QQ[x_(1,1) .. x_(3,9)]
X = transpose genericMatrix(R, 9, 3)
I = ideal(det X_{0,1,2}, det X_{3,4,5}, det X_{6,7,8}, 
    det X_{0,3,6}, det X_{1,4,7}, det X_{2,5,8})
dim I


-----------------------------
M = matrix {
    {x_2 - x_3, -x_1 + x_3, x_1 - x_2,  0, 0, 0,  0, 0, 0,  0, 0, 0},
    {0, 0, 0,  x_5 - x_6, -x_4 + x_6, x_4 - x_5,  0, 0, 0,  0, 0, 0},
    {0, 0, 0,  0, 0, 0,  x_8 - x_9, -x_7 + x_9, x_7 + x_8,  0, 0, 0},
    {0, 0, 0,  0, 0, 0,  0, 0, 0,  x_11 - x_12, -x_10 + x_12, x_10 - x_11},
    {x_4 - x_7, 0, 0,  -x_1 + x_7, 0, 0,  x_1 - x_4, 0, 0,  0, 0, 0},
    {x_4 - x_10, 0, 0,  -x_1 + x_10, 0, 0,  0, 0, 0,  x_1 - x_4, 0, 0},
    {x_7 - x_10, 0, 0,  0, 0, 0,  -x_1 + x_10, 0, 0,  x_1 - x_7, 0, 0},
    {0, 0, 0,  x_7 - x_10, 0, 0,  -x_4 + x_10, 0, 0,  x_4 - x_7, 0, 0},
    {0, x_5 - x_8, 0,   0, -x_2 + x_8, 0,   0, x_2 - x_5, 0,  0, 0, 0},
    {0, x_5 - x_11, 0,  0, -x_2 + x_11, 0,  0, 0, 0,  0, x_2 - x_5, 0},
    {0, x_8 - x_11, 0,  0, 0, 0,  0, -x_2 + x_11, 0,  0, x_2 - x_8, 0},
    {0, 0, 0,  0, x_8 - x_11, 0,  0, -x_5 + x_11, 0,  0, x_5 - x_8, 0},
    {0, 0, x_6 - x_9,  0, 0, -x_3 + x_9,  0, 0, x_3 - x_6,  0, 0, 0},
    {0, 0, x_6 - x_12,  0, 0, -x_3 + x_12,  0, 0, 0,  0, 0, x_3 - x_6},
    {0, 0, x_9 - x_12,  0, 0, 0,  0, 0, -x_3 + x_12,  0, 0, x_3 - x_9},
    {0, 0, 0,  0, 0, x_9 - x_12,  0, 0, -x_6 + x_12,  0, 0, x_6 - x_9}
    } 

f = det M^{0,1,2,3, 4,5, 8,9, 12,13}_{0,1,2,3,4,5,6,7,8,9};
f % P2'

-----------------
load "gridExample.m2"
I = ideal dets
J = radical I;
