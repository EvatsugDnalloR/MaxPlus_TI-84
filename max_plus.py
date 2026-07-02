# MAXPLUS PROGRAM WITH TRUE -INF SUPPORT, CONVOLUTION, AND EIGENVALUE (FIXED)

def parse_vector(input_str):
    """Parse a space-separated vector, treating 'i' as -inf."""
    elements = input_str.strip().split()
    vec = []
    for x in elements:
        xl = x.lower()
        if xl == "i" or "inf" in xl:
            vec.append("-inf")
        else:
            vec.append(float(x))
    return vec

def parse_matrix(input_str):
    """Parse comma-separated rows, each row space-separated; 'i' means -inf."""
    matrix = []
    rows = input_str.split(",")
    for r in rows:
        elements = r.strip().split()
        row = []
        for x in elements:
            xl = x.lower()
            if xl == "i" or "inf" in xl:
                row.append("-inf")
            else:
                row.append(float(x))
        matrix.append(row)
    return matrix

def mp_add(a, b):
    if a == "-inf": return b
    if b == "-inf": return a
    return max(a, b)

def mp_mult(a, b):
    if a == "-inf" or b == "-inf":
        return "-inf"
    return a + b

def mat_add(A, B):
    r, c = len(A), len(A[0])
    C = [["-inf"]*c for _ in range(r)]
    for i in range(r):
        for j in range(c):
            C[i][j] = mp_add(A[i][j], B[i][j])
    return C

def mat_mult(A, B):
    rA, cA = len(A), len(A[0])
    rB, cB = len(B), len(B[0])
    if cA != rB:
        print("Dim Mismatch")
        return None
    C = [["-inf"]*cB for _ in range(rA)]
    for i in range(rA):
        for j in range(cB):
            for k in range(cA):
                term = mp_mult(A[i][k], B[k][j])
                C[i][j] = mp_add(C[i][j], term)
    return C

def mp_convolution(a, b):
    """Max-plus convolution of two vectors (lists of numbers or '-inf')."""
    n, m = len(a), len(b)
    c = ["-inf"] * (n + m - 1)
    for k in range(n + m - 1):
        best = "-inf"
        for i in range(n):
            j = k - i
            if 0 <= j < m:
                term = mp_mult(a[i], b[j])
                best = mp_add(best, term)
        c[k] = best
    return c

def max_cycle_mean(M):
    """
    Compute the maximum cycle mean (largest eigenvalue) of a max-plus matrix M.
    Uses Karp's algorithm correctly: max over v of min over k.
    Returns a float or '-inf' if no cycle exists.
    """
    n = len(M)
    if n == 0:
        return "-inf"

    dp = [["-inf"] * n for _ in range(n + 1)]
    for i in range(n):
        dp[0][i] = 0.0

    for k in range(1, n + 1):
        for j in range(n):
            best = "-inf"
            for i in range(n):
                if dp[k-1][i] != "-inf" and M[i][j] != "-inf":
                    val = dp[k-1][i] + M[i][j]
                    best = mp_add(best, val)
            dp[k][j] = best

    max_mean = "-inf"
    for v in range(n):
        if dp[n][v] == "-inf":
            continue
        min_mean = None
        for k in range(n):
            if dp[k][v] != "-inf":
                mean = (dp[n][v] - dp[k][v]) / (n - k)
                if min_mean is None or mean < min_mean:
                    min_mean = mean
        if min_mean is not None:
            if max_mean == "-inf" or min_mean > max_mean:
                max_mean = min_mean

    return max_mean


def compute_eigenvector(M, lam):
    """
    Compute a normal eigenvector of max-plus matrix M with eigenvalue lam.
    
    Returns a list (vector) of floats / '-inf' with max element == 0.
    """
    n = len(M)
    if n == 0:
        return []

    # ---- Step 1: build lam-shifted adjacency matrix ----
    # Definition C.31: edge from x_j to x_i with weight M[i][j] (when != -inf).
    # So adj[from][to] = adj[j][i] = M[i][j] - lam.
    # All D[i][i] start as -inf (no phantom zero self-loops).
    D = [["-inf"] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if M[i][j] != "-inf":
                D[j][i] = M[i][j] - lam   # edge from j to i

    # ---- Step 2: Floyd-Warshall (max-plus) longest paths ----
    for k in range(n):
        for i in range(n):
            if D[i][k] == "-inf":
                continue
            for j in range(n):
                if D[k][j] == "-inf":
                    continue
                val = D[i][k] + D[k][j]
                if D[i][j] == "-inf" or val > D[i][j]:
                    D[i][j] = val

    # ---- Step 3: find a node on a critical cycle (D[i][i] ≈ 0) ----
    # Without phantom self-loops, D[i][i] is the longest genuine cycle.
    # A critical cycle in the lam-shifted graph has total weight 0.
    root = None
    for i in range(n):
        if D[i][i] != "-inf" and abs(D[i][i]) < 1e-10:
            root = i
            break

    # Fallback: any node that has any cycle (shouldn't normally be needed)
    if root is None:
        for i in range(n):
            if D[i][i] != "-inf":
                root = i
                break

    if root is None:
        return ["-inf"] * n

    # ---- Step 4: eigenvector = longest paths from root ----
    eig = [D[root][j] if j != root else 0.0 for j in range(n)]
    # (D[root][root] should already be ≈0 if root is on a critical cycle,
    #  but we force it to exactly 0 for the empty-path interpretation.)

    # ---- Step 5: normalize so max element == 0 ----
    max_val = None
    for v in eig:
        if v != "-inf":
            if max_val is None or v > max_val:
                max_val = v
    if max_val is not None:
        for j in range(n):
            if eig[j] != "-inf":
                eig[j] = eig[j] - max_val

    return eig


# ---------- Main Menu ----------
print("--- MAX-PLUS ALGEBRA ---")
print("1: Scalar Add & Mult")
print("2: Matrix Add")
print("3: Matrix Mult")
print("4: Vector Convolution")
print("5: Eigenvalue (Max Cycle Mean)")
choice = int(input("Select (1-5): "))

if choice == 1:
    raw_x = input("x = ").lower()
    raw_y = input("y = ").lower()
    x = "-inf" if raw_x == "i" or "inf" in raw_x else float(raw_x)
    y = "-inf" if raw_y == "i" or "inf" in raw_y else float(raw_y)
    print("Add (Max):", mp_add(x, y))
    print("Mult (Plus):", mp_mult(x, y))

elif choice == 2 or choice == 3:
    print("Format: rows separated by commas, entries by spaces; use 'i' for -inf")
    strA = input("Matrix A: ")
    strB = input("Matrix B: ")
    matA = parse_matrix(strA)
    matB = parse_matrix(strB)
    if choice == 2:
        res = mat_add(matA, matB)
        print("Result +:")
    else:
        res = mat_mult(matA, matB)
        print("Result *:")
    if res:
        for row in res:
            print(row)

elif choice == 4:
    print("Enter vectors as space-separated numbers; use 'i' for -inf")
    vecA = input("Vector A: ")
    vecB = input("Vector B: ")
    a = parse_vector(vecA)
    b = parse_vector(vecB)
    res = mp_convolution(a, b)
    print("Convolution (max-plus):")
    print(res)

elif choice == 5:
    print("Enter square matrix (rows separated by commas, entries by spaces; use 'i' for -inf)")
    strM = input("Matrix M: ")
    M = parse_matrix(strM)
    if len(M) == 0 or len(M[0]) != len(M):
        print("Matrix must be square!")
    else:
        lam = max_cycle_mean(M)
        print("Largest eigenvalue (max cycle mean):", lam)
        if lam != "-inf":
            vec = compute_eigenvector(M, lam)
            print("Normal eigenvector:", vec)

else:
    print("Invalid choice")
