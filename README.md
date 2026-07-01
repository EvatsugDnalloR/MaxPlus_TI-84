# MaxPlus_TI-84: Max‑Plus Algebra Toolkit for TI‑84 Plus CE‑T
**MaxPlus_TI-84** is a lightweight Python script that brings max‑plus algebra to your TI‑84 Plus CE‑T calculator (or any Python 3 environment). It handles scalars, matrices, vectors, computes eigenvalues (maximum cycle means).

---

## ✨ Features
| Option | Operation | Description |
|--------|-----------|-------------|
| **1** | Scalar Add & Mult | $a \oplus b := \max\{a, b\}$ and $a \otimes b := a + b$ for all $a , b \in \mathbb{T} := \mathbb{R} \cup (-\infty)$ |
| **2** | Matrix Addition | $(\mathbf{A} \oplus \mathbf{B})_{ij} = A_{ij} \oplus B_{ij}$ for all $\mathbf{A}, \mathbf{B} \in \mathbb{T}^{m \times n}$|
| **3** | Matrix Multiplication | $(\mathbf{A} \otimes \mathbf{B})_{ij} = \bigoplus\limits_{k = 1}^{n} A_{ik} \otimes B_{kj}$ for all $\mathbf{A} \in \mathbb{T}^{m \times n}, \mathbf{B} \in \mathbb{T}^{n \times \ell}$|
| **4** | Vector Convolution | $(\mathbf{v} \otimes \mathbf{w})_k := \bigoplus\limits_{i=1}^n v_i \otimes w_{k-i} $ for all $\mathbf{v}, \mathbf{w} \in \mathbb{T}^n$   |
| **5** | Spectral Radius (Largest Eigenvalue) |**maximum cycle mean** of the precedence graph for a given matrix $\mathbf{M} \in \mathbb{T}^{n \times n}$ (Karp’s algorithm) |

---

## 🚀 Usage

### On the Calculator (TI‑84 Plus CE‑T)
1. **Transfer** the script to the calculator.
2. **Run** it from the calculator’s Python environment.
3. **Follow** the on‑screen menu:
   - Enter scalars, matrices, or vectors in the required format.
   - Use **`i`** (case‑insensitive) as shorthand for `-inf`, e.g., `1 i 3` representing row vector $\begin{pmatrix} 1 &  -\infty &  3\end{pmatrix}$.
   - Matrices: rows separated by commas, entries by spaces.  
     Example: `1 i 3, i 2 4, 5 6 i` representing $3\times 3$-matrix $\begin{pmatrix} 1 & -\infty & 3 \\ -\infty & 2 & 4 \\ 5 & 6 & - \infty  \end{pmatrix}$

### In a Desktop Python Environment
```bash
uv run max_plus.py
```

---

## 📝 To‑Do List

| Progress | Task | Description |
| --- | --- | --- |
| ⏳ *In progress* | **Eigenvector** | Implement computation of a **normal eigenvector** with respect to the largest eigenvalue (critical graph + potential propagation). |
| 🔄 *Planned* | **Power Method** | Add an iterative (Krylov) method for approximate spectral radius and eigenvectors. |
| 🔄 *Planned* | **Input Validation** | Add better error handling and dimension checks. |
| 🔄 *Planned* | **Sparse Representation** | Optimise for larger matrices by storing only finite entries. |
