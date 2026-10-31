"""Velocity-magnitude slices from a Nek5000 field file (std format, axis-aligned box elements).

usage: python tools/nek_slices.py <coord_file (with X)> <field_file> <outdir>
"""
import sys, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


def read_fld(path):
    with open(path, "rb") as f:
        hdr = f.read(132).decode()
        tag = np.frombuffer(f.read(4), "<f4")[0]
        end = "<" if abs(tag - 6.54321) < 1e-4 else ">"
        tok = hdr.split()
        wd, nx, ny, nz, nel = int(tok[1]), int(tok[2]), int(tok[3]), int(tok[4]), int(tok[6])
        time = float(tok[7]); step = int(tok[8])
        vars_ = tok[11]
        assert wd == 4
        ids = np.frombuffer(f.read(4 * nel), end + "i4")
        npt = nx * ny * nz
        out = {}
        def blk(ncomp):
            a = np.frombuffer(f.read(4 * nel * npt * ncomp), end + "f4")
            return a.reshape(nel, ncomp, nz, ny, nx).astype(np.float64)
        if "X" in vars_:
            out["X"] = blk(3)
        if "U" in vars_:
            out["U"] = blk(3)
    order = np.argsort(ids)
    for k in out:
        out[k] = out[k][order]
    return out, nx, time, step


def gll(n):
    # GLL nodes on [-1,1] (n points)
    from numpy.polynomial import legendre as L
    c = np.zeros(n); c[-1] = 1
    inner = L.legroots(L.legder(c))
    return np.concatenate(([-1.0], np.sort(inner), [1.0]))


def lagrange_weights(nodes, r):
    w = np.ones(len(nodes))
    for j, xj in enumerate(nodes):
        for m, xm in enumerate(nodes):
            if m != j:
                w[j] *= (r - xm) / (xj - xm)
    return w


def slice_points(X, V, axis, val, nodes):
    """Interpolate |u| onto plane coord[axis]=val. Arrays are (nel, 3, k, j, i);
    local array axis for global coordinate c is 2-c (i->x, j->y, k->z)."""
    c = X[:, axis]
    lo, hi = c.min(axis=(1, 2, 3)), c.max(axis=(1, 2, 3))
    sel = np.where((lo <= val + 1e-9) & (hi >= val - 1e-9))[0]
    pts, vals = [], []
    la = 3 - axis  # position of that local axis in (nel, k, j, i) after picking component
    for e in sel:
        r = 2 * (val - lo[e]) / (hi[e] - lo[e]) - 1
        w = lagrange_weights(nodes, r)
        xe = np.tensordot(X[e], w, axes=([la], [0]))  # (3, n, n)
        ue = np.tensordot(V[e], w, axes=([la], [0]))
        mag = np.sqrt((ue ** 2).sum(axis=0))
        pts.append(xe.reshape(3, -1).T); vals.append(mag.ravel())
    return np.vstack(pts), np.concatenate(vals)


def main(cfile, ffile, outdir):
    xd, n, _, _ = read_fld(cfile)
    fd, n2, time, step = read_fld(ffile)
    X, V = xd["X"], fd["U"]
    nodes = gll(n)
    # sanity: element local i-axis varies in x only, etc.
    e0 = X[0]
    assert np.ptp(e0[0][0, 0, :]) > 0 and np.ptp(e0[1][0, 0, :]) < 1e-9 and np.ptp(e0[2][0, 0, :]) < 1e-9

    os.makedirs(outdir, exist_ok=True)
    vmax = None
    planes = [(1, 1.6, 0, 2, "x (m)", "z (m)", "Horizontal slice, y = 1.6 m (mid-height)", "slice_y1.6.png"),
              (0, 2.5, 2, 1, "z (m)", "y (m, vertical)", "Vertical slice, x = 2.5 m (through outlet)", "slice_x2.5.png")]
    results = []
    for axis, val, a, b, la, lb, title, fname in planes:
        P, M = slice_points(X, V, axis, val, nodes)
        results.append((P, M))
    for (axis, val, a, b, la, lb, title, fname), (P, M) in zip(planes, results):
        vmax = M.max()  # own scale per slice: early-time magnitudes differ by ~1000x
        A, B = P[:, a], P[:, b]
        ga = np.linspace(A.min(), A.max(), 600)
        gb = np.linspace(B.min(), B.max(), int(600 * np.ptp(B) / np.ptp(A)))
        GA, GB = np.meshgrid(ga, gb)
        G = griddata((A, B), M, (GA, GB), method="linear")
        if axis == 0:  # mask outside room+outlet "T": above ceiling only over the outlet
            G[(GB > 3.2 + 1e-6) & ((GA < 2.2) | (GA > 2.8))] = np.nan
        fig, ax = plt.subplots(figsize=(7, 7 * np.ptp(B) / np.ptp(A) + 1.2))
        cs = ax.contourf(GA, GB, G, levels=np.linspace(0, vmax, 41), cmap="viridis", extend="max")
        ax.set_aspect("equal"); ax.set_xlabel(la); ax.set_ylabel(lb)
        ax.set_title(f"|u| (m/s), {title}\nroom5x5, t = {time:.2f} s (step {step})", fontsize=10)
        fig.colorbar(cs, ax=ax, label="|u| (m/s)", shrink=0.85)
        fig.tight_layout(); fig.savefig(os.path.join(outdir, fname), dpi=150); plt.close(fig)
        print(fname, f"points={len(M)} max|u|={M.max():.3f} mean|u|={M.mean():.4f}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
