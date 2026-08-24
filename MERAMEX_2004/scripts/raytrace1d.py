#!/usr/bin/env python3
"""Two-point ray tracing in the layered 1-D model used for the locations.

The model is a stack of constant-velocity layers, so a ray is a straight segment
inside each layer and the travel-time integrals are exact once the layers are
made thin enough to absorb the earth-flattening transform.

Only what the tomography assessment needs: for a source at depth and a receiver
at the surface, find the first-arriving ray and return its path in the
source-receiver vertical plane.
"""
import numpy as np

R_EARTH = 6371.0


def load_layers(nll_vel_file, phase="P"):
    """(top_depth, velocity) pairs from a NonLinLoc LAYER model."""
    col = 2 if phase == "P" else 4
    out = []
    for line in open(nll_vel_file):
        if line.startswith("LAYER"):
            f = line.split()
            z = float(f[1])
            if z < 0:
                continue
            out.append((z, float(f[col])))
    return out


class Model1D:
    """Finely sampled, earth-flattened slowness profile."""

    def __init__(self, layers, zmax=400.0, dz=0.5, smooth=True):
        self.dz_sph = dz
        self.z_sph = np.arange(0.0, zmax + dz, dz)
        v = np.empty_like(self.z_sph)
        tops = np.array([l[0] for l in layers])
        vels = np.array([l[1] for l in layers])
        if smooth:
            # A stack of constant-velocity layers has no diving waves inside a
            # layer, so shooting cannot produce the head wave that is the true
            # first arrival beyond the crossover: at 150 km the constant-layer
            # model returns a 6.1 km/s crustal path instead of ~8 km/s Pn.
            # Interpolating velocity between layer tops turns every arrival into
            # a turning ray, which is also what LOTOS assumes for its reference
            # model.
            v[:] = np.interp(self.z_sph, tops, vels)
        else:
            idx = np.searchsorted(tops, self.z_sph, side="right") - 1
            idx = np.clip(idx, 0, len(vels) - 1)
            v[:] = vels[idx]
        # below the tabulated model keep the deepest velocity, with a mild
        # gradient so deep rays still turn
        deep = self.z_sph > tops[-1]
        v[deep] = vels[-1] * (1 + 0.0004 * (self.z_sph[deep] - tops[-1]))
        self.v_sph = v
        # earth flattening: a spherical layered model becomes an equivalent flat one
        r = R_EARTH - self.z_sph
        r = np.maximum(r, 1.0)
        self.z = R_EARTH * np.log(R_EARTH / r)
        self.v = v * R_EARTH / r
        self.u = 1.0 / self.v                      # flat slowness
        self.dz_flat = np.diff(self.z)
        self.u_mid = 0.5 * (self.u[:-1] + self.u[1:])

    def depth_index(self, z_sph):
        return int(np.clip(round(z_sph / self.dz_sph), 0, len(self.z_sph) - 2))

    def _integrate(self, i0, i1, p):
        """(dx, dt) accumulated over flat layers i0..i1 (exclusive), plus a flag
        for whether the ray turned before reaching i1."""
        u = self.u_mid[i0:i1]
        dz = self.dz_flat[i0:i1]
        eta2 = u * u - p * p
        turned = np.nonzero(eta2 <= 0)[0]
        n = turned[0] if len(turned) else len(u)
        if n == 0:
            return 0.0, 0.0, i0, True
        eta = np.sqrt(eta2[:n])
        dx = np.sum(p / eta * dz[:n])
        dt = np.sum(u[:n] * u[:n] / eta * dz[:n])
        return dx, dt, i0 + n, len(turned) > 0

    def curves(self, i_src, p_grid):
        """Distance and time of the up-going and down-going branches vs p."""
        up_x = np.zeros_like(p_grid); up_t = np.zeros_like(p_grid)
        dn_x = np.full_like(p_grid, np.nan); dn_t = np.full_like(p_grid, np.nan)
        for k, p in enumerate(p_grid):
            x1, t1, iend, turned = self._integrate(0, i_src, p)
            if not turned:
                up_x[k], up_t[k] = x1, t1
            else:
                up_x[k] = up_t[k] = np.nan
            # down-going: source to turning point, then all the way up
            x2, t2, iturn, turned2 = self._integrate(i_src, len(self.u_mid), p)
            if turned2 and not np.isnan(up_x[k]):
                dn_x[k] = up_x[k] + 2 * x2
                dn_t[k] = up_t[k] + 2 * t2
        return up_x, up_t, dn_x, dn_t

    def path(self, i_src, p, branch):
        """Polyline (x, z_spherical) of the ray, source first."""
        xs, zs = [0.0], [self.z_sph[i_src]]
        x = 0.0
        if branch == "down":
            i = i_src
            while i < len(self.u_mid):
                eta2 = self.u_mid[i] ** 2 - p * p
                if eta2 <= 0:
                    break
                x += p / np.sqrt(eta2) * self.dz_flat[i]
                i += 1
                xs.append(x); zs.append(self.z_sph[i])
            i_top = i
        else:
            i_top = i_src
        for i in range(i_top - 1, -1, -1):
            eta2 = self.u_mid[i] ** 2 - p * p
            if eta2 <= 0:
                break
            x += p / np.sqrt(eta2) * self.dz_flat[i]
            xs.append(x); zs.append(self.z_sph[i])
        return np.array(xs), np.array(zs)


class Tracer:
    """Caches the p-curves per source-depth bin; sources repeat a lot."""

    def __init__(self, model, np_grid=1200):
        self.m = model
        self.p_grid = np.linspace(1e-6, self.m.u.max() * 0.999999, np_grid)
        self._cache = {}

    def _curves(self, i_src):
        if i_src not in self._cache:
            self._cache[i_src] = self.m.curves(i_src, self.p_grid)
        return self._cache[i_src]

    def trace(self, z_src, dist):
        """First-arriving ray to a surface receiver at epicentral `dist` (km).

        Returns (traveltime, p, branch) or None.
        """
        i_src = self.m.depth_index(z_src)
        up_x, up_t, dn_x, dn_t = self._curves(i_src)
        best = None
        for x, t, name in ((up_x, up_t, "up"), (dn_x, dn_t, "down")):
            ok = np.isfinite(x)
            if ok.sum() < 2:
                continue
            xs, ts, ps = x[ok], t[ok], self.p_grid[ok]
            order = np.argsort(xs)
            xs, ts, ps = xs[order], ts[order], ps[order]
            if not (xs[0] <= dist <= xs[-1]):
                continue
            tt = float(np.interp(dist, xs, ts))
            pp = float(np.interp(dist, xs, ps))
            if best is None or tt < best[0]:
                best = (tt, pp, name)
        return best

    def ray_path(self, z_src, dist):
        got = self.trace(z_src, dist)
        if got is None:
            return None
        tt, p, branch = got
        i_src = self.m.depth_index(z_src)
        x, z = self.m.path(i_src, p, branch)
        if x[-1] <= 0:
            return None
        x = x * (dist / x[-1])          # small rescale onto the exact distance
        return x, z, tt
