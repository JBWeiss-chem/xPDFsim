"""Minimal crystal structure: unit cell, sites with their element occupancies, neighbor lists."""

import numpy as np

from .elements import ATOMIC_NUMBERS

MAX_BINS_PER_AXIS = 64          # bins of the cell list along every lattice direction
PAIRS_PER_CHUNK = 2_000_000     # candidate pairs compared at once (memory of the neighbor list)


def cell_matrix(a, b, c, alpha, beta, gamma):
    """Lattice vectors as rows (a along x, b in the xy plane), lengths in Å, angles in degrees.

    Raises:
        ValueError: If the parameters do not describe a cell with a volume.
    """
    if not all(np.isfinite([a, b, c, alpha, beta, gamma])) or min(a, b, c) <= 0:
        raise ValueError("the cell lengths must be positive numbers")
    if not all(0 < x < 180 for x in (alpha, beta, gamma)):
        raise ValueError("the cell angles must be between 0 and 180°")
    al, be, ga = np.radians([alpha, beta, gamma])
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz2 = c * c - cx * cx - cy * cy
    if cz2 <= 1e-6 * c * c:                     # e.g. alpha > beta + gamma
        raise ValueError(f"the cell angles {alpha:g}, {beta:g}, {gamma:g}° give no volume")
    cz = np.sqrt(cz2)
    return np.array([[a, 0.0, 0.0], [b * np.cos(ga), b * np.sin(ga), 0.0], [cx, cy, cz]])


class Structure:
    """A periodic structure.

    Args:
        lattice (numpy array): Lattice vectors as rows (3x3, Å).
        frac (numpy array): Fractional coordinates of the sites (N x 3).
        species (list): Per site a dictionary element symbol -> occupancy;
        sites with several elements or occupancies below 1 are disordered.
        labels (list): Site labels.
        name (string): Name of the structure.
        space_group (string): Space-group symbol (for display), or None.
    """

    def __init__(self, lattice, frac, species, labels=None, name="", space_group=None):
        self.lattice = np.asarray(lattice, float)
        self.frac = np.asarray(frac, float).reshape(-1, 3) % 1.0
        self.species = [dict(sp) for sp in species]
        self.labels = list(labels) if labels is not None else [max(sp, key=sp.get) for sp in self.species]
        self.name = name
        self.space_group = space_group

    def __len__(self):
        return len(self.species)

    @property
    def coords(self):
        """Cartesian coordinates of the sites (N x 3, Å)."""
        return self.frac @ self.lattice

    @property
    def volume(self):
        return abs(np.linalg.det(self.lattice))

    @property
    def abc(self):
        return tuple(np.linalg.norm(self.lattice, axis=1))

    @property
    def angles(self):
        a, b, c = self.lattice
        angle = lambda u, v: np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))
        return angle(b, c), angle(a, c), angle(a, b)

    @property
    def is_ordered(self):
        """True if every site holds one element with full occupancy."""
        return all(len(sp) == 1 and next(iter(sp.values())) > 0.999 for sp in self.species)

    def composition(self):
        """Amount of every element in the cell (occupancies summed), in order of appearance."""
        amounts = {}
        for sp in self.species:
            for el, occ in sp.items():
                amounts[el] = amounts.get(el, 0.0) + occ
        return amounts

    def num_atoms(self):
        """Number of atoms in the cell (sum of all occupancies)."""
        return sum(self.composition().values())

    def site_weights(self, lengths):
        """Per site the sum over its elements of occupancy times a per-element value."""
        return np.array([sum(lengths[el] * occ for el, occ in sp.items()) for sp in self.species])

    def atomic_numbers(self):
        return {el: ATOMIC_NUMBERS[el] for el in self.composition()}

    def neighbor_list(self, r):
        """All pairs of sites closer than r, in both directions.

        The sites are sorted into a grid of bins (a cell list) whose thickness
        along every lattice direction is at least r where the cell allows it,
        so only sites in nearby bins are compared; the bins needed along each
        direction follow from the spacing of the lattice planes (1/|b_k|,
        b: reciprocal vectors), which holds for triclinic cells too. Bins
        beyond the cell edge are its periodic images, which also covers cells
        thinner than r. The work grows linearly with the number of sites and
        is vectorised over all candidate pairs, in chunks of bounded memory.

        Args:
            r (float): Cutoff distance in Å.

        Returns:
            numpy arrays: Indices of the centre sites, indices of the neighbour
            sites, cell offsets of the neighbours (N x 3, integers) and the
            distances.
        """
        L = self.lattice
        cart = self.coords
        n = len(cart)
        spacing = 1 / np.linalg.norm(np.linalg.inv(L).T, axis=1)
        bins = np.clip(np.floor(spacing / r), 1, MAX_BINS_PER_AXIS).astype(int)
        reach = np.ceil(r * bins / spacing).astype(int)
        # bin of every site; frac % 1 can round to exactly 1.0
        b = np.minimum(np.floor(self.frac * bins).astype(int), bins - 1)
        flat = lambda idx: (idx[:, 0] * bins[1] + idx[:, 1]) * bins[2] + idx[:, 2]
        site_bin = flat(b)
        order = np.argsort(site_bin, kind="stable")
        counts = np.bincount(site_bin, minlength=int(np.prod(bins)))
        starts = np.cumsum(counts) - counts

        ranges = [np.arange(-k, k + 1) for k in reach]
        shifts = np.array(np.meshgrid(*ranges, indexing="ij")).reshape(3, -1).T
        centers, neighbors, images, distances = [], [], [], []
        for shift in shifts:
            unwrapped = b + shift
            offset = np.floor_divide(unwrapped, bins)              # cell image of the target bin
            target = flat(unwrapped - offset * bins)
            n_cand = counts[target]
            # chunks of centre sites with at most PAIRS_PER_CHUNK candidate pairs
            cum = np.cumsum(n_cand)
            cuts = np.searchsorted(cum, np.arange(PAIRS_PER_CHUNK, cum[-1] if n else 0, PAIRS_PER_CHUNK))
            for lo, hi in zip(np.r_[0, cuts], np.r_[cuts, n]):
                c = n_cand[lo:hi]
                total = int(c.sum())
                if not total:
                    continue
                i = np.repeat(np.arange(lo, hi), c)
                first = np.cumsum(c) - c
                j = order[np.arange(total) - np.repeat(first - starts[target[lo:hi]], c)]
                off = offset[i]
                d = np.linalg.norm(cart[j] + off @ L - cart[i], axis=1)
                keep = (d <= r) & ((i != j) | off.any(axis=1))      # no site with itself
                centers.append(i[keep])
                neighbors.append(j[keep])
                images.append(off[keep])
                distances.append(d[keep])
        if not centers:
            return np.zeros(0, int), np.zeros(0, int), np.zeros((0, 3), int), np.zeros(0)
        return (np.concatenate(centers), np.concatenate(neighbors),
                np.concatenate(images).astype(int), np.concatenate(distances))
