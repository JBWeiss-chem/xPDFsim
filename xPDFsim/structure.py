"""Minimal crystal structure: unit cell, sites with their element occupancies, neighbor lists."""

import numpy as np

from .elements import ATOMIC_NUMBERS


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

        Every site is compared with every site in all cell images that can
        lie within r; the images needed along each lattice direction follow
        from the spacing of the lattice planes (1/|b_k|, b: reciprocal vectors),
        which holds for triclinic cells too. The work is vectorised over
        all pairs of sites per image.

        Args:
            r (float): Cutoff distance in Å.

        Returns:
            numpy arrays: Indices of the centre sites, indices of the neighbour
            sites, cell offsets of the neighbours (N x 3, integers) and the
            distances.
        """
        L = self.lattice
        cart = self.coords
        spacing = 1 / np.linalg.norm(np.linalg.inv(L).T, axis=1)
        reach = np.ceil(r / spacing).astype(int)
        ranges = [np.arange(-n, n + 1) for n in reach]
        offsets = np.array(np.meshgrid(*ranges, indexing="ij")).reshape(3, -1).T
        n = len(cart)
        diff0 = cart[np.newaxis, :, :] - cart[:, np.newaxis, :]          # j - i
        centers, neighbors, images, distances = [], [], [], []
        for offset in offsets:
            d = np.linalg.norm(diff0 + offset @ L, axis=2)
            mask = d <= r
            if not offset.any():
                mask &= ~np.eye(n, dtype=bool)                           # no site with itself
            i, j = np.nonzero(mask)
            if len(i):
                centers.append(i)
                neighbors.append(j)
                images.append(np.broadcast_to(offset, (len(i), 3)))
                distances.append(d[i, j])
        if not centers:
            return np.zeros(0, int), np.zeros(0, int), np.zeros((0, 3), int), np.zeros(0)
        return (np.concatenate(centers), np.concatenate(neighbors),
                np.concatenate(images).astype(int), np.concatenate(distances))
