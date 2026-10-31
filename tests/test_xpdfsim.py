"""Tests of the CIF reader, the PDF calculation and the API. Run with
``python -m unittest discover tests`` (no dependencies besides numpy)."""

import io
import itertools
import os
import unittest
from unittest import mock

import numpy as np

from xPDFsim import server, structure
from xPDFsim.cif import parse_symop, read_structures
from xPDFsim.main import simulate_pdf
from xPDFsim.structure import Structure, cell_matrix

EXAMPLES = os.path.join(os.path.dirname(__file__), os.pardir, "xPDFsim", "static", "examples")


def example(name):
    with open(os.path.join(EXAMPLES, name), encoding="utf-8") as fh:
        return read_structures(fh.read())[0][0][1]


def simple_cif(atoms, cell="4 4 4 90 90 90", extra=""):
    a, b, c, al, be, ga = cell.split()
    rows = "\n".join(atoms)
    return f"""data_test
_cell_length_a {a}
_cell_length_b {b}
_cell_length_c {c}
_cell_angle_alpha {al}
_cell_angle_beta {be}
_cell_angle_gamma {ga}
{extra}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_note
{rows}
"""


class TestCif(unittest.TestCase):
    def test_nacl_expansion(self):
        s = example("1000041.cif")
        self.assertEqual(len(s), 8)
        self.assertEqual(s.composition(), {"Na": 4.0, "Cl": 4.0})

    def test_quoted_values_are_not_keywords(self):
        cif = simple_cif(["Na1 Na 0 0 0 '_first'", "Cl1 Cl 0.5 0.5 0.5 'data_x'"])
        structures, _ = read_structures(cif)
        self.assertEqual(len(structures), 1)
        self.assertEqual(structures[0][1].composition(), {"Na": 1.0, "Cl": 1.0})

    def test_symop(self):
        R, t = parse_symop("-x+1/2, y-x, 2z+1/4")
        np.testing.assert_allclose(R, [[-1, 0, 0], [-1, 1, 0], [0, 0, 2]])
        np.testing.assert_allclose(t, [0.5, 0, 0.25])
        for bad in ("x,y", "x*y,y,z", "x/0,y,z", "__import__('os'),y,z"):
            with self.assertRaises(ValueError):
                parse_symop(bad)

    def test_loop_with_missing_value(self):
        cif = simple_cif(["Na1 Na 0 0 0 a", "Cl1 Cl 0.5 0.5 0.5 b", "Cl2 Cl 0.5 0.5"])
        structures, warnings = read_structures(cif)
        self.assertEqual(len(structures[0][1]), 2)
        self.assertTrue(any("not a multiple" in w for w in warnings))

    def test_invalid_cell(self):
        for cell in ((5, 5, 5, 90, 90, 0), (5, 5, 5, 120, 120, 120.1), (-1, 5, 5, 90, 90, 90)):
            with self.assertRaises(ValueError):
                cell_matrix(*cell)
        with self.assertRaises(ValueError):
            read_structures(simple_cif(["Na1 Na 0 0 0 ."], cell="4 4 4 90 90 0"))


class TestStructure(unittest.TestCase):
    def test_neighbor_list_triclinic(self):
        rng = np.random.default_rng(1)
        L = cell_matrix(3.1, 4.3, 5.2, 70, 110, 65)
        s = Structure(L, rng.random((5, 3)), [{"C": 1.0}] * 5)
        r = 7.0
        centers, neighbors, offsets, distances = s.neighbor_list(r)
        found = sorted(zip(centers, neighbors, map(tuple, offsets)))
        # brute force over a generous range of images
        expected = []
        for i, j in itertools.product(range(5), repeat=2):
            for n in itertools.product(range(-6, 7), repeat=3):
                if i == j and not any(n):
                    continue
                if np.linalg.norm(s.coords[j] + np.array(n) @ L - s.coords[i]) <= r:
                    expected.append((i, j, n))
        self.assertEqual(found, sorted(expected))

    def test_neighbor_list_cell_list(self):
        # large cells, where the sites are sorted into several bins per direction, compared with
        # all pairs over all images; also with small chunks of candidate pairs
        rng = np.random.default_rng(2)
        for cell, r in (((21, 23, 25, 90, 90, 90), 4.0), ((18, 26, 22, 75, 100, 115), 5.5),
                        ((30, 4, 30, 90, 90, 90), 6.0)):
            L = cell_matrix(*cell)
            frac = rng.random((80, 3))
            frac[0] = [1 - 1e-17, 0, 0.5]                       # rounds to the cell edge
            s = Structure(L, frac, [{"C": 1.0}] * 80)
            ranges = [range(-k, k + 1) for k in np.ceil(r * np.linalg.norm(np.linalg.inv(L), axis=0)).astype(int)]
            expected = []
            for n in itertools.product(*ranges):
                d = np.linalg.norm(s.coords[None, :, :] + np.array(n) @ L - s.coords[:, None, :], axis=2)
                for i, j in zip(*np.nonzero(d <= r)):
                    if i != j or any(n):
                        expected.append((i, j, n))
            for chunk in (structure.PAIRS_PER_CHUNK, 37):
                with mock.patch.object(structure, "PAIRS_PER_CHUNK", chunk):
                    c, nb, off, dist = s.neighbor_list(r)
                self.assertEqual(sorted(zip(c, nb, map(tuple, off))), sorted(expected))
                np.testing.assert_allclose(dist, np.linalg.norm(s.coords[nb] + off @ L - s.coords[c], axis=1))


class TestPdf(unittest.TestCase):
    def test_coordination_number(self):
        table, _, _ = simulate_pdf(example("1000041.cif"), "N_r", 10, 0.01, 0.1, True)
        self.assertAlmostEqual(table["total"][np.searchsorted(table["radii"], 3.0)], 6.0, places=6)

    def test_partials_add_up(self):
        table, _, _ = simulate_pdf(example("9015662.cif"), "G_r", 10, 0.01, 0.1, False)
        partials = sum(v for k, v in table.items() if k not in ("radii", "total"))
        np.testing.assert_allclose(partials, table["total"], atol=1e-9)

    def test_baseline_not_broadened(self):
        # below the first peak the broadened functions follow their low-r limits exactly
        s = example("1000041.cif")
        for conv, limit in (("GK_r", lambda r, rho, b2: -b2 + 0 * r),
                            ("G_r", lambda r, rho, b2: -4 * np.pi * rho * r)):
            table, rho, b2 = simulate_pdf(s, conv, 10, 0.01, 0.1, False)
            low = table["radii"] < 1.5
            np.testing.assert_allclose(table["total"][low], limit(table["radii"][low], rho, b2), atol=1e-9)

    def test_broadening_near_r_max(self):
        # with a wide Gaussian, the function at r_max includes the peaks beyond it
        s = example("1000041.cif")
        short, _, _ = simulate_pdf(s, "G_r", 8, 0.01, 0.4, False)
        long, _, _ = simulate_pdf(s, "G_r", 14, 0.01, 0.4, False)
        np.testing.assert_allclose(short["total"], long["total"][:len(short["total"])], atol=1e-6)

    def test_small_partials_keep_precision(self):
        table, _, _ = simulate_pdf(example("4002052.cif"), "rho_r", 20, 0.01, 0.1, False)
        tail = table["radii"] > 10
        self.assertGreater(table["H-H"][tail].mean(), 1e-6)


class TestServer(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(EXAMPLES, "1000041.cif"), encoding="utf-8") as fh:
            self.text = fh.read()

    def test_load(self):
        out = server.api_load({"text": self.text, "r_max": 10, "dr": 0.01})
        self.assertEqual(out["structures"][0]["structure"]["formula"], "NaCl")

    def request(self, path, method="GET", **headers):
        environ = {"REQUEST_METHOD": method, "PATH_INFO": path, "HTTP_HOST": "127.0.0.1:8766",
                   "CONTENT_LENGTH": "2", "wsgi.input": io.BytesIO(b"{}")}
        environ.update({f"HTTP_{k.upper()}": v for k, v in headers.items()})
        status = []
        server.app(environ, lambda s, h: status.append(s))
        return status[0]

    def test_host_and_origin(self):
        with mock.patch.object(server, "ALLOWED_HOSTS", set(server.LOOPBACK)):
            self.assertTrue(self.request("/api/startup").startswith("200"))
            self.assertTrue(self.request("/api/startup", host="localhost:8766").startswith("200"))
            self.assertTrue(self.request("/api/startup", host="evil.example:8766").startswith("403"))
            self.assertTrue(self.request("/api/pdf", "POST", origin="http://127.0.0.1:8766").startswith("400"))
            self.assertTrue(self.request("/api/pdf", "POST", origin="https://evil.example").startswith("403"))

    def test_bad_block(self):
        for block in (5, "x", None):
            with self.assertRaises(server.ClientError):
                server.api_pdf({"text": self.text, "block": block})


if __name__ == "__main__":
    unittest.main()
