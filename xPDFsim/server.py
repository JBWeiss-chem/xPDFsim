"""Graphical user interface: JSON API + static three.js front end.

Start it with ``xPDFsim-gui`` (or ``python -m xPDFsim.server``); the command
line tool ``xPDFsim`` also uses it to show its result (see PRELOAD). The server
is stateless: the browser keeps the CIF text and sends it again whenever a PDF
or the list of distances has to be calculated. The X-ray weighted R(r)
histograms are calculated here, together with the X-ray and neutron pair
weights; switching the PDF convention and radiation, broadening and plotting
happen in the browser.
"""

import argparse
import functools
import json
import mimetypes
import os
import threading
import time
import webbrowser
from socketserver import ThreadingMixIn
from urllib.parse import urlsplit
from wsgiref.simple_server import WSGIRequestHandler, WSGIServer, make_server

import math

import numpy as np

from .cif import read_structures
from .elements import ATOMIC_NUMBERS, COVALENT_RADII, D_NATURAL, color, is_metal
from .pdf import NEUTRON_B, gen_R_r, pair_weights, partial_labels, scattering_length

STATIC = os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "static"))
MAX_BODY = 50 * 1024 * 1024     # bytes of a request
MAX_R = 100.0                   # Angstrom
MAX_BINS = 200_000
SIGMA_MAX = 0.4                 # Angstrom, largest broadening of the GUI (SIGMA_MAX in app.js)
R_EXTRA = 6 * SIGMA_MAX         # the histograms go beyond r_max by the reach of the Gaussians
BOND_TOLERANCE = 0.45           # Angstrom added to the sum of covalent radii
MIN_BOND = 0.6                  # shorter contacts are alternative (disordered) positions
TYPE_WINDOW = 0.02              # Angstrom, windows in which the strongest distance types are kept
TYPES_PER_WINDOW = 12


class ClientError(Exception):
    """Error caused by the request; its message is shown to the user."""


def read_cif(text):
    """Reads every structure (data block) of a CIF file, see cif.read_structures()."""
    try:
        return read_structures(text)
    except ValueError as e:
        raise ClientError(str(e))


def formula(structure):
    """Formula of the cell, elements in order of appearance, reduced to the smallest whole
    numbers if the amounts are whole numbers (e.g. 'NaCl', 'SrTiO3', 'Co16H99.824C113.28N64')."""
    amounts = structure.composition()
    if all(abs(v - round(v)) < 1e-3 for v in amounts.values()):
        divisor = functools.reduce(math.gcd, [int(round(v)) for v in amounts.values()])
        amounts = {el: round(v) // divisor for el, v in amounts.items()}
    text = lambda v: "" if v == 1 else (str(int(v)) if float(v).is_integer() else f"{v:.6g}")
    return "".join(f"{el}{text(v)}" for el, v in amounts.items())


def structure_payload(structure):
    """Atoms, bonds and cell of a structure for the 3D viewer."""
    atoms = []
    for species, label, xyz in zip(structure.species, structure.labels, structure.coords):
        parts = sorted(species.items(), key=lambda item: -item[1])
        el = parts[0][0]                        # majority element
        atoms.append({
            "el": el,
            "label": label or el,
            "species": ", ".join(f"{sp} {occ:.3g}" for sp, occ in species.items()),
            # occupancy per element, largest first, for drawing partly occupied sites
            "parts": [[sp, round(occ, 4), color(sp)] for sp, occ in parts],
            "occ": round(sum(species.values()), 4),
            "xyz": [round(float(x), 5) for x in xyz],
            "r": COVALENT_RADII.get(el, 1.5),
            "color": color(el),
        })

    # Directed bonds i -> j (+ cell offset); every bond appears in both directions.
    # Covalent radii overestimate the contact distances of large cations, so in
    # structures with nonmetals, metal-metal contacts (e.g. Na-Na in NaF, K-K in
    # KCl) are not drawn as bonds; pure metals and alloys keep them.
    radii = np.array([a["r"] for a in atoms])
    metal = np.array([is_metal(a["el"]) for a in atoms])
    cutoff = 2 * radii.max() + BOND_TOLERANCE
    centers, neighbors, offsets, distances = structure.neighbor_list(cutoff)
    bonded = ((distances < radii[centers] + radii[neighbors] + BOND_TOLERANCE)
              & (distances > MIN_BOND))
    if not metal.all():
        bonded &= ~(metal[centers] & metal[neighbors])
    bonds = [[int(i), int(j), [int(n) for n in offset]]
             for i, j, offset in zip(centers[bonded], neighbors[bonded], offsets[bonded])]

    return {
        "atoms": atoms,
        "bonds": bonds,
        "L": structure.lattice.round(6).tolist(),
        "abc": [round(float(x), 4) for x in structure.abc],
        "angles": [round(float(x), 3) for x in structure.angles],
        "volume": round(float(structure.volume), 3),
        "formula": formula(structure),
        "space_group": structure.space_group,
        "ordered": structure.is_ordered,
    }


def pdf_payload(structure, r_max, dr):
    """R(r) histograms of all element pairs, with the X-ray pair weights and <b_coh>^2 needed to
    convert them into the other functions. For neutrons, the browser calculates the weights
    from the scattering lengths and atomic fractions, as b of H depends on the deuteration
    chosen there (b_H = b(natural H) at d_natural % D, linear up to b_D at 100 %)."""
    if not 0 < r_max <= MAX_R:
        raise ClientError(f"r max must be between 0 and {MAX_R:g} Å.")
    if not 0 < dr or (r_max + R_EXTRA) / dr > MAX_BINS:
        raise ClientError(f"The bin width must be positive and give at most {MAX_BINS} bins.")
    start = time.time()
    table, rho = gen_R_r(structure, r_max, dr, verbose=False, r_extra=R_EXTRA)
    weights, mean_b2 = pair_weights(structure)
    # neutrons: the browser converts the X-ray partials with the ratio of the bb factors
    composition = structure.composition()
    total = sum(composition.values())
    try:
        neutron = {"b": {el: scattering_length(el, "neutron") for el in composition},
                   "fraction": {el: amount / total for el, amount in composition.items()},
                   "b_D": NEUTRON_B["D"], "d_natural": D_NATURAL}
    except ValueError as e:
        neutron = {"error": str(e)}
    return {
        "r_max": r_max,
        "dr": dr,
        "r0": float(table["radii"][0]),
        "n": len(table["radii"]),
        "rho": rho,
        "n_atoms": structure.num_atoms(),
        "mean_b2": mean_b2,
        "weights": weights,
        "neutron": neutron,
        "partials": {label: [float(f"{y:.7g}") for y in table[label]] for label in partial_labels(table)},
        "seconds": round(time.time() - start, 2),
    }


def distance_types(structure, r_max):
    """Groups all atom pairs up to r_max into distance types for highlighting them in the viewer.

    Pairs of the same two site labels at the same distance (to 0.001 Å) are
    symmetry equivalent and form one type, represented by one of its pairs.
    Within every window of TYPE_WINDOW only the TYPES_PER_WINDOW types with
    the largest X-ray weight and those with the largest neutron weight are kept.

    Args:
        structure (Structure): Structure as shown in the viewer.
        r_max (float): Maximum distance.

    Returns:
        dictionary: Per type the distance "d", the representative pair "i", "j"
        (atom indices of the viewer) and "img" (cell offset of j), the number
        of pairs per cell "n" and the weight "w" (sum over the pairs of the
        occupancies times atomic numbers of both atoms). "norm" converts a weight
        into the area of the type in R(r): area = norm * w, with
        norm = 2 / (N <Z>^2) (both directions of every pair, per atom).
        With neutrons, the weight (sum of the signed b_i b_j) depends on the
        scattering length b_H of H, which the deuteration sets in the browser:
        w_n = wn0 + wn1 b_H + wn2 b_H^2; its norm follows from <b> there.
    """
    _, label_index = np.unique(np.array(structure.labels), return_inverse=True)
    weight = structure.site_weights(ATOMIC_NUMBERS)
    # neutron scattering length of every site: u + v b_H, with v the occupancy of H
    try:
        u = np.array([sum(NEUTRON_B[el] * occ for el, occ in sp.items() if el != "H") for sp in structure.species])
        v = np.array([sp.get("H", 0.0) for sp in structure.species])
    except KeyError:
        u = v = None
    centers, neighbors, offsets, distances = structure.neighbor_list(r_max)
    # every pair once (the neighbor list contains both directions)
    image_key = offsets[:, 0] * 10**6 + offsets[:, 1] * 10**3 + offsets[:, 2]
    once = (centers < neighbors) | ((centers == neighbors) & (image_key > 0))
    c, n, off, d = centers[once], neighbors[once], offsets[once], distances[once]
    w = weight[c] * weight[n]

    # group by (label pair, distance to 0.001 Å); np.unique gives the first pair of every group
    la, lb = np.minimum(label_index[c], label_index[n]), np.maximum(label_index[c], label_index[n])
    keys = np.column_stack([la, lb, np.round(d * 1000).astype(np.int64)])
    _, first, inverse, count = np.unique(keys, axis=0, return_index=True, return_inverse=True, return_counts=True)
    inverse = inverse.ravel()
    w_sum = np.bincount(inverse, weights=w)
    if u is not None:
        wn = [np.bincount(inverse, weights=x) for x in (u[c] * u[n], u[c] * v[n] + v[c] * u[n], v[c] * v[n])]
    else:
        wn = [np.zeros(len(first))] * 3
    dk = np.round(d[first], 3)

    # keep the strongest types of every window, by X-ray and by neutron weight (with neutrons,
    # pairs with H can be strong although their X-ray weight is small), with neutrons for
    # hydrogen and for deuterium, then sort by distance
    window = np.floor(dk / TYPE_WINDOW).astype(np.int64)

    def strongest(weights):
        order = np.lexsort((-weights, window))
        starts = np.r_[0, np.nonzero(np.diff(window[order]))[0] + 1]
        rank = np.arange(len(order)) - np.repeat(starts, np.diff(np.r_[starts, len(order)]))
        return order[rank < TYPES_PER_WINDOW]

    wn_at = lambda b_h: wn[0] + wn[1] * b_h + wn[2] * b_h**2
    keep = strongest(w_sum)
    for b_h in (NEUTRON_B["H"], NEUTRON_B["D"]):
        keep = np.union1d(keep, strongest(np.abs(wn_at(b_h))))
    keep = keep[np.argsort(dk[keep], kind="stable")]
    rows = first[keep]

    n_atoms = structure.num_atoms()
    mean_z = weight.sum() / n_atoms
    return {
        "norm": 2 / (n_atoms * mean_z**2),
        # neutrons: coefficients of the signed weights in b_H (None without scattering lengths)
        **{f"wn{k}": [float(f"{x:.6g}") for x in wn[k][keep]] if u is not None else None for k in range(3)},
        "d": dk[keep].tolist(), "n": count[keep].tolist(),
        "w": [float(f"{x:.4g}") for x in w_sum[keep]],
        "i": c[rows].tolist(), "j": n[rows].tolist(), "img": off[rows].tolist(),
    }


def _settings(body):
    try:
        return float(body.get("r_max", 20)), float(body.get("dr", 0.01))
    except (TypeError, ValueError):
        raise ClientError("r max and bin width must be numbers.")


def _block(body, structures):
    """Index of the requested structure of a CIF file."""
    try:
        block = int(body.get("block", 0))
    except (TypeError, ValueError):
        block = -1
    if not 0 <= block < len(structures):
        raise ClientError("Unknown structure.")
    return block


def api_load(body):
    """Reads a CIF file and calculates the PDF of every structure in it."""
    structures, messages = read_cif(str(body.get("text", "")))
    r_max, dr = _settings(body)
    return {
        "warnings": messages,
        "structures": [{"block": k, "name": name,
                        "structure": structure_payload(s),
                        "pdf": pdf_payload(s, r_max, dr)}
                       for k, (name, s) in enumerate(structures)],
    }


def api_pdf(body):
    """Recalculates the PDF of one structure of a CIF file with new settings."""
    structures, _ = read_cif(str(body.get("text", "")))
    block = _block(body, structures)
    r_max, dr = _settings(body)
    return {"pdf": pdf_payload(structures[block][1], r_max, dr)}


def api_distances(body):
    """Distance types of one structure of a CIF file, for highlighting distances in the viewer."""
    structures, _ = read_cif(str(body.get("text", "")))
    block = _block(body, structures)
    r_max, _ = _settings(body)
    if not 0 < r_max <= MAX_R:
        raise ClientError(f"r max must be between 0 and {MAX_R:g} Å.")
    # the broadened contributions of distances up to R_EXTRA beyond r_max reach into the plot
    return {"distances": distance_types(structures[block][1], r_max + R_EXTRA)}


API = {"/api/load": api_load, "/api/pdf": api_pdf, "/api/distances": api_distances}

# file and settings that the GUI opens when it starts (set by the command line tool)
PRELOAD = {}

# Host names under which the server answers (set by serve()); empty: any. A server on the loopback
# interface only answers to loopback names, so a web page cannot reach it through a domain of its
# own that resolves to 127.0.0.1 (DNS rebinding) and read, e.g., the preloaded file.
LOOPBACK = {"localhost", "127.0.0.1", "::1"}
ALLOWED_HOSTS = set()


def _hostname(netloc):
    """'localhost:8766' -> 'localhost', '[::1]:8766' -> '::1'."""
    try:
        return (urlsplit(f"//{netloc}").hostname or "").lower()
    except ValueError:
        return ""


def _forbidden(environ, path):
    """Reason to refuse a request from outside the GUI, or None."""
    host = environ.get("HTTP_HOST", "")
    if ALLOWED_HOSTS and _hostname(host) not in ALLOWED_HOSTS:
        return "unknown host"
    # the calculations can only be started by the GUI itself, not by other web pages
    origin = environ.get("HTTP_ORIGIN")
    if path in API and origin is not None and urlsplit(origin).netloc.lower() != host.lower():
        return "cross-origin request"
    return None


def app(environ, start_response):
    """WSGI application serving the API and the static files."""
    method, path = environ["REQUEST_METHOD"], environ.get("PATH_INFO", "/")

    def reply(status, body, content_type="application/json"):
        start_response(status, [("Content-Type", content_type),
                                ("Content-Length", str(len(body))),
                                ("Cache-Control", "no-store")])
        return [body]

    def js(obj, status="200 OK"):
        return reply(status, json.dumps(obj).encode())

    reason = _forbidden(environ, path)
    if reason:
        return reply("403 Forbidden", reason.encode(), "text/plain")
    if path == "/api/startup":
        return js(PRELOAD)
    if path in API:
        if method != "POST":
            return js({"error": "POST required"}, "405 Method Not Allowed")
        try:
            length = int(environ.get("CONTENT_LENGTH") or 0)
            if length > MAX_BODY:
                raise ClientError("The file is too large.")
            body = json.loads(environ["wsgi.input"].read(length) or b"{}")
            return js(API[path](body))
        except ClientError as e:
            return js({"error": str(e)}, "400 Bad Request")
        except Exception as e:
            return js({"error": f"Calculation failed: {e}"}, "500 Internal Server Error")

    # static files; the resolved path must stay inside the static folder
    rel = "index.html" if path in ("", "/") else path.lstrip("/")
    full = os.path.realpath(os.path.join(STATIC, rel))
    if method != "GET" or not full.startswith(STATIC + os.sep) or not os.path.isfile(full):
        return reply("404 Not Found", b"not found", "text/plain")
    content_type = mimetypes.guess_type(full)[0] or "application/octet-stream"
    if content_type in ("text/javascript", "application/javascript"):
        content_type = "text/javascript; charset=utf-8"
    with open(full, "rb") as fh:
        return reply("200 OK", fh.read(), content_type)


class ThreadingWSGIServer(ThreadingMixIn, WSGIServer):
    daemon_threads = True


class QuietHandler(WSGIRequestHandler):
    def log_message(self, *args):
        pass


def serve(host="127.0.0.1", port=8766, open_browser=True, preload=None):
    """Starts the GUI on a local threaded server.

    Args:
        host (string): Address to listen on.
        port (int): Port; if it is in use, a free one is taken.
        open_browser (Boolean): If True, the GUI is opened in the default browser.
        preload (dictionary): File ("file", "text") and command line settings
        ("settings") that the GUI opens at the start.
    """
    PRELOAD.clear()
    PRELOAD.update(preload or {})
    ALLOWED_HOSTS.clear()
    if host in LOOPBACK or host.startswith("127."):
        ALLOWED_HOSTS.update(LOOPBACK | {host})
    try:
        httpd = make_server(host, port, app, server_class=ThreadingWSGIServer,
                            handler_class=QuietHandler)
    except OSError:
        httpd = make_server(host, 0, app, server_class=ThreadingWSGIServer,
                            handler_class=QuietHandler)
    url = f"http://{host}:{httpd.server_port}/"
    print(f"xPDFsim GUI running at {url}  (Ctrl+C to stop)")
    if open_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass


def main():
    """Command line interface of the GUI."""
    parser = argparse.ArgumentParser(description="xPDFsim graphical user interface")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8766)
    parser.add_argument("--no-browser", action="store_true",
                        help="do not open the GUI in the default browser")
    args = parser.parse_args()
    serve(args.host, args.port, not args.no_browser)


if __name__ == "__main__":
    main()
