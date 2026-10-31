"""Reading structures from CIF files, without dependencies besides numpy.

The tokenizer and the parser of symmetry operations are adapted from
NetsUnhinged (same author). Text fields (';' ... ';') are single tokens, so
files with embedded data, e.g. the refinement fcf file that Olex2 stores in
_iucr_refine_fcf_details, are read correctly.
"""

import re
import numpy as np

from .elements import ATOMIC_NUMBERS
from .structure import Structure, cell_matrix

MERGE_TOL = 0.01                # Å: positions closer than this are the same site
_NEIGHBOR_BINS = np.array(np.meshgrid(*[[-1, 0, 1]] * 3, indexing="ij")).reshape(3, -1).T


# ---------------------------------------------------------------- tokenizer

class _Quoted(str):
    """A quoted string or text field: a value even if it reads like 'loop_', 'data_x' or '_tag'."""


def _tokenize(text):
    """Yields CIF tokens; text fields (;...;) and quoted strings become one _Quoted token."""
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(";"):
            buf = [line[1:]]
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            i += 1
            yield _Quoted(("\n".join(buf)).strip())
            continue
        j, n = 0, len(line)
        while j < n:
            ch = line[j]
            if ch.isspace():
                j += 1
            elif ch == "#":
                break
            elif ch in "'\"":
                k = j + 1
                # a closing quote must be followed by whitespace or the end of the line
                while k < n and not (line[k] == ch and (k + 1 == n or line[k + 1].isspace())):
                    k += 1
                yield _Quoted(line[j + 1:k])
                j = k + 1
            else:
                k = j
                while k < n and not line[k].isspace():
                    k += 1
                yield line[j:k]
                j = k
        i += 1


def _parse_blocks(text):
    """Returns a list of (block name, items, loops, warnings); items: tag -> value,
    loops: (tags, rows)."""
    blocks, cur = [], None
    tokens = list(_tokenize(text))
    i = 0
    while i < len(tokens):
        t = tokens[i]
        tl = t.lower()
        if isinstance(t, _Quoted):
            i += 1                                      # a value without a tag
        elif tl.startswith("data_"):
            cur = (t[5:], {}, [], [])
            blocks.append(cur)
            i += 1
        elif cur is None:
            i += 1
        elif tl == "loop_":
            i += 1
            tags = []
            while i < len(tokens) and _is_tag(tokens[i]):
                tags.append(tokens[i].lower())
                i += 1
            values = []
            while i < len(tokens):
                v = tokens[i]
                vl = v.lower()
                if not isinstance(v, _Quoted) and (v.startswith("_") or vl == "loop_" or vl.startswith("data_")):
                    break
                values.append(v)
                i += 1
            if tags:
                nrow, rest = divmod(len(values), len(tags))
                if rest:
                    # e.g. an unquoted value with a space: the columns after it are shifted
                    cur[3].append(f"The loop of {tags[0]} has {len(values)} values, which is not a "
                                  f"multiple of its {len(tags)} tags; the last {rest} were ignored. "
                                  "Values with spaces must be quoted.")
                cur[2].append((tags, [values[r * len(tags):(r + 1) * len(tags)] for r in range(nrow)]))
        elif t.startswith("_"):
            if i + 1 < len(tokens):
                cur[1][tl] = tokens[i + 1]
            i += 2
        else:
            i += 1
    return blocks


def _is_tag(token):
    return token.startswith("_") and not isinstance(token, _Quoted)


def _num(s):
    """Parses a CIF number, dropping the standard uncertainty '(..)'."""
    s = s.strip()
    if s in (".", "?", ""):
        raise ValueError("missing number")
    return float(re.sub(r"\(.*?\)", "", s))


# ------------------------------------------------------------ symmetry operations

_SYMOP_TOKEN = re.compile(r"\s*(?:(\d+\.\d*|\.\d+|\d+)|([xyz])|([-+*/()]))", re.I)


def _parse_affine(expr):
    """Parses one component of a symmetry operation ('-x+1/2', 'x-y', '1/2+z', ...) into the
    coefficients of x, y, z and a constant, with a small recursive-descent parser: no eval,
    bounded length and nesting, so a crafted file can neither run code nor hang the parser."""
    if len(expr) > 40:
        raise ValueError("expression too long")
    toks, pos = [], 0
    while pos < len(expr):
        m = _SYMOP_TOKEN.match(expr, pos)
        if not m or m.end() == pos:
            raise ValueError("unexpected character")
        num, var, op = m.groups()
        if num is not None:
            toks.append(("n", float(num)))
        elif var is not None:
            toks.append(("v", "xyz".index(var.lower())))
        else:
            toks.append(("o", op))
        pos = m.end()
        if expr[pos:].strip() == "":
            break
    k = [0]

    def peek():
        return toks[k[0]] if k[0] < len(toks) else (None, None)

    def take():
        tok = peek()
        k[0] += 1
        return tok

    def factor(depth):
        kind, val = take()
        if kind == "n":
            out = (np.zeros(3), val)
            if peek()[0] == "v":                      # implicit product such as '2x'
                _, idx = take()
                vec = np.zeros(3)
                vec[idx] = val
                out = (vec, 0.0)
            return out
        if kind == "v":
            vec = np.zeros(3)
            vec[val] = 1.0
            return vec, 0.0
        if kind == "o" and val in "+-":
            vec, c = factor(depth)
            return (-vec, -c) if val == "-" else (vec, c)
        if kind == "o" and val == "(" and depth < 3:
            out = expr_(depth + 1)
            if take() != ("o", ")"):
                raise ValueError("unbalanced parentheses")
            return out
        raise ValueError("syntax error")

    def term(depth):
        vec, c = factor(depth)
        while peek() in (("o", "*"), ("o", "/")):
            _, op = take()
            v2, c2 = factor(depth)
            if op == "*":
                if vec.any() and v2.any():
                    raise ValueError("product of variables")
                vec, c = vec * c2 + v2 * c, c * c2
            else:
                if v2.any() or c2 == 0:
                    raise ValueError("division by a variable or zero")
                vec, c = vec / c2, c / c2
        return vec, c

    def expr_(depth):
        vec, c = term(depth)
        while peek() in (("o", "+"), ("o", "-")):
            _, op = take()
            v2, c2 = term(depth)
            vec, c = (vec + v2, c + c2) if op == "+" else (vec - v2, c - c2)
        return vec, c

    vec, c = expr_(0)
    if k[0] != len(toks):
        raise ValueError("trailing characters")
    if np.abs(vec).max() > 3 or abs(c) > 3:
        raise ValueError("coefficient out of range")
    return vec, c


def parse_symop(op):
    """'-x+1/2, y, z+1/2' -> (rotation 3x3, translation 3)."""
    parts = [p.strip() for p in op.split(",")]
    if len(parts) != 3:
        raise ValueError(f"bad symmetry operation '{op[:60]}'")
    R, t = np.zeros((3, 3)), np.zeros(3)
    for row, expr in enumerate(parts):
        try:
            R[row], t[row] = _parse_affine(expr)
        except ValueError as exc:
            raise ValueError(f"bad symmetry operation '{op[:60]}' ({exc})") from None
    return R, t


# ---------------------------------------------------------------- reader

_ELEMENT_RE = re.compile(r"^([A-Z][a-z]?)")


def _element(type_symbol, label):
    """Element symbol from _atom_site_type_symbol (e.g. 'Na', 'Na+', 'O2-') or the label."""
    for src in (type_symbol, label):
        if not src or src in (".", "?"):
            continue
        s = src.strip()
        s = s[0].upper() + s[1:].lower() if len(s) > 1 else s.upper()
        m = _ELEMENT_RE.match(s)
        if m:
            if m.group(1) in ATOMIC_NUMBERS:
                return m.group(1)
            if m.group(1)[0] in ATOMIC_NUMBERS:
                return m.group(1)[0]
    raise ValueError(f"cannot determine the element of atom '{str(label)[:32]}'")


def _formula_amounts(formula):
    """'C5 H Cu N4' -> {'C': 5, 'H': 1, 'Cu': 1, 'N': 4}."""
    amounts = {}
    for el, n in re.findall(r"([A-Z][a-z]?)\s*(\d+\.?\d*|\.\d+)?", formula):
        amounts[el] = amounts.get(el, 0.0) + (float(n) if n else 1.0)
    return amounts


def _symmetry_operations(items, loops):
    for tags, rows in loops:
        for key in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz",
                    "_space_group_symop.operation_xyz"):
            if key in tags:
                col = tags.index(key)
                ops = [r[col] for r in rows if r[col].strip() not in ("?", ".", "")]
                if ops:
                    return [parse_symop(op) for op in ops], []
    for key in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"):
        if items.get(key, "?") not in ("?", "."):
            return [parse_symop(items[key])], []
    hm = _space_group(items)
    if hm and hm.replace(" ", "") != "P1":
        raise ValueError(f"The CIF lists no symmetry operations for space group {hm}. Please export "
                         "the CIF with symmetry operations (e.g. from Mercury or VESTA).")
    return [parse_symop("x,y,z")], ([] if hm else ["No symmetry operations found - assuming P1."])


def _space_group(items):
    for key in ("_space_group_name_h-m_alt", "_symmetry_space_group_name_h-m"):
        value = items.get(key, "?")
        if value not in ("?", "."):
            return value.strip()
    return None


def _read_block(name, items, loops):
    """One structure from a data block with atom sites (see read_structures)."""
    try:
        cell = [_num(items[k]) for k in ("_cell_length_a", "_cell_length_b", "_cell_length_c",
                                         "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")]
    except KeyError as exc:
        raise ValueError(f"missing cell parameter {exc}") from None
    L = cell_matrix(*cell)
    ops, warnings = _symmetry_operations(items, loops)

    tags, rows = next((tg, r) for tg, r in loops if "_atom_site_fract_x" in tg)
    col = lambda tag: tags.index(tag) if tag in tags else None
    c_label, c_type = col("_atom_site_label"), col("_atom_site_type_symbol")
    c_xyz = [col(f"_atom_site_fract_{k}") for k in "xyz"]
    c_occ, c_calc = col("_atom_site_occupancy"), col("_atom_site_calc_flag")

    # asymmetric unit; every atom is expanded by symmetry, its copies on special positions are
    # merged; atoms of different labels on the same position form one site whose occupancies
    # are the sums per element (substitutional disorder, split sites)
    sites = []                  # [frac, {element: occupancy}, label]
    # sites are found by their bin in a grid whose bins are at least MERGE_TOL thick, so a
    # close site lies in the same or a neighbouring bin (periodic)
    spacing = 1 / np.linalg.norm(np.linalg.inv(L).T, axis=1)
    grid = np.maximum(np.floor(spacing / MERGE_TOL), 1).astype(int)
    by_bin = {}                 # bin -> indices of sites

    def nearby(f):
        b = np.minimum(np.floor(f * grid).astype(int), grid - 1)
        keys = {tuple((b + d) % grid) for d in _NEIGHBOR_BINS}
        return tuple(b), [k for key in keys for k in by_bin.get(key, ())]

    for r in rows:
        label = r[c_label] if c_label is not None else (r[c_type] if c_type is not None else "X")
        if c_calc is not None and r[c_calc].lower() == "dum":
            continue
        try:
            frac = np.array([_num(r[c]) for c in c_xyz])
        except (ValueError, TypeError):
            continue
        occ = 1.0
        if c_occ is not None:
            try:
                occ = _num(r[c_occ])
            except ValueError:
                pass
        el = _element(r[c_type] if c_type is not None else None, label)
        images = np.array([(R @ frac + t) % 1.0 for R, t in ops])
        diff = images[:, np.newaxis, :] - images[np.newaxis, :, :]
        close = np.linalg.norm((diff - np.round(diff)) @ L, axis=2) < MERGE_TOL
        unique = []                             # copies on special positions coincide
        for k in range(len(images)):
            if not close[k, unique].any():
                unique.append(k)
        for f in images[unique]:
            key, candidates = nearby(f)
            if candidates:
                d = np.array([sites[k][0] for k in candidates]) - f
                dist = np.linalg.norm((d - np.round(d)) @ L, axis=1)
                best = int(np.argmin(dist))
                if dist[best] < MERGE_TOL:
                    site = sites[candidates[best]]
                    site[1][el] = site[1].get(el, 0.0) + occ
                    continue
            by_bin.setdefault(key, []).append(len(sites))
            sites.append([f, {el: occ}, label])
    if not sites:
        raise ValueError("no atom sites")

    # a site cannot be more than fully occupied (rounded occupancies of split sites)
    n_scaled = 0
    for s in sites:
        total = sum(s[1].values())
        if total > 1.0001:
            s[1] = {el: occ / total for el, occ in s[1].items()}
            n_scaled += 1
    if n_scaled:
        warnings.append(f"{n_scaled} sites with a total occupancy above 1 were scaled to 1.")

    structure = Structure(L, [s[0] for s in sites], [s[1] for s in sites], [s[2] for s in sites],
                          name=name, space_group=_space_group(items))

    # compare with the formula given in the file (catches missing H atoms and the like)
    if items.get("_chemical_formula_sum", "?") not in ("?", "."):
        given, found = _formula_amounts(items["_chemical_formula_sum"]), structure.composition()
        ref = next((el for el in given if el in found and given[el] > 0), None)
        if ref:
            scale = found[ref] / given[ref]
            missing = [el for el in given if el not in found]
            off = [el for el in given if el in found and abs(found[el] / scale - given[el]) > 0.02 * given[el] + 0.01]
            if missing:
                warnings.append(f"Elements of the formula in the CIF ({items['_chemical_formula_sum']}) "
                                f"are missing in the atom sites: {', '.join(missing)}.")
            if off:
                warnings.append(f"The atom sites do not match the formula in the CIF "
                                f"({items['_chemical_formula_sum']}) for {', '.join(off)}.")
    return structure, warnings


def read_structures(text):
    """Reads every structure (data block with atom sites) of a CIF file.

    Args:
        text (string): Content of the CIF file.

    Returns:
        list: (block name, Structure) of every readable block.
        list: Warnings.

    Raises:
        ValueError: If the file contains no readable structure.
    """
    structures, messages, errors = [], [], []
    for name, items, loops, block_warnings in _parse_blocks(text):
        if not any("_atom_site_fract_x" in tags for tags, _ in loops):
            continue
        try:
            structure, warnings = _read_block(name, items, loops)
        except ValueError as e:
            errors.append(f"{name}: {e}")
            continue
        structures.append((name, structure))
        messages += [w for w in block_warnings + warnings if w not in messages]
    if not structures:
        raise ValueError(errors[0] if errors else "No structure (data block with atom sites) found in the file.")
    return structures, messages
