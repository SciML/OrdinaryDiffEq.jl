#!/usr/bin/env python3
"""
Compute which sublibrary tests need to run based on changed files.

Reads the dependency graph from lib/*/Project.toml [deps] sections,
identifies internal dependencies by matching dep names against known
sublibrary directory names, computes the transitive reverse-dependency
map, then given a list of changed files outputs affected sublibraries.

Usage:
    # Pass changed files on stdin (one per line)
    git diff --name-only origin/master...HEAD | python3 compute_affected_sublibraries.py /path/to/repo

    # Output: JSON array of affected sublibrary names, e.g.
    # ["OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", ...]
"""

import json
import os
import re
import sys
from collections import defaultdict


def parse_deps_from_toml(toml_path, known_sublibs):
    """Parse [deps] and [extras] sections from a Project.toml to find internal dependencies.

    A dependency is 'internal' if its name matches a known sublibrary directory.
    We include [extras] because test-only dependencies (used by Pkg.test()) can
    also break when upstream sublibraries change.
    Returns a list of sublibrary names this package depends on.
    """
    deps = set()
    in_relevant_section = False
    with open(toml_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped in ("[deps]", "[extras]"):
                in_relevant_section = True
                continue
            if in_relevant_section:
                if stripped.startswith("["):
                    in_relevant_section = False
                    continue
                # Match lines like: OrdinaryDiffEqCore = "uuid-here"
                m = re.match(r'^(\w+)\s*=', stripped)
                if m:
                    dep_name = m.group(1)
                    if dep_name in known_sublibs:
                        deps.add(dep_name)
    return list(deps)


def build_dependency_graph(lib_dir):
    """Build forward dependency graph: pkg -> [list of internal deps]."""
    # First, collect all sublibrary names
    known_sublibs = set()
    for entry in sorted(os.listdir(lib_dir)):
        if os.path.isfile(os.path.join(lib_dir, entry, "Project.toml")):
            known_sublibs.add(entry)

    # Then parse each sublibrary's deps
    graph = {}
    for entry in sorted(known_sublibs):
        toml = os.path.join(lib_dir, entry, "Project.toml")
        graph[entry] = parse_deps_from_toml(toml, known_sublibs)

    return graph


def compute_reverse_deps(graph):
    """Build reverse dependency map: pkg -> set of packages that depend on it (transitively)."""
    # Direct reverse deps
    rev = defaultdict(set)
    for pkg, deps in graph.items():
        for dep in deps:
            rev[dep].add(pkg)

    # Compute transitive closure
    def get_all_rdeps(pkg, visited=None):
        if visited is None:
            visited = set()
        if pkg in visited:
            return visited
        visited.add(pkg)
        for rdep in rev.get(pkg, set()):
            get_all_rdeps(rdep, visited)
        return visited

    transitive = {}
    for pkg in graph:
        rdeps = get_all_rdeps(pkg)
        rdeps.discard(pkg)  # don't include self
        transitive[pkg] = rdeps

    return transitive


def compute_affected(changed_files, graph, reverse_deps):
    """Given changed files, return set of sublibraries that need testing."""
    affected = set()

    for filepath in changed_files:
        filepath = filepath.strip()
        if not filepath:
            continue

        parts = filepath.split("/")

        # Check if the change is inside lib/<sublibrary>/
        if len(parts) >= 2 and parts[0] == "lib" and parts[1] in graph:
            pkg = parts[1]
            # This sublibrary itself is affected
            affected.add(pkg)
            # All transitive reverse dependencies are also affected
            affected.update(reverse_deps.get(pkg, set()))

    return affected


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <repo_root>", file=sys.stderr)
        sys.exit(1)

    repo_root = sys.argv[1]
    lib_dir = os.path.join(repo_root, "lib")

    if not os.path.isdir(lib_dir):
        print(f"Error: {lib_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    graph = build_dependency_graph(lib_dir)
    reverse_deps = compute_reverse_deps(graph)

    changed_files = sys.stdin.read().strip().split("\n")
    affected = compute_affected(changed_files, graph, reverse_deps)

    # Sort for deterministic output
    result = sorted(affected)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
