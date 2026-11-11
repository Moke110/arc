"""Frozen integer gene identities at the native-source/display boundary."""
import csv
import hashlib
from numbers import Integral
from collections import defaultdict
from pathlib import Path

SOURCE_NAMES = {"temp_env": "tenv", "screen_ctox": "sctox"}


AXIS_SHA256 = '19ce44acfd615a601d5bc6617074941e13d51cd6e527c8d519174e5085385f42'


class IdentityError(ValueError):
    """Missing, ambiguous or malformed gene identity."""


def read_rows(path, delimiter=','):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


class GeneIdentity:
    """Resolve identities without reordering the supplied frozen axis.

    Native identifiers require a source. Global identifiers are HGNC IDs,
    canonical symbols/Ensembl/Entrez IDs and annotation aliases. Every candidate
    is retained: collisions fail instead of choosing an arbitrary gene.
    """
    def __init__(self, axis, crosswalk=(), annotation=()):
        self.axis = [dict(row, gene_index=index) for index, row in enumerate(axis)]
        self._global = defaultdict(set)
        self._source = defaultdict(set)
        self._native = defaultdict(set)
        self.evidence = []
        by_hgnc = {}
        if not self.axis:
            raise IdentityError('empty gene axis')
        for field in ('symbol', 'hgnc_id'):
            values = [row.get(field, '') for row in self.axis]
            if not all(values) or len(set(values)) != len(values):
                raise IdentityError(f'missing or duplicate axis {field}')
        for row in self.axis:
            hgnc = row['hgnc_id']
            if not hgnc or hgnc in by_hgnc:
                raise IdentityError('missing or duplicate HGNC axis identity')
            by_hgnc[hgnc] = row['gene_index']
            for field in ('symbol', 'hgnc_id', 'ensembl_id', 'entrez_id'):
                value = row.get(field, '')
                if value:
                    self._global[value].add(row['gene_index'])
        for row in crosswalk:
            hgnc = row.get('hgnc_id', '')
            if not hgnc:
                continue
            index = by_hgnc.get(hgnc, hgnc)
            source, identifier = row['dataset_id'], row['source_gene_id']
            source = SOURCE_NAMES.get(source, source)
            if not source or not identifier:
                raise IdentityError('missing source gene identity')
            self._source[source, identifier].add(index)
            if isinstance(index, int):
                self._native[index, source].add(identifier)
                self.evidence.append(dict(row, gene_index=index))
        # Keep off-axis annotation candidates as well, so an alias shared with
        # an off-axis gene cannot silently resolve to an on-axis gene.
        for row in annotation:
            hgnc = row['hgnc_id']
            candidate = by_hgnc.get(hgnc, hgnc)
            for field in ('alias_symbol', 'prev_symbol'):
                for alias in row.get(field, '').split('|'):
                    if alias:
                        self._global[alias].add(candidate)

    @classmethod
    def load(cls, axis, crosswalk, annotation=None, expected_size=18164, expected_sha256=AXIS_SHA256, native_mappings=()):
        if expected_sha256 and hashlib.sha256(Path(axis).read_bytes()).hexdigest() != expected_sha256:
            raise IdentityError('frozen gene axis digest mismatch')
        rows = read_rows(axis)
        if len(rows) != expected_size:
            raise IdentityError(f'expected {expected_size} frozen genes, got {len(rows)}')
        return cls(rows, [*read_rows(crosswalk), *native_mappings],
                   read_rows(annotation, '\t') if annotation else ())

    def indices(self, identifiers, source=None):
        result = []
        for identifier in identifiers:
            identifier = str(identifier)
            candidates = (self._source.get((source, identifier), set()) if source is not None
                          else self._global.get(identifier, set()))
            if len(candidates) != 1:
                reason = 'ambiguous' if candidates else 'unmapped'
                raise IdentityError(f'{reason} gene identity: {source}:{identifier}')
            value = next(iter(candidates))
            if not isinstance(value, int):
                raise IdentityError(f'off-axis gene identity: {identifier}')
            result.append(value)
        return result

    def labels(self, indices, field='symbol'):
        result = []
        for index in indices:
            if isinstance(index, bool) or not isinstance(index, Integral) or not 0 <= index < len(self.axis):
                raise IdentityError(f'invalid gene index: {index}')
            result.append(self.axis[index][field])
        return result

    def native_ids(self, index, source):
        self.labels([index])
        values = sorted(self._native.get((index, source), ()))
        if not values:
            raise IdentityError(f'no native identity for {source}:{index}')
        return values


def native_mapping_rows(shared):
    """Expose every native namespace from the independently replayed shared release."""
    mappings = []
    for source in ('tahoe', 'tenv'):
        for row in read_rows(Path(shared) / f'{source}_mapping.csv'):
            for field in ('native_symbol', 'native_ensembl_id', 'native_token_id'):
                if row[field]:
                    mappings.append({'dataset_id': source, 'source_gene_id': row[field],
                                     'symbol': row['symbol'], 'hgnc_id': row['hgnc_id'],
                                     'match_type': row['resolution'],
                                     'evidence': f'{source}_mapping.csv:{field}'})
    return mappings
