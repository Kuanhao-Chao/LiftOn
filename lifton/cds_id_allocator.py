"""Deterministic, file-scoped allocation of emitted CDS identifiers."""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable


_COPY_SUFFIX = re.compile(r"_(\d+)$")


class CdsIdAllocator:
    """Reserve stable IDs and collision-resolve only synthesized CDS IDs.

    ``reserved_source_ids`` contains explicit source identifiers in LiftOn's
    synthetic ``cds-`` namespace.  Reserving them before the first output row
    prevents an early ID-less transcript from taking an identifier that a
    later transcript must preserve.
    """

    def __init__(
        self,
        reserved_source_ids: Iterable[str] = (),
        *,
        exact_reserved_source_ids: Iterable[str] = (),
        source_id_exists: Callable[[str], bool] | None = None,
        exact_source_id_exists: Callable[[str], bool] | None = None,
    ) -> None:
        self._reserved_source_ids = {
            str(identifier)
            for identifier in reserved_source_ids
            if identifier
        }
        self._exact_reserved_source_ids = {
            str(identifier)
            for identifier in exact_reserved_source_ids
            if identifier
        }
        self._source_id_exists = source_id_exists
        self._exact_source_id_exists = (
            exact_source_id_exists or source_id_exists
        )
        self._source_lookup_cache: dict[str, bool] = {}
        self._exact_source_lookup_cache: dict[str, bool] = {}
        self._claims: dict[str, tuple[str, str]] = {}

    @staticmethod
    def _copy_bases(identifier: str):
        """Yield an ID and every base obtained by removing copy suffixes."""

        current = identifier
        yield current
        while True:
            match = _COPY_SUFFIX.search(current)
            if match is None:
                return
            current = current[:match.start()]
            yield current

    def _is_source_reserved(self, identifier: str) -> bool:
        for candidate in self._copy_bases(identifier):
            if candidate in self._reserved_source_ids:
                return True
            if self._source_id_exists is None:
                continue
            cached = self._source_lookup_cache.get(candidate)
            if cached is None:
                try:
                    cached = bool(self._source_id_exists(candidate))
                except Exception:
                    cached = False
                self._source_lookup_cache[candidate] = cached
            if cached:
                return True
        return False

    def _is_exact_source_reserved(self, identifier: str) -> bool:
        if (
            identifier in self._reserved_source_ids
            or identifier in self._exact_reserved_source_ids
        ):
            return True
        if self._exact_source_id_exists is None:
            return False
        cached = self._exact_source_lookup_cache.get(identifier)
        if cached is None:
            try:
                cached = bool(self._exact_source_id_exists(identifier))
            except Exception:
                cached = False
            self._exact_source_lookup_cache[identifier] = cached
        return cached

    def reserve_explicit(
        self,
        identifier: str,
        parent_id: str,
        *,
        source_identifier: str | None = None,
    ) -> bool:
        """Claim a preserved ID, returning false on a cross-parent conflict."""

        identifier = str(identifier)
        parent_id = str(parent_id)
        source_identifier = (
            identifier if source_identifier is None
            else str(source_identifier)
        )
        # An original source ID is entitled to its own reserved name. A copy
        # transformation is not: ``source_id + _N`` may be another feature's
        # exact stable source ID, which must win the collision.
        if (
            identifier != source_identifier
            and self._is_exact_source_reserved(identifier)
        ):
            return False
        claim = self._claims.get(identifier)
        if claim is None:
            self._claims[identifier] = ("explicit", parent_id)
            return True
        return claim == ("explicit", parent_id)

    def allocate_synthetic(
        self,
        base: str,
        parent_id: str,
        *,
        local_reserved_ids: Iterable[str] = (),
    ) -> str:
        """Return and claim a deterministic collision-free synthetic ID."""

        base = str(base)
        parent_id = str(parent_id)
        local_reserved = {
            str(identifier) for identifier in local_reserved_ids if identifier
        }
        index = 0
        while True:
            candidate = base if index == 0 else f"{base}-{index}"
            claim = self._claims.get(candidate)
            if claim == ("synthetic", parent_id):
                return candidate
            if (
                candidate not in local_reserved
                and claim is None
                and not self._is_source_reserved(candidate)
            ):
                self._claims[candidate] = ("synthetic", parent_id)
                return candidate
            index += 1
