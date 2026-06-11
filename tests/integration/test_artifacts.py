"""Tier B — load each produced artifact and assert shape.

Initial smoke version: one bus-count assertion to prove the fixture
plumbing works in CI. Later PRs flesh this out with per-stage assertions.
"""

import pypsa
import pytest


@pytest.mark.integration
def test_post_cluster_simpl_bus_count(built):
    n = pypsa.Network(str(built.elec_s))
    assert len(n.buses) == 20, f"cluster_simpl produced {len(n.buses)} buses, expected 20 (config.test.yaml simpl=20)"
