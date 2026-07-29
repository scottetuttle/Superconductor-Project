from shs.boundaries import (
    BoundaryCondition,
    BoundarySide,
    BoundaryType,
    BoundarySet,
)


def test_boundary_creation():

    boundary = BoundaryCondition(
        side=BoundarySide.LEFT,
        type=BoundaryType.FIXED_TEMPERATURE,
        temperature=4.2,
    )

    assert boundary.side == BoundarySide.LEFT

    assert boundary.type == BoundaryType.FIXED_TEMPERATURE

    assert boundary.temperature == 4.2


def test_boundary_set():

    left = BoundaryCondition(
        side=BoundarySide.LEFT,
        type=BoundaryType.FIXED_TEMPERATURE,
        temperature=4.2,
    )

    right = BoundaryCondition(
        side=BoundarySide.RIGHT,
        type=BoundaryType.FIXED_TEMPERATURE,
        temperature=4.2,
    )

    top = BoundaryCondition(
        side=BoundarySide.TOP,
        type=BoundaryType.INSULATING,
    )

    bottom = BoundaryCondition(
        side=BoundarySide.BOTTOM,
        type=BoundaryType.INSULATING,
    )

    boundaries = BoundarySet(
        boundaries={
            BoundarySide.LEFT: left,
            BoundarySide.RIGHT: right,
            BoundarySide.TOP: top,
            BoundarySide.BOTTOM: bottom,
        }
    )

    assert boundaries.get(
        BoundarySide.LEFT
    ).temperature == 4.2

    assert boundaries.get(
        BoundarySide.TOP
    ).type == BoundaryType.INSULATING

def test_boundary_add():

    boundaries = BoundarySet()

    left = BoundaryCondition(
        side=BoundarySide.LEFT,
        type=BoundaryType.FIXED_TEMPERATURE,
        temperature=4.2,
    )

    boundaries.add(left)

    assert BoundarySide.LEFT in boundaries

    assert boundaries.get(
        BoundarySide.LEFT
    ).temperature == 4.2