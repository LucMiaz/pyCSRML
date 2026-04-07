import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--recompute-clinventory",
        action="store_true",
        default=False,
        help=(
            "Force recomputation of CLinventory concordance matrices, "
            "ignoring any existing .npz cache files."
        ),
    )
