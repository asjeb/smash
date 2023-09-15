from __future__ import annotations

from smash._constant import (
    STRUCTURE_NAME,
    OPR_PARAMETERS,
    OPR_STATES,
    STRUCTURE_OPR_PARAMETERS,
    STRUCTURE_OPR_STATES,
    STRUCTURE_COMPUTE_CI,
    FEASIBLE_OPR_PARAMETERS,
    FEASIBLE_OPR_INITIAL_STATES,
    DEFAULT_OPR_PARAMETERS,
    DEFAULT_OPR_INITIAL_STATES,
    DEFAULT_BOUNDS_OPR_PARAMETERS,
    DEFAULT_BOUNDS_OPR_INITIAL_STATES,
    OPTIMIZABLE_OPR_PARAMETERS,
    OPTIMIZABLE_OPR_INITIAL_STATES,
    INPUT_DATA_FORMAT,
    RATIO_PET_HOURLY,
    DATASET_NAME,
    PEAK_QUANT,
    MAX_DURATION,
)

import numpy as np


def test_structure_model():
    # % Check structure name
    assert STRUCTURE_NAME == [
        "gr4-lr",
        "gr4-kw",
        "gr5-lr",
        "gr5-kw",
        "loieau-lr",
        "grd-lr",
    ]

    # % Check opr parameters
    assert OPR_PARAMETERS == [
        "ci",
        "cp",
        "ct",
        "kexc",
        "aexc",
        "ca",
        "cc",
        "kb",
        "llr",
        "akw",
        "bkw",
    ]

    # % Check opr states
    assert OPR_STATES == ["hi", "hp", "ht", "ha", "hc", "hlr"]

    # % Check structure opr parameter
    assert list(STRUCTURE_OPR_PARAMETERS.values()) == [
        ["ci", "cp", "ct", "kexc", "llr"],  # % gr4-lr
        ["ci", "cp", "ct", "kexc", "akw", "bkw"],  # % gr4-kw
        ["ci", "cp", "ct", "kexc", "aexc", "llr"],  # % gr5-lr
        ["ci", "cp", "ct", "kexc", "aexc", "akw", "bkw"],  # % gr5-kw
        ["ca", "cc", "kb", "llr"],  # % loieau-lr
        ["cp", "ct", "llr"],  # % grd-lr
    ]

    # % Check structure opr state
    assert list(STRUCTURE_OPR_STATES.values()) == [
        ["hi", "hp", "ht", "hlr"],  # % gr4-lr
        ["hi", "hp", "ht"],  # % gr4-kw
        ["hi", "hp", "ht", "hlr"],  # % gr5-lr
        ["hi", "hp", "ht"],  # % gr5-kw
        ["ha", "hc", "hlr"],  # % loieau-lr
        ["hp", "ht", "hlr"],  # % grd-lr
    ]


def test_feasible_domain():
    # % Feasible opr parameters
    assert list(FEASIBLE_OPR_PARAMETERS.values()) == [
        (0, np.inf),  # % ci
        (0, np.inf),  # % cp
        (0, np.inf),  # % ct
        (-np.inf, np.inf),  # % kexc
        (0, 1),  # % aexc
        (0, np.inf),  # % ca
        (0, np.inf),  # % cc
        (0, np.inf),  # % kb
        (0, np.inf),  # % llr
        (0, np.inf),  # % akw
        (0, np.inf),  # % bkw
    ]

    # % Feasible opr states
    assert list(FEASIBLE_OPR_INITIAL_STATES.values()) == [
        (0, 1),  # % hi
        (0, 1),  # % hp
        (0, 1),  # % ht
        (0, 1),  # % ha
        (0, 1),  # % hc
        (0, np.inf),  # % hlr
    ]


def test_default_parameters():
    # % Check default opr parameters
    assert list(DEFAULT_OPR_PARAMETERS.values()) == [
        1e-6,  # % ci
        200,  # % cp
        500,  # % ct
        0,  # % kexc
        0.1,  # % aexc
        200,  # % ca
        500,  # % cc
        1,  # % kb
        5,  # % llr
        5,  # % akw
        0.6,  # % bkw
    ]

    # % Check default opr states
    assert list(DEFAULT_OPR_INITIAL_STATES.values()) == [
        1e-2,  # % hi
        1e-2,  # % hp
        1e-2,  # % ht
        1e-2,  # % ha
        1e-2,  # % hc
        1e-6,  # % hlr
    ]


def test_default_bounds_parameters():
    # % Check default bounds opr parameters
    assert list(DEFAULT_BOUNDS_OPR_PARAMETERS.values()) == [
        (1e-6, 1e2),  # % ci
        (1e-6, 1e3),  # % cp
        (1e-6, 1e3),  # % ct
        (-50, 50),  # % kexc
        (1e-6, 0.999999),  # % aexc
        (1e-6, 1e3),  # % ca
        (1e-6, 1e3),  # % cc
        (1e-6, 4),  # % kb
        (1e-6, 1e3),  # % llr
        (1e-3, 50),  # % akw
        (1e-3, 1),  # % bkw
    ]

    # % Check default bounds opr states
    assert list(DEFAULT_BOUNDS_OPR_INITIAL_STATES.values()) == [
        (1e-6, 0.999999),  # % hi
        (1e-6, 0.999999),  # % hp
        (1e-6, 0.999999),  # % ht
        (1e-6, 0.999999),  # % ha
        (1e-6, 0.999999),  # % hc
        (1e-6, 1e3),  # % hlr
    ]


def test_read_input_data():
    # % Check input data format
    assert INPUT_DATA_FORMAT == ["tif", "nc"]

    # % Check ratio_pet_hourly
    assert np.array_equal(
        RATIO_PET_HOURLY,
        np.array(
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.035,
                0.062,
                0.079,
                0.097,
                0.11,
                0.117,
                0.117,
                0.11,
                0.097,
                0.079,
                0.062,
                0.035,
                0,
                0,
                0,
                0,
                0,
            ],
            dtype=np.float32,
        ),
    )


def test_event_seg():
    assert PEAK_QUANT == 0.995
    assert MAX_DURATION == 240


def test_dataset_name():
    assert DATASET_NAME == ["flwdir", "cance", "lez", "france"]
