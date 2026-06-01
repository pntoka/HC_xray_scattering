"""Default parameter bounds and scaling for iobs_ngc fitting."""

# (min_val, max_val) for min-max scaling
DEFAULT_PARAM_BOUNDS = {
    'mu':     (1.2, 8.0),
    'beta':   (0.01, 8.0),
    'a3':     (2.7, 4.0),
    'da3':    (0.01, 1.0),
    'sig3':   (0.01, 2.0),
    'eta':    (0.6, 1.0),
    'alpha':  (0.01, 1.5),
    'sig1':   (0.0, 0.6),
    'lcc':    (1.2, 1.5),
    'q':      (0.0, 0.75),
    'g':      (-1.0, 1.0),
    'u3':     (0.0, 0.4),
    'const1': (-500.0, 750.0),
    'const2': (-2.0, 5.0),
    'k':      (1.0, 1000.0),
}

# Default sequential fitting stages
DEFAULT_STAGE_DEFINITIONS = [
    {
        'name': 'normalization',
        'param_names': ['q', 'k', 'const1', 'const2', 'g'],
    },
    {
        'name': 'interlayer',
        'param_names': ['mu', 'beta', 'a3', 'da3', 'sig3', 'eta'],
    },
    {
        'name': 'intralayer',
        'param_names': ['alpha', 'sig1'],
    },
    {
        'name': 'lcc',
        'param_names': ['lcc'],
    },
    {
        'name': 'all_parameters',
        'param_names': list(DEFAULT_PARAM_BOUNDS.keys()),
    },
]