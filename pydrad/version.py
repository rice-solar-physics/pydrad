# NOTE: First try _dev.scm_version if it exists and setuptools_scm is installed
# This file is not included in wheels/tarballs, so otherwise it will
# fall back on the generated _version module.
try:
    try:
        from pydrad._dev.scm_version import version
    except ImportError:
        from pydrad._version import version
except Exception:
    import warnings

    warnings.warn(
        f'could not determine {__name__.split(".")[0]} package version; this indicates a broken installation'
    )
    del warnings

    version = '0.0.0'
