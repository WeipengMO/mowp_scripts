import numpy as np
from scipy.stats import median_abs_deviation


def is_outlier(adata, metric: str, nmads: int):
    """Detect outliers in a metric.

    Parameters
    ----------
    adata
        Annotated data matrix.
    metric
        Metric to detect outliers in.
    nmads
        Number of median absolute deviations to use as threshold.
    """
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier