import numpy
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import t
import math


def pitman_morgan_test(x, y, alternative, **kwargs):
    ratio = 1
    conf_level = 0.95

    alpha = 1 - conf_level
    n = len(x)
    df = n / 2

    #r = pearsonr(x, y)
    r = spearmanr(x, y)
    if isinstance(r, tuple) and r[1] <= 0.1:
        r = r[0]
    else:
        r = 0

    var1 = numpy.var(x, dtype=numpy.float64)
    var2 = numpy.var(y, dtype=numpy.float64)

    w = var1/var2


    stat_t = ((w - ratio) * numpy.sqrt(n - 2)) / numpy.sqrt(4 * (1 - pow(r,2) ) * w * ratio)

    ##############
    # quantile function "qt" in R = percent-point function "ppf" in scipy
    # "pt" in R = "cdf" in scipy

    if alternative == "two.sided":
        k = t.ppf(1-alpha, n-2)
        K = 1 + (2 * (1 - pow(r ,2 )) * pow(k, 2) ) / (n - 2)
        low = w * (K - numpy.sqrt( pow(K, 2) - 1))
        up = w * (K + numpy.sqrt( pow(K, 2) - 1))
        pval = 2 * (1 - t.cdf(abs(stat_t), df))
    elif alternative == "less":
        k = t.ppf(1-alpha, n-2)
        K = 1 + (2 * (1 - pow(r, 2 )) * pow(k, 2)) / (n - 2)
        low = 0
        up = w * (K + numpy.sqrt(pow(K, 2) - 1))
        pval = t.cdf(stat_t, df)
    elif alternative == "greater":
        k = t.ppf(1 - alpha, n - 2)
        K = 1 + (2 * (1 - pow(r, 2)) * pow(k, 2)) / (n - 2)
        low = w * (K - numpy.sqrt(pow(K, 2) - 1))
        up = math.inf
        pval = 1 - t.cdf(stat_t, df)

    estimate = {"variance x": var1, "variance y ": var2}
    cint = {"low ": low, "up ": up}
    method = "Paired Pitman-Morgan test"

    v = kwargs.get("verbose", False)
    if v:
        print(method)
        print("Statistic:" + str(stat_t))
        print("Confidence interval: "+str(cint))
        print("Estimate : "+str(estimate))
        print("p-value : "+str(pval))

    return pval
