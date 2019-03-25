# from rpy2.robjects.packages import importr
from rpy2 import robjects
import scipy.stats as scstat
robjects.r('library(Barnard)')
robjects.r('''barnard_test <- function(a, b, c, d) {invisible(capture.output(y<-barnardw.test(a,b,c,d)$p.value)); y}''')
barnard_test = robjects.globalenv['barnard_test']

barnard_memoize = {}
fisher_memoize = {}


def get_barnard(a, b, c, d):
    if (a, b, c, d) not in barnard_memoize:
        barnard_memoize[(a, b, c, d)] = barnard_test(a, b, c, d)[0]
    return barnard_memoize[(a, b, c, d)]


def get_fisher(a, b, c, d):
    if (a, b, c, d) not in fisher_memoize:
        fisher_memoize[(a, b, c, d)] = scstat.fisher_exact([[a, b], [c, d]])[-1]
    return fisher_memoize[(a, b, c, d)]

# def barnard_test(a, b, c, d):
# 	p_values = barnard.barnardw_test(a,b,c,d, dp = 0.001)[5]
# 	one_sided, two_sided = p_values
# 	return round(one_sided, 4)