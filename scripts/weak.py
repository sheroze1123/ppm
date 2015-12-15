import csv
import itertools
import sys
import matplotlib

import matplotlib.pyplot as plt

################################################################################
# parsing
################################################################################
FFTW_OMP = "FFTW OMP"
FULL_OMP = "Full OMP"
PARTIAL_OMP = "Partial OMP"
SERIAL = "Serial"

def flatten(ls):
    return [x for l in ls for x in l]

def typecast(d, types):
    return {k: types[k](v) for (k, v) in d.iteritems()}

def parse_file(filename):
    def name(s):
        names = {
            "ppm_omp": FFTW_OMP,
            "ppm_omp_full": FULL_OMP,
            "ppm_omp_partial": PARTIAL_OMP,
            "serial": SERIAL,
        }
        return names.get(s, s)

    types = {
        "name": name,
        "L": float,
        "N": int,
        "N_p": int,
        "dt": float,
        "T": int,
        "num_threads": int,
        "total_time": float,
        "ave_time": float,
    }

    with open(filename, "r") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=" ")
        return [typecast(row, types) for row in reader]

################################################################################
# querying
################################################################################
def sortuniq(iterator):
    return list(sorted(set(iterator)))

def matches(d, predicate):
    return all(d[k] == v for (k, v) in predicate.iteritems())

def domain(data, k):
    return sortuniq(d[k] for d in data)

def get(data, predicate, (k, ks), f):
    """
    predicate = {"N_p": 300, "N": 4}
    (k, ks) = ("num_threads", [1, 2, 3, 4, 5])
    f = "ave_time"
    get(data, predicate, ks, f)
    """
    res = {}
    for d in data:
        if matches(d, predicate) and d[k] in ks:
            if d[k] in res:
                print "predicate", predicate
                print "(k, ks)", (k, ks)
                print "f", f
                print "res", res
                print "d", d
                raise Exception("data not unique")
            else:
                res[d[k]] = d[f]
    try:
        return [res[k_] for k_ in ks]
    except KeyError, e:
        print "predicate", predicate
        print "(k, ks)", (k, ks)
        print "f", f
        print "res", res
        print "d", d
        raise e

################################################################################
# plotting
################################################################################
def linestyles(cnum=7, pnum=4):
    colors = ["b","g","r","c","m","y","k"][:cnum]
    patterns = ["-", "--", ":", "-."][:pnum]
    # for (p, c) in itertools.cycle(itertools.product(patterns, colors)):
    for (c, p) in itertools.cycle(itertools.product(colors, patterns)):
        yield c + p

def weak_scaling(data):
    programs = [FFTW_OMP, FULL_OMP, PARTIAL_OMP]
    num_threads = domain(data, "num_threads")
    N_ps = domain(data, "N_p")

    for N_p in N_ps:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        styles = linestyles(3, 1)
        for program in programs:
            predicate = {"name": program, "N_p":N_p}
            times = get(data, predicate, ("num_threads", num_threads), "ave_time")
            t0 = times[0]
            times = [t0 / t_ for t_ in times]
            ax.plot(num_threads, times, next(styles),
                    label="{}".format(program))

        plt.xlabel("number of threads, N_p={}".format(N_p))
        plt.ylabel("efficiency")
        plt.grid()
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, fancybox=True, shadow=True,
                        loc='upper center', bbox_to_anchor=(0.5, -0.1),
                        ncol = 2)
        fig.savefig("weak-N_p={}.pdf".format(N_p),
                    bbox_extra_artists=(lgd,),
                    bbox_inches='tight')
        plt.close()

def main(filenames):
    data = flatten(parse_file(filename) for filename in filenames)
    weak_scaling(data)

if __name__ == "__main__":
    matplotlib.rc('font', size=15)
    main(sys.argv[1:])
