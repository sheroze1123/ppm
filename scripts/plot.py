import csv
import itertools
import sys

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

def time_vs_num_threads(data):
    programs = [FFTW_OMP, FULL_OMP, PARTIAL_OMP]
    num_threads = domain(data, "num_threads")

    for N in [128, 256, 512, 1024, 2048, 4096]:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        styles = linestyles(3, 1)
        N_p = 12800
        for program in programs:
            predicate = {"name": program, "N":N, "N_p":N_p}
            times = get(data, predicate, ("num_threads", num_threads), "ave_time")
            ax.plot(num_threads, times, next(styles),
                    label="{} (N = {})".format(program, N))

        plt.xlabel("number of threads, N_p={}".format(N_p))
        plt.ylabel("average time per time step (ms)")
        plt.grid()
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, fancybox=True, shadow=True,
                        loc='center left', bbox_to_anchor=(1, 0.5))
        fig.savefig("time_vs_num_threads-N={}.pdf".format(N),
                    bbox_extra_artists=(lgd,),
                    bbox_inches='tight')
        plt.close()

def speedup_vs_num_threads(data):
    programs = [FFTW_OMP, FULL_OMP, PARTIAL_OMP]
    num_threads = domain(data, "num_threads")

    for N in [128, 256, 512, 1024, 2048, 4096]:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        styles = linestyles(3, 1)
        N_p = 12800

        predicate = {"name": SERIAL, "num_threads": 1, "N_p": N_p, "N": N}
        serial_time, = get(data, predicate, ("num_threads", [1]), "ave_time")

        for program in programs:
            predicate = {"name": program, "N":N, "N_p":N_p}
            times = get(data, predicate, ("num_threads", num_threads), "ave_time")
            times = [serial_time / t_ for t_ in times]
            ax.plot(num_threads, times, next(styles),
                    label="{} (N = {})".format(program, N))

        plt.xlabel("number of threads, N_p={}".format(N_p))
        plt.ylabel("speedup compared to serial version")
        plt.grid()
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, fancybox=True, shadow=True,
                        loc='center left', bbox_to_anchor=(1, 0.5))
        fig.savefig("speedup_vs_num_threads-N={}.pdf".format(N),
                    bbox_extra_artists=(lgd,),
                    bbox_inches='tight')
        plt.close()

def time_vs_N(data):
    programs = [FFTW_OMP, FULL_OMP, PARTIAL_OMP]
    Ns = domain(data, "N")

    for t in [1, 2, 3, 4, 5, 6, 7, 8, 12, 16, 20, 24]:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        styles = linestyles(3, 1)
        N_p = 12800

        for program in programs:
            predicate = {"name": program, "num_threads":t, "N_p":N_p}
            times = get(data, predicate, ("N", Ns), "ave_time")
            ax.plot(Ns, times, next(styles),
                      label="{} ({} threads)".format(program, t))

        plt.xlabel("N, N_p={}".format(N_p))
        plt.ylabel("average time per time step (ms)")
        plt.grid()
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, fancybox=True, shadow=True,
                        loc='center left', bbox_to_anchor=(1, 0.5))
        fig.savefig("time_vs_N-t={}.pdf".format(t),
                    bbox_extra_artists=(lgd,),
                    bbox_inches='tight')
        plt.close()

def speedup_vs_N(data):
    programs = [FFTW_OMP, FULL_OMP, PARTIAL_OMP]
    Ns = domain(data, "N")

    for t in [1, 2, 3, 4, 5, 6, 7, 8, 12, 16, 20, 24]:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        styles = linestyles(3, 1)
        N_p = 12800

        predicate = {"name": SERIAL, "num_threads": 1, "N_p": N_p}
        serial_times = get(data, predicate, ("N", Ns), "ave_time")

        for program in programs:
            predicate = {"name": program, "num_threads":t, "N_p":N_p}
            times = get(data, predicate, ("N", Ns), "ave_time")
            times = [serial_time / t_ for (t_, serial_time) in zip(times, serial_times)]
            ax.plot(Ns, times, next(styles),
                      label="{} ({} threads)".format(program, t))

        plt.xlabel("N, N_p={}".format(N_p))
        plt.ylabel("speedup compared to serial implementation")
        plt.grid()
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, fancybox=True, shadow=True,
                        loc='center left', bbox_to_anchor=(1, 0.5))
        fig.savefig("speedup_vs_N-t={}.pdf".format(t),
                    bbox_extra_artists=(lgd,),
                    bbox_inches='tight')
        plt.close()

def time_vs_p(data):
    programs = [FFTW_OMP, FULL_OMP, PARTIAL_OMP]
    ps = domain(data, "N_p")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    styles = linestyles(3, 4)
    N = 128

    # for t in [1, 2, 3, 4, 5, 6, 7, 8, 12, 16, 20, 24]:
    for program in programs:
        for t in [1, 2, 3, 4]:
            predicate = {"name": program, "N":N, "num_threads": t}
            times = get(data, predicate, ("N_p", ps), "ave_time")
            ax.plot(ps, times, next(styles),
                    label="{} ({} threads)".format(program, t))

        plt.xlabel("number of particles, N={}".format(N))
        plt.ylabel("average time per time step (ms)")
        plt.grid()
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, fancybox=True, shadow=True,
                        loc='center left', bbox_to_anchor=(1, 0.5))
        fig.savefig("time_vs_p.pdf",
                    bbox_extra_artists=(lgd,),
                    bbox_inches='tight')
        plt.close()

def main(filenames):
    data = flatten(parse_file(filename) for filename in filenames)
    time_vs_num_threads(data)
    speedup_vs_num_threads(data)
    time_vs_N(data)
    speedup_vs_N(data)
    time_vs_p(data)

if __name__ == "__main__":
    main(sys.argv[1:])
