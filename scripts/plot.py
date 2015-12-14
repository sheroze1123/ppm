import csv
import itertools
import sys

import matplotlib.pyplot as plt

################################################################################
# parsing
################################################################################
def flatten(ls):
    return [x for l in ls for x in l]

def typecast(d, types):
    return {k: types[k](v) for (k, v) in d.iteritems()}

def parse_file(filename):
    types = {
        "name": str,
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
    return [res[k] for k in ks]

################################################################################
# plotting
################################################################################
def linestyles():
    colors = "bgrcmyk"
    patterns = ["-", "--", "-.", ":"]
    for (p, c) in itertools.cycle(itertools.product(patterns, colors)):
        yield c + p

def time_vs_num_threads(data):
    programs = domain(data, "name")
    num_threads = domain(data, "num_threads")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    styles = linestyles()
    N_p = 6400

    for program in programs:
        for N in [128, 256, 512, 1024, 2048]:
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
    fig.savefig("time_vs_num_threads.pdf",
                bbox_extra_artists=(lgd,),
                bbox_inches='tight')

def time_vs_N(data):
    programs = domain(data, "name")
    Ns = domain(data, "N")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    styles = linestyles()
    N_p = 6400

    for program in programs:
        for t in [1, 2, 4, 12, 24]:
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
    fig.savefig("time_vs_N.pdf",
                bbox_extra_artists=(lgd,),
                bbox_inches='tight')

def time_vs_p(data):
    programs = domain(data, "name")
    ps = domain(data, "N_p")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    styles = linestyles()
    N = 2048

    for program in programs:
        for t in [1, 2, 4, 12, 24]:
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

def main(filenames):
    data = flatten(parse_file(filename) for filename in filenames)
    time_vs_num_threads(data)
    time_vs_N(data)
    time_vs_p(data)

if __name__ == "__main__":
    main(sys.argv[1:])
