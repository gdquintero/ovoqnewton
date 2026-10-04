import subprocess, re

# Sensitivity sweep for the Bard experiment of Section 4.2 over delta and sigmin,
# with gamma fixed at the value used in the paper. For each (delta, sigmin) and
# each method it reports f(x*) for o = 0..6, the drop factor between o = 3 and
# o = 4, whether the o = 4 run identifies exactly the four added outliers
# (indices 16-19, appended to the 15 clean observations), and the mean iteration
# and evaluation counts summed over o >= 1.

gamma = "5.d0"
deltas = ["1.0d-1", "1.0d-2", "1.0d-3", "1.0d-4"]
sigmins = ["1.0d-1", "1.0d-2", "1.0d-3"]
true_outliers = {16, 17, 18, 19}
methods = [("first-order", 1), ("second-order", 0)]

RE_BEST = re.compile(r"esta\s+\d+ &\s+([\d.E+-]+)")
RE_AVG = re.compile(r"avg\s+\d+\s+&\s+([\d.]+)\s+\+-\s+[\d.]+\s+&\s+([\d.]+)")


def run(delta, sigmin, o, flag):
    with open("param.txt", "w") as f:
        f.write("%s %s %s %d %d\n" % (delta, sigmin, gamma, o, flag))
    out = subprocess.run(["./bard"], capture_output=True, text=True).stdout
    m = RE_BEST.search(out)
    fovo = float(m.group(1)) if m else float("nan")
    a = RE_AVG.search(out)
    it_mean, ev_mean = (float(a.group(1)), float(a.group(2))) if a else (float("nan"),) * 2
    outl = None
    if o == 4:
        with open("../output/outliers_bard.txt") as fh:
            nums = [int(x) for x in fh.read().split()]
        outl = set(nums[1:1 + nums[0]])   # first entry is the count
    return fovo, it_mean, ev_mean, outl


hdr = "%-8s %-8s %-13s | %-9s %-9s %-9s | %-8s | %-6s | %s" % (
    "delta", "sigmin", "method", "f(o=2)", "f(o=3)", "f(o=4)",
    "drop", "detect", "sum(o>=1) it / fcnt")
print(hdr)
print("-" * len(hdr))

for delta in deltas:
    for sigmin in sigmins:
        for name, flag in methods:
            fovos, its, evs, ok = [], [], [], None
            for o in range(7):
                fovo, it_mean, ev_mean, outl = run(delta, sigmin, o, flag)
                fovos.append(fovo)
                its.append(it_mean)
                evs.append(ev_mean)
                if o == 4:
                    ok = (outl == true_outliers)
            drop = fovos[3] / fovos[4] if fovos[4] > 0 else float("inf")
            print("%-8s %-8s %-13s | %9.3e %9.3e %9.3e | %7.0fx | %-6s | %7.1f / %7.1f" % (
                delta, sigmin, name, fovos[2], fovos[3], fovos[4], drop,
                "yes" if ok else "NO",
                sum(its[1:]), sum(evs[1:])))

# Each ./bard call overwrites ../output/{solution,outliers}_bard.txt. Leave them
# holding the run the paper and the figure refer to (o = 4, first-order, with the
# parameters of Section 4.2) instead of whatever the last sweep entry produced.
run("1.0d-1", "1.0d-1", 4, 1)
print("\nRestored ../output/*_bard.txt to the o=4 first-order run of Section 4.2.")
