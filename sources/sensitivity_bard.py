import subprocess, re, os

# Sensitivity sweep for Bard over delta and sigmin (gamma fixed, first-order).
# For each (delta, sigmin) it reports f(x*) for o=0..6 and whether the o=4 run
# identifies exactly the injected outliers {3,7,11,15}.

gamma = "2.d0"
deltas = ["1.0d-1", "1.0d-2", "1.0d-3", "1.0d-4"]
sigmins = ["1.0d-1", "1.0d-2", "1.0d-3"]
true_outliers = {3, 7, 11, 15}

def run(delta, sigmin, o):
    with open("param.txt", "w") as f:
        f.write("%s %s %s %d 1\n" % (delta, sigmin, gamma, o))
    out = subprocess.run(["./bard"], capture_output=True, text=True).stdout
    m = re.search(r"esta\s+\d+ &\s+([\d.E+-]+)", out)
    fovo = float(m.group(1)) if m else float("nan")
    outl = None
    if o == 4:
        with open("../output/outliers_bard.txt") as fh:
            nums = [int(x) for x in fh.read().split()]
        outl = set(nums[1:1+nums[0]])   # first entry is the count
    return fovo, outl

print("%-8s %-8s | %s | drop(o3->o4) | o4 outliers OK" % ("delta", "sigmin", "  ".join("o=%d"%o for o in range(7))))
print("-"*110)
for delta in deltas:
    for sigmin in sigmins:
        fovos = []
        o4_ok = None
        for o in range(7):
            fovo, outl = run(delta, sigmin, o)
            fovos.append(fovo)
            if o == 4:
                o4_ok = (outl == true_outliers)
        drop = fovos[3]/fovos[4] if fovos[4] > 0 else float("inf")
        row = "  ".join("%.1e" % v for v in fovos)
        print("%-8s %-8s | %s | %10.1fx | %s%s" % (
            delta, sigmin, row, drop,
            "YES" if o4_ok else "no ",
            "" if o4_ok else "  <-- fails"))
