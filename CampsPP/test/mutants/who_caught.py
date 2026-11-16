"""Which tests caught each mutant: failing checks grouped by test file and testset.

usage: python3 -I who_caught.py <results dir> ID...

Each failing check is tagged "old" when its testset runs on the old CliffordMPS stack or
compares against the old code (testset names with "old CAMPS", ", CAMPS:", "old code",
"old vs new" or "reproduce the old"). Those checks disappear with the Phase 6 deletions, so a
mutant whose failures are all "old" would survive them. The tag goes by testset name: a testset
that compares against both dense and old code (e.g. "vs dense and old code") counts as old,
which overstates the old share.
"""
import os
import re
import sys
from collections import Counter, defaultdict

out, ids = sys.argv[1], sys.argv[2:]
head = re.compile(r"^(?P<ts>.*): (?P<kind>Test Failed|Error During Test) at (?P<path>\S+):(?P<line>\d+)")
OLD = re.compile(r"old CAMPS|, CAMPS$|old code|old vs new|reproduce the old")

for mid in ids:
    path = os.path.join(out, f"{mid}.log")
    if not os.path.exists(path):
        print(f"=== {mid}: no log")
        continue
    groups, sites = Counter(), defaultdict(set)
    for l in open(path, encoding="utf-8", errors="replace"):
        m = head.match(l)
        if not m:
            continue
        rel = m["path"].split("/test/", 1)[-1]
        groups[(rel, m["ts"].strip(), m["kind"][0])] += 1
        sites[(rel, m["ts"].strip(), m["kind"][0])].add(int(m["line"]))
    tot = sum(groups.values())
    old = sum(n for (_, ts, _), n in groups.items() if OLD.search(ts))
    print(f"=== {mid}: {tot} failing checks, {tot - old} outside old-code tests"
          + ("   <-- all old: would survive Phase 6" if tot and old == tot else ""))
    for (rel, ts, kind), n in sorted(groups.items(), key=lambda kv: -kv[1]):
        tag = " [old]" if OLD.search(ts) else ""
        print(f"  {kind} {n:4d}  {rel} :: {ts}  (lines {sorted(sites[(rel, ts, kind)])[:6]}){tag}")
