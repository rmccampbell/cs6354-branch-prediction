#!/usr/bin/env python3
import sys, os, re

def avg_mppki(dir):
    mppkis = []
    for fn in os.listdir(dir):
        fn = os.path.join(dir, fn)
        with open(fn) as f:
            txt = f.read()
        m = re.search(r'Final Score Run1_Conditional_MPPKI:\s*(.*)', txt)
        mppkis.append(float(m.group(1)))
    return sum(mppkis) / len(mppkis)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit('usage: avg_mppki.py <run_dir>')
    print(avg_mppki(sys.argv[1]))
