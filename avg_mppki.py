#!/usr/bin/env python3
import os, re, argparse

def avg_mppki(dir):
    mppkis = []
    for fn in os.listdir(dir):
        if os.path.splitext(fn)[1] != '.out':
            continue
        fn = os.path.join(dir, fn)
        with open(fn) as f:
            txt = f.read()
        m = re.search(r'Final Score Run1_Conditional_MPPKI:\s*(.*)', txt)
        if m:
            mppkis.append(float(m.group(1)))
    return sum(mppkis) / max(len(mppkis), 1)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('run_dir')
    args = p.parse_args()
    print(avg_mppki(args.run_dir))
