import sys

if __name__ == '__main__':
    lines = list()
    for line in sys.stdin: lines.append(line.strip().split(',')[0:2])
    lines.sort()
    for l in lines:
        f, km = l
        sys.stdout.write("{},{}\n".format(f,km))
