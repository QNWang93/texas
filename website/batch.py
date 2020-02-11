import lc_ex
import sys

if __name__ == "__main__":
    filename = sys.argv[1]
    f = open(filename, 'r')
    for line in f:
        para = line.split()
        print(para)
        para[3] = para[3][:4]
        lc_ex.main(para)
