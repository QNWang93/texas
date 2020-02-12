import lc_ex
import sys

def main(filename):
#    filename = arg[0]
    f = open(filename, 'r')
    for line in f:
        para = line.split()
        print(para)
        para[3] = para[3][:4]
        lc_ex.main(para)

if __name__ == "__main__":
    main(sys.argv[1])
