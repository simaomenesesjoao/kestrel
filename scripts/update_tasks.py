import sys

key = sys.argv[2].split("=")[0]
value = sys.argv[2].split("=")[1]
filename = "/home/simao/Sync/article_magnetic/program/kestrel/scripts/tasks.txt"

file1 = open(filename,'r')
f = file1.readlines()
file1.close()

print_string = ""
for line in f:
    splitted = line[:-1].split(" ")

    Lx          = splitted[0][3:]
    Ly          = splitted[1][3:]
    seed        = splitted[2][5:]
    w           = splitted[3][2:]
    nx          = splitted[4][3:]
    ny          = splitted[5][3:]
    B           = splitted[6][2:]
    energies    = splitted[7][3:-1]
    NR          = splitted[8][3:]
    ND          = splitted[9][3:]
    N_max       = int(splitted[10][2:])
    N_slice     = int(splitted[11][3:])
    N_current   = int(splitted[12][3:])
    status      = splitted[13][7:]


    if(sys.argv[1] == "B={0}".format(B)):
        # print("found")
        if(key == "status"):
            status = value
        elif(key =="Na"):
            N_current = value


    print_string +=  "Lx={0} Ly={1} seed={2} W={3} nx={4} ny={5} B={6} E=[{7}] NR={8} ND={9} N={10} dN={11} Na={12} status={13}\n".format(Lx, Ly, seed, w, nx, ny, B, energies, NR, ND, N_max, N_slice, N_current, status)

file2 = open(filename, "w")
file2.write(print_string)
file2.close()
